{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE MagicHash #-}
{-# LANGUAGE UnboxedTuples #-}
{-# LANGUAGE BangPatterns #-}

-- |
-- Module      : Numeric.LinearAlgebra.Massiv.Eigen.Symmetric
-- Copyright   : (c) Nadia Chambers 2026
-- License     : BSD-3-Clause
-- Maintainer  : nadia.chambers@iohk.io
-- Stability   : experimental
--
-- Eigenvalue algorithms specialised to real symmetric matrices, following
-- Golub & Van Loan, /Matrix Computations/, 4th edition (GVL4), Chapter 8,
-- pp. 449--512.
module Numeric.LinearAlgebra.Massiv.Eigen.Symmetric
  ( tridiagonalize
  , symmetricEigen
  , symmetricEigenP
  , symmetricEigenPPar
  , jacobiEigen
  ) where

import qualified Data.Massiv.Array as M
import Data.Massiv.Array (Ix2(..), Ix1, unwrapByteArray, unwrapByteArrayOffset, unwrapMutableByteArray, unwrapMutableByteArrayOffset)
import GHC.TypeNats (KnownNat)
import Control.Monad (when, forM_)
import Control.Monad.ST (ST, stToIO)
import Control.Concurrent (forkIO, newEmptyMVar, putMVar, takeMVar)
import System.IO.Unsafe (unsafePerformIO)
import GHC.Exts
import GHC.ST (ST(..))
import Data.Primitive.ByteArray (MutableByteArray(..), newByteArray, unsafeFreezeByteArray)

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.Orthogonal.Givens (givensRotation)
import Numeric.LinearAlgebra.Massiv.Internal.Kernel
  ( rawMutApplyGivensColumns
  , rawMutSumSqColumn
  , rawMutSymMatvecSub
  , rawMutSymRank2Update
  , rawMutTridiagQAccum
  , rawGemmKernel
  )

-- | Reduce a symmetric matrix to tridiagonal form (GVL4 Algorithm 8.3.1).
tridiagonalize :: forall n r e. (KnownNat n, M.Manifest r e, Floating e, Ord e)
               => Matrix n n r e
               -> (Matrix n n r e, Vector n r e, Vector n r e)
tridiagonalize a =
  let nn = dimVal @n

      -- Phase 1: In-place tridiagonalisation via symmetric rank-2 updates.
      (betaList, tArr) = M.withMArrayST (unMatrix a) $ \mt -> do
        betas <- mapM (tridiagStep mt nn) [0..nn-3]
        pure betas

      -- Phase 2: Accumulate Q from stored Householder vectors.
      qMat = createMatrix @n @n @r $ \mq -> do
        forM_ [0..nn-1] $ \i -> forM_ [0..nn-1] $ \j ->
          M.write_ mq (i :. j) (if i == j then 1 else 0)
        -- Forward accumulation: Q <- Q · H_k for k = 0..n-3
        forM_ (zip [0..] betaList) $ \(k, beta_k) ->
          when (beta_k /= 0) $
            forM_ [0..nn-1] $ \i -> do
              qik1 <- M.readM mq (i :. (k+1))
              rest <- sumQV mq tArr i (k+1) nn k
              let wi = beta_k * (qik1 + rest)
              M.write_ mq (i :. (k+1)) (qik1 - wi)
              forM_ [k+2..nn-1] $ \l -> do
                let vl = M.index' tArr (l :. k)
                qil <- M.readM mq (i :. l)
                M.write_ mq (i :. l) (qil - wi * vl)

      diag_ = makeVector @n @r $ \i -> M.index' tArr (i :. i)
      subdiag = makeVector @n @r $ \i ->
        if i < nn - 1 then M.index' tArr ((i+1) :. i) else 0

  in (qMat, diag_, subdiag)

-- | One step of Householder tridiagonalisation.
tridiagStep :: (M.Manifest r e, Floating e, Ord e)
            => M.MArray s r Ix2 e -> Int -> Int -> ST s e
tridiagStep mt nn k = do
  x0 <- M.readM mt ((k+1) :. k)
  sigma <- sumSqBelow mt (k+1) nn k
  if sigma == 0 && x0 >= 0
    then pure 0
    else do
      let mu = sqrt (x0 * x0 + sigma)
          v0 = if x0 <= 0 then x0 - mu else -sigma / (x0 + mu)
          beta = 2 * v0 * v0 / (sigma + v0 * v0)
      -- Build v as a list: v(k+1)=1, v(i)=T(i,k)/v0 for i>k+1
      vList <- mapM (\i -> do
        tik <- M.readM mt (i :. k)
        pure (tik / v0)
        ) [k+2..nn-1]
      let fullV = 1 : vList  -- indices k+1, k+2, ..., n-1
      -- p = beta * T * v (rows k+1..n-1)
      pList <- mapM (\i -> do
        s <- dotTV mt i fullV (k+1) nn
        pure (beta * s)
        ) [k+1..nn-1]
      let ptv = sum $ zipWith (*) pList fullV
          alpha_ = beta * ptv / 2
          wList = zipWith (\pi_ vi -> pi_ - alpha_ * vi) pList fullV
      -- Symmetric rank-2 update: T(i,j) -= v(i)*w(j) + w(i)*v(j)
      forM_ (zip3 [k+1..nn-1] fullV wList) $ \(i, vi, wi) ->
        forM_ (zip3 [k+1..nn-1] fullV wList) $ \(j, vj, wj) -> do
          tij <- M.readM mt (i :. j)
          M.write_ mt (i :. j) (tij - vi * wj - wi * vj)
      -- Store Householder vector in below-subdiagonal of column k
      forM_ (zip [k+2..nn-1] vList) $ \(i, vi) ->
        M.write_ mt (i :. k) vi
      -- Set subdiagonal
      M.write_ mt ((k+1) :. k) mu
      M.write_ mt (k :. (k+1)) mu
      pure beta

-- Helpers for tridiagonalize
sumSqBelow :: (M.Manifest r e, Num e) => M.MArray s r Ix2 e -> Int -> Int -> Int -> ST s e
sumSqBelow mt start end col = go (start + 1) 0
  where go i !acc | i >= end = pure acc
                  | otherwise = do v <- M.readM mt (i :. col); go (i+1) (acc + v*v)

dotTV :: (M.Manifest r e, Num e) => M.MArray s r Ix2 e -> Int -> [e] -> Int -> Int -> ST s e
dotTV mt i vList start end = go start vList 0
  where go _ [] !acc = pure acc
        go j (v:vs) !acc | j >= end = pure acc
                         | otherwise = do t <- M.readM mt (i :. j); go (j+1) vs (acc + t*v)

sumQV :: (M.Manifest r1 e, M.Manifest r2 e, Num e)
      => M.MArray s r1 Ix2 e -> M.Array r2 Ix2 e -> Int -> Int -> Int -> Int -> ST s e
sumQV mq tArr row start end col = go (start + 1) 0
  where go l !acc | l >= end = pure acc
                  | otherwise = do
                      q <- M.readM mq (row :. l)
                      let v = M.index' tArr (l :. col)
                      go (l+1) (acc + q*v)

-- | Raw-primop tridiagonalisation specialised for @P Double@.
-- Two-phase: (1) in-place Householder via raw ByteArray# kernels,
-- (2) Q accumulation using rawMutTridiagQAccum.
tridiagonalizeP :: forall n. KnownNat n
                => Matrix n n M.P Double
                -> (Matrix n n M.P Double, Vector n M.P Double, Vector n M.P Double)
tridiagonalizeP a =
  let nn = dimVal @n

      -- Phase 1: In-place Householder tridiagonalisation
      (betaList, tArr) = M.withMArrayST (unMatrix a) $ \mt -> do
        let !mbaT = unwrapMutableByteArray mt
            !offT = unwrapMutableByteArrayOffset mt
        -- Allocate reusable temporary vectors (n Doubles = n*8 bytes)
        mbaV <- newByteArray (nn * 8)
        mbaP <- newByteArray (nn * 8)
        mbaW <- newByteArray (nn * 8)
        betas <- mapM (\k -> tridiagStepP mbaT offT nn mbaV mbaP mbaW k) [0..nn-3]
        pure betas

      -- Get underlying ByteArray from frozen T for Q accumulation
      !tBA  = unwrapByteArray tArr
      !tOff = unwrapByteArrayOffset tArr

      -- Phase 2: Accumulate Q using raw primops
      qMat = createMatrix @n @n @M.P $ \mq -> do
        let !mbaQ = unwrapMutableByteArray mq
            !offQ = unwrapMutableByteArrayOffset mq
        -- Set Q = I
        forM_ [0..nn-1] $ \i -> forM_ [0..nn-1] $ \j ->
          writeRawD mbaQ offQ (i*nn+j) (if i == j then 1 else 0)
        -- Forward accumulation: Q <- Q · H_k for k = 0..n-3
        forM_ (zip [0..] betaList) $ \(k, beta_k) ->
          when (beta_k /= 0) $
            forM_ [0..nn-1] $ \row ->
              rawMutTridiagQAccum mbaQ offQ nn tBA tOff nn beta_k (k+1) k nn row

      -- Read diagonal and subdiagonal from frozen T
      diag_   = makeVector @n @M.P $ \i -> M.index' tArr (i :. i)
      subdiag = makeVector @n @M.P $ \i ->
        if i < nn - 1 then M.index' tArr ((i+1) :. i) else 0

  in (qMat, diag_, subdiag)
{-# NOINLINE tridiagonalizeP #-}

-- | One step of raw-primop Householder tridiagonalisation.
tridiagStepP :: MutableByteArray s -> Int -> Int
             -> MutableByteArray s -> MutableByteArray s -> MutableByteArray s
             -> Int -> ST s Double
tridiagStepP mbaT offT nn mbaV mbaP mbaW k = do
  -- 1. Read x0 = T[k+1,k]
  x0 <- readRawD mbaT offT ((k+1)*nn + k)
  -- 2. Compute sigma = Σ T[i,k]^2 for i=k+2..nn-1
  sigma <- rawMutSumSqColumn mbaT offT nn (k+2) nn k
  if sigma == 0 && x0 >= 0
    then pure 0
    else do
      let mu = sqrt (x0 * x0 + sigma)
          v0 = if x0 <= 0 then x0 - mu else -sigma / (x0 + mu)
          beta = 2 * v0 * v0 / (sigma + v0 * v0)
          subSize = nn - k - 1

      -- 3. Build v in mbaV: v[0]=1, v[i]=T[k+1+i,k]/v0
      writeRawD mbaV 0 0 1.0
      forM_ [1..subSize-1] $ \i -> do
        tik <- readRawD mbaT offT ((k+1+i)*nn + k)
        writeRawD mbaV 0 i (tik / v0)

      -- 4. p = beta * T_sub * v
      rawMutSymMatvecSub mbaT offT nn mbaV 0 mbaP 0 (k+1) nn
      forM_ [0..subSize-1] $ \i -> do
        pi_ <- readRawD mbaP 0 i
        writeRawD mbaP 0 i (beta * pi_)

      -- 5. Dot product p^T v
      ptv <- mutDotVec mbaP 0 mbaV 0 subSize
      let alpha_ = beta * ptv / 2

      -- 6. w = p - alpha*v
      forM_ [0..subSize-1] $ \i -> do
        pi_ <- readRawD mbaP 0 i
        vi  <- readRawD mbaV 0 i
        writeRawD mbaW 0 i (pi_ - alpha_ * vi)

      -- 7. Rank-2 update: T -= vw^T + wv^T
      rawMutSymRank2Update mbaT offT nn mbaV 0 mbaW 0 (k+1) nn

      -- 8. Store Householder vector in column k subdiagonal
      forM_ [1..subSize-1] $ \i -> do
        vi <- readRawD mbaV 0 i
        writeRawD mbaT offT ((k+1+i)*nn + k) vi

      -- 9. Set subdiagonal element
      writeRawD mbaT offT ((k+1)*nn + k) mu
      writeRawD mbaT offT (k*nn + (k+1)) mu

      pure beta

-- | Dot product of two mutable vectors.
mutDotVec :: MutableByteArray s -> Int -> MutableByteArray s -> Int -> Int -> ST s Double
mutDotVec mbaA offA mbaB offB n = go 0 0
  where
    go !i !acc
      | i >= n = pure acc
      | otherwise = do
          ai <- readRawD mbaA offA i
          bi <- readRawD mbaB offB i
          go (i+1) (acc + ai * bi)

-- | Symmetric eigenvalue decomposition (GVL4 Algorithm 8.3.3).
symmetricEigen :: forall n r e. (KnownNat n, M.Manifest r e, Floating e, Ord e)
               => Matrix n n r e -> Int -> e -> (Vector n r e, Matrix n n r e)
symmetricEigen a maxIter tol =
  let nn = dimVal @n
      (q0, diag_, subdiag) = tridiagonalize a
      (dArr, qArr) = M.withMArrayST (unMatrix q0) $ \mq -> do
        md <- M.thawS (unVector diag_)
        msd <- M.thawS (unVector subdiag)
        tridiagQRLoop md msd mq nn maxIter tol
        dFrozen <- M.freezeS md
        pure (MkVector dFrozen)
  in (dArr, MkMatrix qArr)

-- | In-place QR iteration on tridiagonal (d, sd) with mutable Q.
-- Uses both top and bottom deflation to shrink the active range [lo..hi],
-- effectively achieving divide-and-conquer behaviour.
tridiagQRLoop :: (M.Manifest r e, Floating e, Ord e)
              => M.MArray s r Ix1 e -> M.MArray s r Ix1 e -> M.MArray s r Ix2 e
              -> Int -> Int -> e -> ST s ()
tridiagQRLoop md msd mq nn maxIter tol = go 0 0 (nn - 1)
  where
    go !iter !lo !hi
      | iter >= maxIter = pure ()
      | lo >= hi = pure ()
      | otherwise = do
          -- Bottom deflation
          sdhi <- M.readM msd (hi - 1)
          dhi1 <- M.readM md (hi - 1)
          dhi  <- M.readM md hi
          if abs sdhi <= tol * (abs dhi1 + abs dhi)
            then do
              M.write_ msd (hi - 1) 0
              go iter lo (hi - 1)
            else do
              -- Top deflation
              sdlo <- M.readM msd lo
              dlo  <- M.readM md lo
              dlo1 <- M.readM md (lo + 1)
              if abs sdlo <= tol * (abs dlo + abs dlo1)
                then do
                  M.write_ msd lo 0
                  go iter (lo + 1) hi
                else do
                  -- Interior deflation: find split point
                  split <- findSplit md msd lo hi tol
                  case split of
                    Just q -> do
                      -- Split into two subproblems [lo..q] and [q+1..hi]
                      M.write_ msd q 0
                      go iter lo q
                      go iter (q + 1) hi
                    Nothing -> do
                      -- No split found: apply QR step on [lo..hi]
                      let sp1 = sdhi
                          delta = (dhi1 - dhi) / 2
                          sgn = if delta >= 0 then 1 else -1
                          shift = dhi - sp1*sp1 / (delta + sgn * sqrt (delta*delta + sp1*sp1))
                      implicitQRStepInPlace md msd mq nn shift lo hi
                      go (iter + 1) lo hi

-- | Find an interior split point where the subdiagonal is negligible.
findSplit :: (M.Manifest r e, Floating e, Ord e)
          => M.MArray s r Ix1 e -> M.MArray s r Ix1 e -> Int -> Int -> e -> ST s (Maybe Int)
findSplit md msd lo hi tol = scan (lo + 1)
  where
    scan q
      | q >= hi - 1 = pure Nothing
      | otherwise = do
          sdq <- M.readM msd q
          dq  <- M.readM md q
          dq1 <- M.readM md (q + 1)
          if abs sdq <= tol * (abs dq + abs dq1)
            then pure (Just q)
            else scan (q + 1)

-- | One implicit symmetric QR step via bulge-chasing Givens rotations.
-- Operates on the active sub-range [lo..hi] of the tridiagonal.
implicitQRStepInPlace :: (M.Manifest r e, Floating e, Ord e)
                      => M.MArray s r Ix1 e -> M.MArray s r Ix1 e -> M.MArray s r Ix2 e
                      -> Int -> e -> Int -> Int -> ST s ()
implicitQRStepInPlace md msd mq nn shift lo hi = do
  dlo <- M.readM md lo
  sdlo <- M.readM msd lo
  chase lo (dlo - shift) sdlo
  where
    chase k x z = do
      let (c, s) = givensRotation x z
      when (k > lo) $
        M.write_ msd (k-1) (c * x - s * z)
      dk  <- M.readM md k
      ek  <- M.readM msd k
      dk1 <- M.readM md (k+1)
      M.write_ md k     (c*c*dk - 2*c*s*ek + s*s*dk1)
      M.write_ md (k+1) (s*s*dk + 2*c*s*ek + c*c*dk1)
      M.write_ msd k    (c*s*(dk - dk1) + (c*c - s*s)*ek)
      applyGivensRightQ mq c s k (k+1) nn
      if k + 1 < hi
        then do
          ek1 <- M.readM msd (k+1)
          let z' = -s * ek1
          M.write_ msd (k+1) (c * ek1)
          ek_new <- M.readM msd k
          chase (k+1) ek_new z'
        else pure ()

-- | Apply Givens rotation from the right to Q: Q <- Q · G(ci, ck)
-- For P Double, uses raw ByteArray# primops; generic fallback otherwise.
applyGivensRightQ :: (M.Manifest r e, Num e)
                  => M.MArray s r Ix2 e -> e -> e -> Int -> Int -> Int -> ST s ()
applyGivensRightQ mq c s ci ck nn =
  forM_ [0..nn-1] $ \row -> do
    qrc <- M.readM mq (row :. ci)
    qrk <- M.readM mq (row :. ck)
    M.write_ mq (row :. ci) (c * qrc - s * qrk)
    M.write_ mq (row :. ck) (s * qrc + c * qrk)

-- | Specialised symmetric eigenvalue decomposition for @P Double@.
-- Uses raw ByteArray# primops for the entire QR iteration, including
-- diagonal/subdiagonal reads and writes plus Givens rotation on Q.
symmetricEigenP :: forall n. KnownNat n
                => Matrix n n M.P Double -> Int -> Double -> (Vector n M.P Double, Matrix n n M.P Double)
symmetricEigenP a maxIter tol =
  let nn = dimVal @n
      (q0, diag_, subdiag) = tridiagonalizeP a
      (dArr, qArr) = M.withMArrayST (unMatrix q0) $ \mq -> do
        md <- M.thawS (unVector diag_)
        msd <- M.thawS (unVector subdiag)
        let !mbaD  = unwrapMutableByteArray md
            !offD  = unwrapMutableByteArrayOffset md
            !mbaSD = unwrapMutableByteArray msd
            !offSD = unwrapMutableByteArrayOffset msd
            !mbaQ  = unwrapMutableByteArray mq
            !offQ  = unwrapMutableByteArrayOffset mq
        rawTridiagQRLoop mbaD offD mbaSD offSD mbaQ offQ nn maxIter tol
        dFrozen <- M.freezeS md
        pure (MkVector dFrozen)
  in (dArr, MkMatrix qArr)
{-# NOINLINE symmetricEigenP #-}

-- | Parallel specialised symmetric eigenvalue decomposition for @P Double@.
-- Uses raw-primop tridiagonalisation and forks independent sub-problems
-- when the QR loop finds a split point.
symmetricEigenPPar :: forall n. KnownNat n
                   => Matrix n n M.P Double -> Int -> Double
                   -> (Vector n M.P Double, Matrix n n M.P Double)
symmetricEigenPPar a maxIter tol = unsafePerformIO $ do
  let nn = dimVal @n
      (q0, diag_, subdiag) = tridiagonalizeP a
  -- Thaw into IO (s = RealWorld) for parallel QR iteration
  mq  <- M.thawS (unMatrix q0)
  md  <- M.thawS (unVector diag_)
  msd <- M.thawS (unVector subdiag)
  let !mbaD  = unwrapMutableByteArray md
      !offD  = unwrapMutableByteArrayOffset md
      !mbaSD = unwrapMutableByteArray msd
      !offSD = unwrapMutableByteArrayOffset msd
      !mbaQ  = unwrapMutableByteArray mq
      !offQ  = unwrapMutableByteArrayOffset mq
  rawTridiagQRLoopPar mbaD offD mbaSD offSD mbaQ offQ nn maxIter tol
  dFrozen <- M.freezeS md
  qFrozen <- M.freezeS mq
  pure (MkVector dFrozen, MkMatrix qFrozen)
{-# NOINLINE symmetricEigenPPar #-}

-- | Parallel QR loop: forks independent sub-problems when a split is found.
-- Operates in IO to enable forkIO for non-overlapping sub-problem ranges.
rawTridiagQRLoopPar :: MutableByteArray RealWorld -> Int
                    -> MutableByteArray RealWorld -> Int
                    -> MutableByteArray RealWorld -> Int
                    -> Int -> Int -> Double -> IO ()
rawTridiagQRLoopPar mbaD offD mbaSD offSD mbaQ offQ nn maxIter tol = go 0 0 (nn - 1)
  where
    rd mba off i = stToIO (readRawD mba off i)
    wr mba off i v = stToIO (writeRawD mba off i v)

    go !iter !lo !hi
      | iter >= maxIter = pure ()
      | lo >= hi = pure ()
      | otherwise = do
          sdhi <- rd mbaSD offSD (hi - 1)
          dhi1 <- rd mbaD offD (hi - 1)
          dhi  <- rd mbaD offD hi
          if abs sdhi <= tol * (abs dhi1 + abs dhi)
            then do wr mbaSD offSD (hi - 1) 0; go iter lo (hi - 1)
            else do
              sdlo <- rd mbaSD offSD lo
              dlo  <- rd mbaD offD lo
              dlo1 <- rd mbaD offD (lo + 1)
              if abs sdlo <= tol * (abs dlo + abs dlo1)
                then do wr mbaSD offSD lo 0; go iter (lo + 1) hi
                else do
                  split <- stToIO $ rawFindSplit mbaD offD mbaSD offSD lo hi tol
                  case split of
                    Just q -> do
                      wr mbaSD offSD q 0
                      done <- newEmptyMVar
                      _ <- forkIO $ do
                        go iter lo q
                        putMVar done ()
                      go iter (q + 1) hi
                      takeMVar done
                    Nothing -> do
                      let sp1 = sdhi
                          delta = (dhi1 - dhi) / 2
                          sgn = if delta >= 0 then 1 else -1
                          shift = dhi - sp1*sp1 / (delta + sgn * sqrt (delta*delta + sp1*sp1))
                      stToIO $ rawImplicitQRStep mbaD offD mbaSD offSD mbaQ offQ nn shift lo hi
                      go (iter + 1) lo hi

-- | Read a Double from a raw MutableByteArray at element index.
readRawD :: MutableByteArray s -> Int -> Int -> ST s Double
readRawD (MutableByteArray mba) (I# off) (I# i) = ST $ \s ->
  case readDoubleArray# mba (off +# i) s of
    (# s', v #) -> (# s', D# v #)
{-# INLINE readRawD #-}

-- | Write a Double to a raw MutableByteArray at element index.
writeRawD :: MutableByteArray s -> Int -> Int -> Double -> ST s ()
writeRawD (MutableByteArray mba) (I# off) (I# i) (D# v) = ST $ \s ->
  case writeDoubleArray# mba (off +# i) v s of
    s' -> (# s', () #)
{-# INLINE writeRawD #-}

-- | Read an Int from a raw MutableByteArray at element index.
readRawI :: MutableByteArray s -> Int -> Int -> ST s Int
readRawI (MutableByteArray mba) (I# off) (I# i) = ST $ \s ->
  case readIntArray# mba (off +# i) s of
    (# s', v #) -> (# s', I# v #)
{-# INLINE readRawI #-}

-- | Write an Int to a raw MutableByteArray at element index.
writeRawI :: MutableByteArray s -> Int -> Int -> Int -> ST s ()
writeRawI (MutableByteArray mba) (I# off) (I# i) (I# v) = ST $ \s ->
  case writeIntArray# mba (off +# i) v s of
    s' -> (# s', () #)
{-# INLINE writeRawI #-}

-- | Raw primop QR loop: all diagonal/subdiagonal access via raw ByteArray# primops.
rawTridiagQRLoop :: MutableByteArray s -> Int   -- ^ diagonal array + offset
                 -> MutableByteArray s -> Int   -- ^ subdiagonal array + offset
                 -> MutableByteArray s -> Int   -- ^ Q matrix + offset
                 -> Int -> Int -> Double -> ST s ()
rawTridiagQRLoop mbaD offD mbaSD offSD mbaQ offQ nn maxIter tol = go 0 0 (nn - 1)
  where
    go !iter !lo !hi
      | iter >= maxIter = pure ()
      | lo >= hi = pure ()
      | otherwise = do
          -- Bottom deflation
          sdhi <- readRawD mbaSD offSD (hi - 1)
          dhi1 <- readRawD mbaD offD (hi - 1)
          dhi  <- readRawD mbaD offD hi
          if abs sdhi <= tol * (abs dhi1 + abs dhi)
            then do writeRawD mbaSD offSD (hi - 1) 0; go iter lo (hi - 1)
            else do
              -- Top deflation
              sdlo <- readRawD mbaSD offSD lo
              dlo  <- readRawD mbaD offD lo
              dlo1 <- readRawD mbaD offD (lo + 1)
              if abs sdlo <= tol * (abs dlo + abs dlo1)
                then do writeRawD mbaSD offSD lo 0; go iter (lo + 1) hi
                else do
                  -- Interior deflation: find split point
                  split <- rawFindSplit mbaD offD mbaSD offSD lo hi tol
                  case split of
                    Just q -> do
                      writeRawD mbaSD offSD q 0
                      go iter lo q
                      go iter (q + 1) hi
                    Nothing -> do
                      let sp1 = sdhi
                          delta = (dhi1 - dhi) / 2
                          sgn = if delta >= 0 then 1 else -1
                          shift = dhi - sp1*sp1 / (delta + sgn * sqrt (delta*delta + sp1*sp1))
                      rawImplicitQRStep mbaD offD mbaSD offSD mbaQ offQ nn shift lo hi
                      go (iter + 1) lo hi

-- | Raw primop interior split search.
rawFindSplit :: MutableByteArray s -> Int -> MutableByteArray s -> Int
             -> Int -> Int -> Double -> ST s (Maybe Int)
rawFindSplit mbaD offD mbaSD offSD lo hi tol = scan (lo + 1)
  where
    scan q
      | q >= hi - 1 = pure Nothing
      | otherwise = do
          sdq <- readRawD mbaSD offSD q
          dq  <- readRawD mbaD offD q
          dq1 <- readRawD mbaD offD (q + 1)
          if abs sdq <= tol * (abs dq + abs dq1)
            then pure (Just q)
            else scan (q + 1)

-- | Raw primop implicit QR step via bulge-chasing Givens rotations.
rawImplicitQRStep :: MutableByteArray s -> Int -> MutableByteArray s -> Int
                  -> MutableByteArray s -> Int
                  -> Int -> Double -> Int -> Int -> ST s ()
rawImplicitQRStep mbaD offD mbaSD offSD mbaQ offQ nn shift lo hi = do
  dlo <- readRawD mbaD offD lo
  sdlo <- readRawD mbaSD offSD lo
  chase lo (dlo - shift) sdlo
  where
    chase k x z = do
      let (c, s) = givensRotation x z
      when (k > lo) $
        writeRawD mbaSD offSD (k-1) (c * x - s * z)
      dk  <- readRawD mbaD offD k
      ek  <- readRawD mbaSD offSD k
      dk1 <- readRawD mbaD offD (k+1)
      writeRawD mbaD offD k     (c*c*dk - 2*c*s*ek + s*s*dk1)
      writeRawD mbaD offD (k+1) (s*s*dk + 2*c*s*ek + c*c*dk1)
      writeRawD mbaSD offSD k   (c*s*(dk - dk1) + (c*c - s*s)*ek)
      -- Raw primop Givens rotation on Q
      rawMutApplyGivensColumns mbaQ offQ nn c s k (k+1) nn
      if k + 1 < hi
        then do
          ek1 <- readRawD mbaSD offSD (k+1)
          let z' = -s * ek1
          writeRawD mbaSD offSD (k+1) (c * ek1)
          ek_new <- readRawD mbaSD offSD k
          chase (k+1) ek_new z'
        else pure ()

-- --------------------------------------------------------------------------
-- Divide-and-conquer tridiagonal eigensolver (GVL4 Section 8.4)
-- Optimised: pre-allocated workspace, QR fallback for small subproblems,
-- unsafeFreezeByteArray to avoid O(n²) copies.
-- --------------------------------------------------------------------------

-- | Optimised divide-and-conquer eigensolver for a symmetric tridiagonal matrix.
-- Pre-allocates all temporary arrays once (eliminating per-level GC pressure),
-- falls back to QR for subproblems ≤ 25, and uses unsafeFreezeByteArray for
-- O(1) GEMM input preparation.
dcEigenTridiagOpt :: MutableByteArray s -> Int   -- ^ d + offset
                  -> MutableByteArray s -> Int   -- ^ e + offset
                  -> MutableByteArray s -> Int   -- ^ Q + offset (fullN × fullN row-major)
                  -> Int                         -- ^ fullN (Q dimension)
                  -> Int -> Int                  -- ^ lo, hi (active range, inclusive)
                  -> Double                      -- ^ tolerance
                  -> ST s ()
dcEigenTridiagOpt mbaD offD mbaE offE mbaQ offQ fullN lo0 hi0 tol
  | hi0 <= lo0 = pure ()
  | otherwise = do
      let !maxN = hi0 - lo0 + 1
      -- Pre-allocate all workspace at maximum needed size (once, not per-level)
      wsLam    <- newByteArray (maxN * 8)
      wsZ      <- newByteArray (maxN * 8)
      wsDSort  <- newByteArray (maxN * 8)
      wsZSort  <- newByteArray (maxN * 8)
      wsIdx    <- newByteArray (maxN * 8)
      wsPerm   <- newByteArray (maxN * 8)  -- deflation permutation (Int indices)
      wsW      <- newByteArray (maxN * maxN * 8)
      wsQsub   <- newByteArray (fullN * maxN * 8)
      wsResult <- newByteArray (fullN * maxN * 8)
      wsQtemp  <- newByteArray (maxN * maxN * 8)

      let !dcThreshold = 25

          -- Main recursive function (captures workspace via closure)
          dcGo !lo !hi
            | lo >= hi = pure ()
            | hi == lo + 1 = do
                -- 2×2 direct eigensolve
                d0 <- readRawD mbaD offD lo
                d1 <- readRawD mbaD offD hi
                e0 <- readRawD mbaE offE lo
                let !tr = d0 + d1
                    !det_ = d0 * d1 - e0 * e0
                    !disc = sqrt (max 0 (tr * tr - 4 * det_))
                    !lam1 = (tr - disc) / 2
                    !lam2 = (tr + disc) / 2
                    (!c, !s) = if abs e0 < tol * (abs d0 + abs d1)
                               then (1, 0)
                               else let !theta = (d1 - d0) / (2 * e0)
                                        !t_ = if theta >= 0
                                              then 1 / (theta + sqrt (1 + theta * theta))
                                              else 1 / (theta - sqrt (1 + theta * theta))
                                        !c_ = 1 / sqrt (1 + t_ * t_)
                                    in (c_, t_ * c_)
                writeRawD mbaD offD lo lam1
                writeRawD mbaD offD hi lam2
                writeRawD mbaE offE lo 0
                rawMutApplyGivensColumns mbaQ offQ fullN c s lo hi fullN
            | hi - lo + 1 <= dcThreshold = do
                -- QR fallback for small subproblems
                let !k = hi - lo + 1
                -- Initialise wsQtemp as k×k identity
                forM_ [0..k*k-1] $ \i -> writeRawD wsQtemp 0 i 0
                forM_ [0..k-1] $ \i -> writeRawD wsQtemp 0 (i * k + i) 1
                -- Run QR iteration on d[lo..hi], e[lo..hi-1]
                rawTridiagQRLoop mbaD (offD + lo) mbaE (offE + lo) wsQtemp 0 k (30 * k) tol
                -- Apply rotation: Q[:, lo..lo+k-1] = Q[:, lo..lo+k-1] * wsQtemp
                applyRotToQ lo k wsQtemp
            | otherwise = do
                -- D&C merge for larger subproblems
                let !k  = (lo + hi) `div` 2
                    !n1 = k - lo + 1
                    !n2 = hi - k
                    !nn = hi - lo + 1

                -- Read and modify the coupling element
                beta <- readRawD mbaE offE k
                dk   <- readRawD mbaD offD k
                dk1  <- readRawD mbaD offD (k + 1)
                let !absBeta = abs beta
                    !rho = beta
                writeRawD mbaD offD k     (dk - absBeta)
                writeRawD mbaD offD (k+1) (dk1 - absBeta)
                writeRawD mbaE offE k 0

                -- Recurse on T1 [lo..k] and T2 [k+1..hi]
                dcGo lo k
                dcGo (k + 1) hi

                -- === Merge phase ===
                -- Extract z vector from Q rows
                forM_ [0..n1-1] $ \i -> do
                  qv <- readRawD mbaQ 0 (offQ + k * fullN + (lo + i))
                  writeRawD wsZ 0 i qv
                forM_ [0..n2-1] $ \i -> do
                  qv <- readRawD mbaQ 0 (offQ + (k+1) * fullN + (lo + n1 + i))
                  writeRawD wsZ 0 (n1 + i) qv

                -- Copy d[lo..hi] and z into sortable arrays with indices
                forM_ [0..nn-1] $ \i -> do
                  di <- readRawD mbaD offD (lo + i)
                  writeRawD wsDSort 0 i di
                  writeRawD wsZSort 0 i =<< readRawD wsZ 0 i
                  writeRawD wsIdx 0 i (fromIntegral i)

                -- Sort by d values (insertion sort)
                forM_ [1..nn-1] $ \i -> do
                  di   <- readRawD wsDSort 0 i
                  zi   <- readRawD wsZSort 0 i
                  idxi <- readRawD wsIdx 0 i
                  let insertAt !j
                        | j < 0 = do
                            writeRawD wsDSort 0 0 di
                            writeRawD wsZSort 0 0 zi
                            writeRawD wsIdx 0 0 idxi
                        | otherwise = do
                            dj <- readRawD wsDSort 0 j
                            if dj > di
                              then do
                                writeRawD wsDSort 0 (j+1) dj
                                writeRawD wsZSort 0 (j+1) =<< readRawD wsZSort 0 j
                                writeRawD wsIdx 0 (j+1) =<< readRawD wsIdx 0 j
                                insertAt (j - 1)
                              else do
                                writeRawD wsDSort 0 (j+1) di
                                writeRawD wsZSort 0 (j+1) zi
                                writeRawD wsIdx 0 (j+1) idxi
                  insertAt (i - 1)

                -- Deflation: partition into non-deflated and deflated
                zn2 <- sumZSq wsZSort 0 nn
                let !deflTol = tol * sqrt zn2
                kND <- deflatePartition wsZSort 0 wsPerm 0 nn deflTol

                -- Extract Q_sub with permuted columns into wsQsub (fullN × nn)
                forM_ [0..nn-1] $ \sortedJ -> do
                  origIdx <- readRawD wsIdx 0 sortedJ
                  let !origJ = round origIdx :: Int
                      !srcCol = lo + origJ
                  forM_ [0..fullN-1] $ \row -> do
                    qv <- readRawD mbaQ 0 (offQ + row * fullN + srcCol)
                    writeRawD wsQsub 0 (row * nn + sortedJ) qv

                if kND == 0
                  then do
                    -- All deflated: eigenvalues = sorted d, eigenvectors = sorted Q cols
                    forM_ [0..nn-1] $ \i ->
                      forM_ [0..fullN-1] $ \row -> do
                        qv <- readRawD wsQsub 0 (row * nn + i)
                        writeRawD mbaQ 0 (offQ + row * fullN + (lo + i)) qv
                    forM_ [0..nn-1] $ \i -> do
                      di <- readRawD wsDSort 0 i
                      writeRawD mbaD offD (lo + i) di

                  else if kND == nn
                    then do
                      -- No deflation: full secular solve + full GEMM (original path)
                      secularSolve wsLam 0 wsDSort 0 wsZSort 0 rho nn deflTol
                      dcEigenvectors wsW 0 wsDSort 0 wsZSort 0 wsLam 0 rho nn
                      baQsub <- unsafeFreezeByteArray wsQsub
                      baW    <- unsafeFreezeByteArray wsW
                      forM_ [0..fullN * nn - 1] $ \i -> writeRawD wsResult 0 i 0
                      rawGemmKernel baQsub 0 baW 0 wsResult 0 fullN nn nn
                      forM_ [0..nn-1] $ \i ->
                        forM_ [0..fullN-1] $ \row -> do
                          rv <- readRawD wsResult 0 (row * nn + i)
                          writeRawD mbaQ 0 (offQ + row * fullN + (lo + i)) rv
                      forM_ [0..nn-1] $ \i -> do
                        lam <- readRawD wsLam 0 i
                        writeRawD mbaD offD (lo + i) lam

                    else do
                      -- Partial deflation: reduced secular solve + reduced GEMM
                      -- Build compressed d_nd[0..kND-1] and z_nd[0..kND-1]
                      -- Store in wsQtemp: d_nd at offset 0, z_nd at offset kND
                      forM_ [0..kND-1] $ \j -> do
                        pi_ <- readRawI wsPerm 0 j
                        dpi <- readRawD wsDSort 0 pi_
                        zpi <- readRawD wsZSort 0 pi_
                        writeRawD wsQtemp 0 j dpi
                        writeRawD wsQtemp 0 (kND + j) zpi

                      -- Solve kND secular equations on compressed system
                      secularSolve wsLam 0 wsQtemp 0 wsQtemp kND rho kND deflTol

                      -- Compute kND×kND eigenvector matrix W_nd
                      dcEigenvectors wsW 0 wsQtemp 0 wsQtemp kND wsLam 0 rho kND

                      -- Copy deflated columns from wsQsub to mbaQ
                      forM_ [kND..nn-1] $ \j -> do
                        pi_ <- readRawI wsPerm 0 j
                        forM_ [0..fullN-1] $ \row -> do
                          qv <- readRawD wsQsub 0 (row * nn + pi_)
                          writeRawD mbaQ 0 (offQ + row * fullN + (lo + j)) qv

                      -- Extract Q_nd (fullN×kND) from non-deflated columns
                      forM_ [0..kND-1] $ \j -> do
                        pi_ <- readRawI wsPerm 0 j
                        forM_ [0..fullN-1] $ \row -> do
                          qv <- readRawD wsQsub 0 (row * nn + pi_)
                          writeRawD wsResult 0 (row * kND + j) qv

                      -- GEMM: wsQsub(fullN×kND) = Q_nd(fullN×kND) × W_nd(kND×kND)
                      baQnd <- unsafeFreezeByteArray wsResult
                      baW   <- unsafeFreezeByteArray wsW
                      forM_ [0..fullN * kND - 1] $ \i -> writeRawD wsQsub 0 i 0
                      rawGemmKernel baQnd 0 baW 0 wsQsub 0 fullN kND kND

                      -- Copy GEMM result (non-deflated columns) to mbaQ
                      forM_ [0..kND-1] $ \j ->
                        forM_ [0..fullN-1] $ \row -> do
                          rv <- readRawD wsQsub 0 (row * kND + j)
                          writeRawD mbaQ 0 (offQ + row * fullN + (lo + j)) rv

                      -- Write eigenvalues: non-deflated from wsLam, deflated from wsDSort
                      forM_ [0..kND-1] $ \i -> do
                        lam <- readRawD wsLam 0 i
                        writeRawD mbaD offD (lo + i) lam
                      forM_ [kND..nn-1] $ \j -> do
                        pi_ <- readRawI wsPerm 0 j
                        di <- readRawD wsDSort 0 pi_
                        writeRawD mbaD offD (lo + j) di

          -- Apply a k×k rotation matrix to Q columns [lo..lo+k-1] via GEMM
          applyRotToQ !lo !k rotMat = do
            -- Extract Q[:, lo..lo+k-1] into wsQsub (fullN × k)
            forM_ [0..k-1] $ \j ->
              forM_ [0..fullN-1] $ \row -> do
                qv <- readRawD mbaQ 0 (offQ + row * fullN + (lo + j))
                writeRawD wsQsub 0 (row * k + j) qv
            -- O(1) freeze for GEMM inputs
            baQsub <- unsafeFreezeByteArray wsQsub
            baRot  <- unsafeFreezeByteArray rotMat
            -- Zero result
            forM_ [0..fullN * k - 1] $ \i -> writeRawD wsResult 0 i 0
            -- GEMM: result(fullN×k) = Qsub(fullN×k) * Rot(k×k)
            rawGemmKernel baQsub 0 baRot 0 wsResult 0 fullN k k
            -- Copy result back to Q
            forM_ [0..k-1] $ \j ->
              forM_ [0..fullN-1] $ \row -> do
                rv <- readRawD wsResult 0 (row * k + j)
                writeRawD mbaQ 0 (offQ + row * fullN + (lo + j)) rv

      dcGo lo0 hi0

-- | Solve the secular equation: f(λ) = 1 + ρ * Σ z[i]² / (d[i] - λ) = 0
-- for all nn roots. Roots are stored in mbaLam.
-- d must be sorted in ascending order.
secularSolve :: MutableByteArray s -> Int    -- ^ output eigenvalues
             -> MutableByteArray s -> Int    -- ^ sorted d (poles)
             -> MutableByteArray s -> Int    -- ^ sorted z
             -> Double -> Int -> Double      -- ^ rho, n, deflation tolerance
             -> ST s ()
secularSolve mbaLam offLam mbaD offD mbaZ offZ rho nn deflTol = do
  forM_ [0..nn-1] $ \i -> do
    zi <- readRawD mbaZ offZ i
    if abs zi <= deflTol
      then do
        -- Deflated: eigenvalue = d[i]
        di <- readRawD mbaD offD i
        writeRawD mbaLam offLam i di
      else do
        lam <- secularSolveOne mbaD offD mbaZ offZ rho i nn
        writeRawD mbaLam offLam i lam

-- | Solve one root of the secular equation between d[j] and d[j+1]
-- (or between d[n-1] and +infinity for the last root when rho > 0,
-- or between -infinity and d[0] for the first root when rho < 0).
-- Uses the Gragg/Borges fixed-weight quadratic method (cf. LAPACK dlasd4):
-- splits f(λ) at the two closest poles, approximates far terms as constant,
-- and solves the resulting quadratic for rapid convergence (2–4 iterations).
secularSolveOne :: MutableByteArray s -> Int -> MutableByteArray s -> Int
                -> Double -> Int -> Int -> ST s Double
secularSolveOne mbaD offD mbaZ offZ rho j nn = do
  dj <- readRawD mbaD offD j
  zj <- readRawD mbaZ offZ j
  let !zj2 = zj * zj
  -- Determine bracket and second pole
  if rho > 0
    then if j < nn - 1
      then do
        dj1 <- readRawD mbaD offD (j+1)
        zj1 <- readRawD mbaZ offZ (j+1)
        -- Interior root between d[j] and d[j+1]
        let !gap = dj1 - dj
            !mid = dj + gap * 0.5
        fixedWeightLoop 0 mid dj dj1 gap zj2 (zj1 * zj1)
      else do
        -- Last root when rho > 0: between d[n-1] and d[n-1] + rho*||z||²
        zn2 <- sumZSq mbaZ offZ nn
        let !hi_ = dj + abs rho * zn2
            !mid = (dj + hi_) * 0.5
        -- For edge root, use Newton+bisection (no second close pole)
        newtonLoop 0 mid dj hi_
    else if j > 0
      then do
        dj0 <- readRawD mbaD offD (j-1)
        zj0 <- readRawD mbaZ offZ (j-1)
        -- Interior root between d[j-1] and d[j] (rho < 0)
        let !gap = dj - dj0
            !mid = dj0 + gap * 0.5
        fixedWeightLoop 0 mid dj0 dj gap (zj0 * zj0) zj2
      else do
        -- First root when rho < 0
        zn2 <- sumZSq mbaZ offZ nn
        let !lo_ = dj - abs rho * zn2
            !mid = (lo_ + dj) * 0.5
        newtonLoop 0 mid lo_ dj
  where
    !maxIter_ = 30 :: Int

    -- Fixed-weight quadratic iteration for interior roots.
    -- dLo, dHi are the two closest poles; gap = dHi - dLo.
    -- z2Lo, z2Hi are z²[lo_pole] and z²[hi_pole].
    fixedWeightLoop !iter !lam !lb !ub !gap !z2Lo !z2Hi
      | iter >= maxIter_ = pure lam
      | otherwise = do
          -- Evaluate f(λ) with split at dLo and dHi
          (psiSum, phiSum) <- secularFuncSplit mbaD offD mbaZ offZ nn lam lb ub
          let !f = 1 + rho * (psiSum + phiSum)
          if abs f < 1e-15 * (1 + abs lam)
            then pure lam
            else do
              -- Extract close-pole contributions
              let !deltaLo = lb - lam   -- d[lo_pole] - λ (negative for interior root)
                  !deltaHi = ub - lam   -- d[hi_pole] - λ (positive for interior root)
                  -- Protect against division by zero near poles
                  !aClose = if abs deltaLo > 1e-300 then z2Lo / deltaLo else 0
                  !bClose = if abs deltaHi > 1e-300 then z2Hi / deltaHi else 0
                  -- "Far" residual: W = f - ρ*(aClose + bClose)
                  !w = f - rho * (aClose + bClose)
                  -- Quadratic in τ = dLo - λ (= deltaLo):
                  -- W*τ² - (W*gap + ρ*z2Lo + ρ*z2Hi)*τ + ρ*z2Lo*gap = 0
                  !qa = w
                  !qb = -(w * gap + rho * z2Lo + rho * z2Hi)
                  !qc = rho * z2Lo * gap
                  !disc = qb * qb - 4 * qa * qc
              if disc < 0 || abs qa < 1e-300
                then do
                  -- Degenerate: fall back to bisection
                  let !lamNew = (lb + ub) * 0.5
                      !(lb', ub') = if f * rho > 0 then (lb, lam) else (lam, ub)
                  fixedWeightLoop (iter + 1) lamNew lb' ub' gap z2Lo z2Hi
                else do
                  let !sqrtDisc = sqrt disc
                      -- Two roots for τ = dLo - λ, i.e. λ = dLo - τ
                      -- Use the numerically stable form
                      !tauA = if qb <= 0
                              then (-qb + sqrtDisc) / (2 * qa)
                              else 2 * qc / (-qb + sqrtDisc)
                      !tauB = if qb <= 0
                              then 2 * qc / (-qb + sqrtDisc)
                              else (-qb + sqrtDisc) / (2 * qa)
                      -- λ = dLo - τ; pick the root in (dLo, dHi) i.e. τ in (-gap, 0)
                      !lamA = lb - tauA
                      !lamB = lb - tauB
                      !lamNew0 = if lamA > lb && lamA < ub then lamA
                                 else if lamB > lb && lamB < ub then lamB
                                 else (lb + ub) * 0.5  -- bisection fallback
                      -- Update bracket
                      !(lb', ub') = if f * rho > 0 then (lb, lam) else (lam, ub)
                      -- Ensure lamNew is in updated bracket
                      !lamNew = if lamNew0 > lb' && lamNew0 < ub'
                                then lamNew0
                                else (lb' + ub') * 0.5
                  if abs (lamNew - lam) < 1e-15 * (1 + abs lam)
                    then pure lamNew
                    else fixedWeightLoop (iter + 1) lamNew lb' ub' gap z2Lo z2Hi

    -- Newton+bisection fallback for edge roots (first/last eigenvalue).
    newtonLoop !iter !lam !lb !ub
      | iter >= maxIter_ = pure lam
      | otherwise = do
          (f, fp) <- secularFuncAndDeriv mbaD offD mbaZ offZ rho nn lam
          if abs f < 1e-15 * (1 + abs lam)
            then pure lam
            else do
              let !step = f / fp
                  !lamNew0 = lam - step
                  !lamNew = if lamNew0 <= lb || lamNew0 >= ub
                            then (lb + ub) * 0.5
                            else lamNew0
                  !(lb', ub') = if f > 0 then (lb, lam) else (lam, ub)
              if abs (lamNew - lam) < 1e-15 * (1 + abs lam)
                then pure lamNew
                else newtonLoop (iter + 1) lamNew lb' ub'

-- | Evaluate the secular function split at the two bracket poles.
-- Returns (ψ, φ) where f(λ) = 1 + ρ*(ψ + φ).
-- ψ = Σ_{d[i] ≤ dLo} z[i]²/(d[i] - λ), φ = Σ_{d[i] ≥ dHi} z[i]²/(d[i] - λ)
secularFuncSplit :: MutableByteArray s -> Int -> MutableByteArray s -> Int
                 -> Int -> Double -> Double -> Double -> ST s (Double, Double)
secularFuncSplit mbaD offD mbaZ offZ nn lam dLo dHi = go 0 0 0
  where
    go i !psiAcc !phiAcc
      | i >= nn = pure (psiAcc, phiAcc)
      | otherwise = do
          di <- readRawD mbaD offD i
          zi <- readRawD mbaZ offZ i
          let !diff = di - lam
              !zi2 = zi * zi
          if abs diff < 1e-300
            then go (i+1) psiAcc phiAcc
            else let !term = zi2 / diff
                 in if di <= dLo
                    then go (i+1) (psiAcc + term) phiAcc
                    else go (i+1) psiAcc (phiAcc + term)

-- | Evaluate the secular function f(λ) = 1 + ρ * Σ z[i]² / (d[i] - λ)
-- and its derivative f'(λ) = ρ * Σ z[i]² / (d[i] - λ)².
secularFuncAndDeriv :: MutableByteArray s -> Int -> MutableByteArray s -> Int
                    -> Double -> Int -> Double -> ST s (Double, Double)
secularFuncAndDeriv mbaD offD mbaZ offZ rho nn lam = do
  (fSum, fpSum) <- go 0 0 0
  pure (1 + rho * fSum, rho * fpSum)
  where
    go i !fAcc !fpAcc
      | i >= nn = pure (fAcc, fpAcc)
      | otherwise = do
          di <- readRawD mbaD offD i
          zi <- readRawD mbaZ offZ i
          let diff = di - lam
              zi2 = zi * zi
          if abs diff < 1e-300
            then go (i+1) fAcc fpAcc  -- skip near-pole
            else go (i+1) (fAcc + zi2 / diff) (fpAcc + zi2 / (diff * diff))

-- | Sum of squares of z vector.
sumZSq :: MutableByteArray s -> Int -> Int -> ST s Double
sumZSq mbaZ offZ nn = go 0 0
  where
    go i !acc
      | i >= nn = pure acc
      | otherwise = do
          zi <- readRawD mbaZ offZ i
          go (i+1) (acc + zi * zi)

-- | Partition sorted indices into non-deflated (|z[i]| > deflTol) and deflated.
-- Returns k (non-deflated count).
-- perm[0..k-1] = sorted indices of non-deflated entries (in sorted order).
-- perm[k..nn-1] = sorted indices of deflated entries (in sorted order).
deflatePartition :: MutableByteArray s -> Int    -- ^ sorted z + offset
                 -> MutableByteArray s -> Int    -- ^ output perm (Int array) + offset
                 -> Int -> Double                -- ^ nn, deflTol
                 -> ST s Int
deflatePartition mbaZ offZ mbaPerm offPerm nn deflTol = do
    k <- goND 0 0
    goDF 0 k
    pure k
  where
    goND !i !kND
      | i >= nn = pure kND
      | otherwise = do
          zi <- readRawD mbaZ offZ i
          if abs zi > deflTol
            then do
              writeRawI mbaPerm offPerm kND i
              goND (i+1) (kND+1)
            else goND (i+1) kND
    goDF !i !pos
      | i >= nn = pure ()
      | otherwise = do
          zi <- readRawD mbaZ offZ i
          if abs zi <= deflTol
            then do
              writeRawI mbaPerm offPerm pos i
              goDF (i+1) (pos+1)
            else goDF (i+1) pos

-- | Compute eigenvector matrix W from secular equation solutions.
-- W[j,i] = z[j] / (d[j] - lambda[i]), each column normalised.
-- Single-pass: writes unnormalised entries and accumulates norm² simultaneously,
-- then normalises each column with SIMD.
dcEigenvectors :: MutableByteArray s -> Int     -- ^ W (nn × nn output)
               -> MutableByteArray s -> Int     -- ^ d (sorted poles)
               -> MutableByteArray s -> Int     -- ^ z
               -> MutableByteArray s -> Int     -- ^ lambda (eigenvalues)
               -> Double -> Int                 -- ^ rho, nn
               -> ST s ()
dcEigenvectors mbaW offW mbaD offD mbaZ offZ mbaLam offLam _rho nn = do
  forM_ [0..nn-1] $ \i -> do
    lami <- readRawD mbaLam offLam i
    -- Single pass: write unnormalised W[j,i] and accumulate norm²
    norm2 <- writeAndNorm lami i 0 0
    let !invNorm = if norm2 > 0 then 1 / sqrt norm2 else 1
    -- SIMD normalisation: W[j,i] *= invNorm for all j
    -- W is stored row-major W[j,i] at offset j*nn+i, so column i has stride nn.
    -- Since column stride is nn (not 1), we use scalar normalisation.
    forM_ [0..nn-1] $ \j -> do
      wji <- readRawD mbaW offW (j * nn + i)
      writeRawD mbaW offW (j * nn + i) (wji * invNorm)
  where
    writeAndNorm !lami !i !j !acc
      | j >= nn = pure acc
      | otherwise = do
          zj <- readRawD mbaZ offZ j
          dj <- readRawD mbaD offD j
          let !diff = dj - lami
              !w = if abs diff < 1e-300 then 0 else zj / diff
          writeRawD mbaW offW (j * nn + i) w
          writeAndNorm lami i (j + 1) (acc + w * w)

-- | Classical Jacobi eigenvalue method (GVL4 Section 8.5).
jacobiEigen :: forall n r e. (KnownNat n, M.Manifest r e, Floating e, Ord e)
            => Matrix n n r e -> Int -> e -> (Vector n r e, Matrix n n r e)
jacobiEigen a maxSweeps tol =
  let nn = dimVal @n
      (eigvals, qArr) = M.withMArrayST (unMatrix (identityMatrix @n @r)) $ \mq -> do
        ma <- M.thawS (unMatrix a)
        jacobiLoop ma mq nn maxSweeps tol
        evs <- mapM (\i -> M.readM ma (i :. i)) [0..nn-1]
        pure (makeVector @n @r $ \i -> evs !! i)
  in (eigvals, MkMatrix qArr)

jacobiLoop :: (M.Manifest r e, Floating e, Ord e)
           => M.MArray s r Ix2 e -> M.MArray s r Ix2 e -> Int -> Int -> e -> ST s ()
jacobiLoop ma mq nn maxSweeps tol = go 0
  where
    go !sweep
      | sweep >= maxSweeps = pure ()
      | otherwise = do
          offNorm <- offDiagNormST ma nn
          if offNorm < tol then pure ()
          else do
            forM_ [(p_, q_) | p_ <- [0..nn-2], q_ <- [p_+1..nn-1]] $ \(p_, q_) -> do
              apq <- M.readM ma (p_ :. q_)
              when (abs apq > tol * 1e-3) $ do
                app <- M.readM ma (p_ :. p_)
                aqq <- M.readM ma (q_ :. q_)
                let (c, s) = jacobiRotation app apq aqq
                applyJacobiInPlace ma c s p_ q_ nn
                applyGivensRightQ mq c s p_ q_ nn
            go (sweep + 1)

offDiagNormST :: (M.Manifest r e, Floating e) => M.MArray s r Ix2 e -> Int -> ST s e
offDiagNormST ma nn = do
  s <- go 0 0 0
  pure (sqrt s)
  where go !i !j !acc
          | i >= nn = pure acc
          | j >= nn = go (i+1) 0 acc
          | i == j = go i (j+1) acc
          | otherwise = do v <- M.readM ma (i :. j); go i (j+1) (acc + v*v)

jacobiRotation :: (Floating e, Ord e) => e -> e -> e -> (e, e)
jacobiRotation app apq aqq
  | apq == 0 = (1, 0)
  | otherwise =
    let tau = (aqq - app) / (2 * apq)
        t = if tau >= 0
            then 1 / (tau + sqrt (1 + tau * tau))
            else 1 / (tau - sqrt (1 + tau * tau))
        c = 1 / sqrt (1 + t * t)
        s = t * c
    in (c, s)

applyJacobiInPlace :: (M.Manifest r e, Num e)
                   => M.MArray s r Ix2 e -> e -> e -> Int -> Int -> Int -> ST s ()
applyJacobiInPlace ma c s p q nn = do
  app <- M.readM ma (p :. p)
  apq_ <- M.readM ma (p :. q)
  aqq <- M.readM ma (q :. q)
  M.write_ ma (p :. p) (c*c*app - 2*s*c*apq_ + s*s*aqq)
  M.write_ ma (q :. q) (s*s*app + 2*s*c*apq_ + c*c*aqq)
  M.write_ ma (p :. q) 0
  M.write_ ma (q :. p) 0
  forM_ [0..nn-1] $ \i -> when (i /= p && i /= q) $ do
    aip <- M.readM ma (i :. p)
    aiq <- M.readM ma (i :. q)
    let aip_new = c * aip - s * aiq
        aiq_new = s * aip + c * aiq
    M.write_ ma (i :. p) aip_new
    M.write_ ma (p :. i) aip_new
    M.write_ ma (i :. q) aiq_new
    M.write_ ma (q :. i) aiq_new
