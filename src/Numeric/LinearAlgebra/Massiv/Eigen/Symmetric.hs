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
import Data.Primitive.ByteArray (MutableByteArray(..), ByteArray(..), newByteArray, freezeByteArray)

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
-- --------------------------------------------------------------------------

-- | Divide-and-conquer eigensolver for a symmetric tridiagonal matrix.
-- Operates on mutable diagonal d[lo..hi], subdiagonal e[lo..hi-1],
-- and eigenvector matrix Q (fullN × fullN). On entry, Q[lo..hi, lo..hi]
-- should be identity. On exit, d[lo..hi] contains eigenvalues and
-- Q[:, lo..hi] contains corresponding eigenvectors.
dcEigenTridiag :: MutableByteArray s -> Int   -- ^ d + offset
               -> MutableByteArray s -> Int   -- ^ e + offset
               -> MutableByteArray s -> Int   -- ^ Q + offset (fullN × fullN row-major)
               -> Int                         -- ^ fullN (Q dimension)
               -> Int -> Int                  -- ^ lo, hi (active range, inclusive)
               -> Double                      -- ^ tolerance
               -> ST s ()
dcEigenTridiag mbaD offD mbaE offE mbaQ offQ fullN lo hi tol
  | lo >= hi = pure ()  -- 1×1: eigenvalue is d[lo], eigenvector is e_lo (already identity)
  | hi == lo + 1 = do
      -- 2×2 direct solve
      d0 <- readRawD mbaD offD lo
      d1 <- readRawD mbaD offD hi
      e0 <- readRawD mbaE offE lo
      let tr = d0 + d1
          det_ = d0 * d1 - e0 * e0
          disc = sqrt (max 0 (tr * tr - 4 * det_))
          lam1 = (tr - disc) / 2   -- smaller eigenvalue
          lam2 = (tr + disc) / 2   -- larger eigenvalue
      -- Eigenvector rotation
      let (c, s) = if abs e0 < tol * (abs d0 + abs d1)
                   then (1, 0)
                   else let theta = (d1 - d0) / (2 * e0)
                            t_ = if theta >= 0
                                 then 1 / (theta + sqrt (1 + theta * theta))
                                 else 1 / (theta - sqrt (1 + theta * theta))
                            c_ = 1 / sqrt (1 + t_ * t_)
                        in (c_, t_ * c_)
      writeRawD mbaD offD lo lam1
      writeRawD mbaD offD hi lam2
      writeRawD mbaE offE lo 0
      -- Apply Givens rotation to Q columns lo and hi
      rawMutApplyGivensColumns mbaQ offQ fullN c s lo hi fullN
  | otherwise = do
      let k = (lo + hi) `div` 2  -- split point
          n1 = k - lo + 1        -- size of T1
          n2 = hi - k            -- size of T2
          nn = hi - lo + 1       -- total active size

      -- Read the coupling element beta = e[k]
      beta <- readRawD mbaE offE k

      -- Modify diagonal: d[k] -= |beta|, d[k+1] -= |beta|
      dk  <- readRawD mbaD offD k
      dk1 <- readRawD mbaD offD (k + 1)
      let absBeta = abs beta
          rho = beta  -- sign matters for z vector
      writeRawD mbaD offD k     (dk - absBeta)
      writeRawD mbaD offD (k+1) (dk1 - absBeta)
      writeRawD mbaE offE k 0  -- zero out the coupling

      -- Recurse on T1 [lo..k] and T2 [k+1..hi]
      dcEigenTridiag mbaD offD mbaE offE mbaQ offQ fullN lo k tol
      dcEigenTridiag mbaD offD mbaE offE mbaQ offQ fullN (k+1) hi tol

      -- === Merge phase ===
      -- After recursion, d[lo..k] has eigenvalues of T1, d[k+1..hi] has eigenvalues of T2.
      -- Q columns lo..k and k+1..hi have been updated with eigenvectors.
      -- Form z = Q^T * v where v has 1 at position k (last row of Q1 block) and
      -- 1 at position k+1 (first row of Q2 block).
      -- z[i] = Q[k, lo+i] (for i in T1 range) or Q[k+1, lo+i] (for i in T2 range)

      -- Allocate temporary arrays
      mbaLam <- newByteArray (nn * 8)    -- new eigenvalues
      mbaZ   <- newByteArray (nn * 8)    -- z vector
      mbaDSort <- newByteArray (nn * 8)  -- sorted d values
      mbaZSort <- newByteArray (nn * 8)  -- sorted z values
      mbaIdx   <- newByteArray (nn * 8)  -- original index (stored as Double)

      -- Extract z vector: z[i] = Q[k, lo+i] for i=0..n1-1, Q[k+1, lo+i] for i=n1..nn-1
      forM_ [0..n1-1] $ \i -> do
        let qIdx = offQ + k * fullN + (lo + i)
        qv <- readRawD mbaQ 0 qIdx
        writeRawD mbaZ 0 i qv
      forM_ [0..n2-1] $ \i -> do
        let qIdx = offQ + (k+1) * fullN + (lo + n1 + i)
        qv <- readRawD mbaQ 0 qIdx
        writeRawD mbaZ 0 (n1 + i) qv

      -- Copy d[lo..hi] into a sortable array with indices
      forM_ [0..nn-1] $ \i -> do
        di <- readRawD mbaD offD (lo + i)
        writeRawD mbaDSort 0 i di
        writeRawD mbaZSort 0 i =<< readRawD mbaZ 0 i
        writeRawD mbaIdx 0 i (fromIntegral i)

      -- Sort by d values (insertion sort — fine for moderate nn)
      forM_ [1..nn-1] $ \i -> do
        di <- readRawD mbaDSort 0 i
        zi <- readRawD mbaZSort 0 i
        idxi <- readRawD mbaIdx 0 i
        let insertAt j
              | j < 0 = do
                  writeRawD mbaDSort 0 0 di
                  writeRawD mbaZSort 0 0 zi
                  writeRawD mbaIdx 0 0 idxi
              | otherwise = do
                  dj <- readRawD mbaDSort 0 j
                  if dj > di
                    then do
                      writeRawD mbaDSort 0 (j+1) dj
                      writeRawD mbaZSort 0 (j+1) =<< readRawD mbaZSort 0 j
                      writeRawD mbaIdx 0 (j+1) =<< readRawD mbaIdx 0 j
                      insertAt (j-1)
                    else do
                      writeRawD mbaDSort 0 (j+1) di
                      writeRawD mbaZSort 0 (j+1) zi
                      writeRawD mbaIdx 0 (j+1) idxi
        insertAt (i-1)

      -- Deflation: set near-zero z entries to zero
      let zNorm2 = sum <$> mapM (\i -> do zi <- readRawD mbaZSort 0 i; pure (zi*zi)) [0..nn-1]
      zn2 <- zNorm2
      let zNorm = sqrt zn2
          deflTol = tol * zNorm

      -- Count non-deflated entries and solve secular equation for them
      -- For simplicity, solve for all entries, treating deflated ones specially
      secularSolve mbaLam 0 mbaDSort 0 mbaZSort 0 rho nn deflTol

      -- Compute eigenvectors: W[j,i] = z[j] / (d[j] - lambda[i]), normalised
      mbaW <- newByteArray (nn * nn * 8)  -- nn×nn eigenvector matrix
      dcEigenvectors mbaW 0 mbaDSort 0 mbaZSort 0 mbaLam 0 rho nn

      -- Compose Q columns: Q_new[:, lo+i] = sum_j Q_old[:, lo+perm[j]] * W[j, i]
      -- This is Q_new = Q_old[:, lo..hi] * P * W where P is the sort permutation
      -- First, build Q_sub = Q_old[:, lo..hi] with columns permuted by sort order
      mbaQsub <- newByteArray (fullN * nn * 8)  -- fullN × nn
      forM_ [0..nn-1] $ \sortedJ -> do
        origIdx <- readRawD mbaIdx 0 sortedJ
        let origJ = round origIdx :: Int
            srcCol = lo + origJ
        forM_ [0..fullN-1] $ \row -> do
          qv <- readRawD mbaQ 0 (offQ + row * fullN + srcCol)
          writeRawD mbaQsub 0 (row * nn + sortedJ) qv

      -- Freeze Q_sub and W for GEMM
      baQsub <- freezeByteArray mbaQsub 0 (fullN * nn * 8)
      baW <- freezeByteArray mbaW 0 (nn * nn * 8)

      -- Allocate result: fullN × nn
      mbaResult <- newByteArray (fullN * nn * 8)
      -- Zero it
      forM_ [0..fullN * nn - 1] $ \i -> writeRawD mbaResult 0 i 0

      -- GEMM: result = Q_sub * W  (fullN × nn = fullN × nn * nn × nn)
      rawGemmKernel (ByteArray (unBA baQsub)) 0 (ByteArray (unBA baW)) 0 mbaResult 0 fullN nn nn

      -- Write result back to Q columns lo..hi
      forM_ [0..nn-1] $ \i ->
        forM_ [0..fullN-1] $ \row -> do
          rv <- readRawD mbaResult 0 (row * nn + i)
          writeRawD mbaQ 0 (offQ + row * fullN + (lo + i)) rv

      -- Write eigenvalues back to d[lo..hi]
      forM_ [0..nn-1] $ \i -> do
        lam <- readRawD mbaLam 0 i
        writeRawD mbaD offD (lo + i) lam

-- Helper to unwrap ByteArray for rawGemmKernel
unBA :: ByteArray -> ByteArray#
unBA (ByteArray ba) = ba
{-# INLINE unBA #-}

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
-- Uses the middle-way method with Newton-like iteration.
secularSolveOne :: MutableByteArray s -> Int -> MutableByteArray s -> Int
                -> Double -> Int -> Int -> ST s Double
secularSolveOne mbaD offD mbaZ offZ rho j nn = do
  dj <- readRawD mbaD offD j
  -- Determine bracket
  (lo_, hi_) <- if rho > 0
    then if j < nn - 1
      then do dj1 <- readRawD mbaD offD (j+1); pure (dj, dj1)
      else do
        -- Last root: bracket (d[n-1], d[n-1] + ||z||²)
        zn2 <- sumZSq mbaZ offZ nn
        pure (dj, dj + abs rho * zn2)
    else if j > 0
      then do dj0 <- readRawD mbaD offD (j-1); pure (dj0, dj)
      else do
        zn2 <- sumZSq mbaZ offZ nn
        pure (dj - abs rho * zn2, dj)

  -- Initial guess: midpoint
  let mid = (lo_ + hi_) / 2
  -- Newton-like iteration on the secular function
  iterate_ 0 mid lo_ hi_
  where
    maxIter_ = 80 :: Int

    iterate_ !iter lam lb ub
      | iter >= maxIter_ = pure lam
      | otherwise = do
          (f, fp) <- secularFuncAndDeriv mbaD offD mbaZ offZ rho nn lam
          if abs f < 1e-15 * (1 + abs lam)
            then pure lam
            else do
              -- Newton step with bisection fallback
              let step = f / fp
                  lamNew0 = lam - step
                  -- Ensure we stay in bracket
                  lamNew = if lamNew0 <= lb || lamNew0 >= ub
                           then (lb + ub) / 2
                           else lamNew0
                  -- Update bracket
                  (lb', ub') = if f > 0 then (lb, lam) else (lam, ub)
              if abs (lamNew - lam) < 1e-15 * (1 + abs lam)
                then pure lamNew
                else iterate_ (iter + 1) lamNew lb' ub'

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

-- | Compute eigenvector matrix W from secular equation solutions.
-- W[j,i] = z[j] / (d[j] - lambda[i]), each column normalised.
dcEigenvectors :: MutableByteArray s -> Int     -- ^ W (nn × nn output)
               -> MutableByteArray s -> Int     -- ^ d (sorted poles)
               -> MutableByteArray s -> Int     -- ^ z
               -> MutableByteArray s -> Int     -- ^ lambda (eigenvalues)
               -> Double -> Int                 -- ^ rho, nn
               -> ST s ()
dcEigenvectors mbaW offW mbaD offD mbaZ offZ mbaLam offLam _rho nn = do
  forM_ [0..nn-1] $ \i -> do
    lami <- readRawD mbaLam offLam i
    -- Compute column i: W[j,i] = z[j] / (d[j] - lami)
    -- First pass: compute unnormalised entries and norm
    norm2 <- go_norm i lami 0 0
    let invNorm = if norm2 > 0 then 1 / sqrt norm2 else 1
    -- Second pass: write normalised entries
    forM_ [0..nn-1] $ \j -> do
      zj <- readRawD mbaZ offZ j
      dj <- readRawD mbaD offD j
      let diff = dj - lami
          wji = if abs diff < 1e-300 then 0 else (zj / diff) * invNorm
      writeRawD mbaW offW (j * nn + i) wji
  where
    go_norm i lami j !acc
      | j >= nn = pure acc
      | otherwise = do
          zj <- readRawD mbaZ offZ j
          dj <- readRawD mbaD offD j
          let diff = dj - lami
              w = if abs diff < 1e-300 then 0 else zj / diff
          go_norm i lami (j+1) (acc + w * w)

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
