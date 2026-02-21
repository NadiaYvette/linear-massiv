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
  , symmetricEigenPDC
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
import Data.Primitive.ByteArray (ByteArray(..), MutableByteArray(..), newByteArray, unsafeFreezeByteArray)

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
      -- For n < panelCrossover: per-column Level-2 (rank-2 update per step)
      -- For n >= panelCrossover: DLATRD-style panel factorisation (Level-3 SYR2K)
      panelCrossover = 256
      (betaList, tArr) = M.withMArrayST (unMatrix a) $ \mt -> do
        let !mbaT = unwrapMutableByteArray mt
            !offT = unwrapMutableByteArrayOffset mt
        mbaV <- newByteArray (nn * 8)
        mbaP <- newByteArray (nn * 8)
        mbaW <- newByteArray (nn * 8)
        if nn < panelCrossover
          then do
            betas <- mapM (\k -> tridiagStepP mbaT offT nn mbaV mbaP mbaW k) [0..nn-3]
            pure betas
          else do
            let !nb = 32
                !numRef = nn - 2  -- number of Householder reflectors
            -- V_panel (nn × nb) and W_panel (nn × nb) for deferred rank-2 updates
            mbaVp <- newByteArray (nn * nb * 8)
            mbaWp <- newByteArray (nn * nb * 8)
            -- Temporary for GEMM-based trailing update
            mbaTemp <- newByteArray (nn * nb * 8)
            let go !k0 !accBetas
                  | k0 > numRef - 1 = pure (reverse accBetas)
                  | otherwise = do
                      let !bs = min nb (numRef - k0)
                      panelBetas <- panelTridiagP mbaT offT nn mbaV mbaP mbaW
                                                  mbaVp mbaWp mbaTemp k0 bs
                      go (k0 + bs) (reverse panelBetas ++ accBetas)
            go 0 []

      -- Get underlying ByteArray from frozen T for Q accumulation
      !tBA  = unwrapByteArray tArr
      !tOff = unwrapByteArrayOffset tArr

      -- Phase 2: Q accumulation.
      -- For n < 200: per-row Householder updates (minimal work, avoids GEMM overhead).
      -- For n >= 200: blocked WY with Level-3 GEMM (better cache/SIMD utilization).
      qMat = createMatrix @n @n @M.P $ \mq -> do
        let !mbaQ = unwrapMutableByteArray mq
            !offQ = unwrapMutableByteArrayOffset mq
        -- Set Q = I
        forM_ [0..nn-1] $ \i -> forM_ [0..nn-1] $ \j ->
          writeRawD mbaQ offQ (i*nn+j) (if i == j then 1 else 0)
        if nn < 256
          then
            -- Per-row approach: Q <- Q · H_k for k = 0..n-3
            forM_ (zip [0..] betaList) $ \(k, beta_k) ->
              when (beta_k /= 0) $
                forM_ [0..nn-1] $ \row ->
                  rawMutTridiagQAccum mbaQ offQ nn tBA tOff nn beta_k (k+1) k nn row
          else do
            -- Blocked WY approach: Q <- Q * (I - Y * T * Y^T) per block
            let !numRef = nn - 2
                !nb = min 32 numRef
            mbaBetas <- newByteArray (numRef * 8)
            forM_ (zip [0..] betaList) $ \(i, b) -> writeRawD mbaBetas 0 i b
            mbaY  <- newByteArray (nn * nb * 8)
            mbaTf <- newByteArray (nb * nb * 8)
            mbaW1 <- newByteArray (nn * nb * 8)
            mbaW2 <- newByteArray (nn * nb * 8)
            mbaYT <- newByteArray (nb * nn * 8)

            forM_ [0, nb .. numRef - 1] $ \k0 -> do
              let !bs = min nb (numRef - k0)

              -- Pack Y (n × bs) from stored Householder vectors
              forM_ [0..nn * bs - 1] $ \idx -> writeRawD mbaY 0 idx 0
              forM_ [0..bs-1] $ \j -> do
                let !k = k0 + j
                writeRawD mbaY 0 ((k+1) * bs + j) 1.0
                forM_ [k+2..nn-1] $ \l ->
                  writeRawD mbaY 0 (l * bs + j) (indexRawD tBA tOff (l * nn + k))

              -- Build T factor (bs × bs upper-triangular)
              forM_ [0..bs*bs-1] $ \idx -> writeRawD mbaTf 0 idx 0
              forM_ [0..bs-1] $ \j -> do
                betaj <- readRawD mbaBetas 0 (k0 + j)
                writeRawD mbaTf 0 (j * bs + j) betaj
                when (j > 0 && betaj /= 0) $ do
                  forM_ [0..j-1] $ \i -> do
                    let dotLoop !l !acc
                          | l >= nn = pure acc
                          | otherwise = do
                              yi <- readRawD mbaY 0 (l * bs + i)
                              yj <- readRawD mbaY 0 (l * bs + j)
                              dotLoop (l+1) (acc + yi * yj)
                    dot <- dotLoop (k0 + j + 1) 0
                    writeRawD mbaW1 0 i dot
                  forM_ [0..j-1] $ \i -> do
                    let triLoop !l !acc
                          | l >= j = pure acc
                          | otherwise = do
                              til <- readRawD mbaTf 0 (i * bs + l)
                              dl  <- readRawD mbaW1 0 l
                              triLoop (l+1) (acc + til * dl)
                    z <- triLoop i 0
                    writeRawD mbaTf 0 (i * bs + j) (negate betaj * z)

              -- W1 = Q · Y (GEMM n×n * n×bs → n×bs)
              baQ <- unsafeFreezeByteArray mbaQ
              baY <- unsafeFreezeByteArray mbaY
              forM_ [0..nn*bs-1] $ \idx -> writeRawD mbaW1 0 idx 0
              rawGemmKernel baQ offQ baY 0 mbaW1 0 nn nn bs

              -- W2 = W1 · T (GEMM n×bs * bs×bs → n×bs)
              baW1 <- unsafeFreezeByteArray mbaW1
              baTf <- unsafeFreezeByteArray mbaTf
              forM_ [0..nn*bs-1] $ \idx -> writeRawD mbaW2 0 idx 0
              rawGemmKernel baW1 0 baTf 0 mbaW2 0 nn bs bs

              -- Negate W2 in-place
              forM_ [0..nn*bs-1] $ \idx -> do
                v <- readRawD mbaW2 0 idx
                writeRawD mbaW2 0 idx (negate v)

              -- Transpose Y → Y^T (bs × n)
              forM_ [0..nn-1] $ \row ->
                forM_ [0..bs-1] $ \col ->
                  writeRawD mbaYT 0 (col * nn + row) (indexRawD baY 0 (row * bs + col))

              -- Q += (-W2) · Y^T (GEMM n×bs * bs×n → n×n)
              baNW2 <- unsafeFreezeByteArray mbaW2
              baYT <- unsafeFreezeByteArray mbaYT
              rawGemmKernel baNW2 0 baYT 0 mbaQ offQ nn bs nn

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

-- | DLATRD-style panel tridiagonalisation.
-- Processes columns k0..k0+bs-1, building V_panel and W_panel matrices
-- that represent the deferred rank-2 updates. After processing all columns
-- in the panel, applies a single Level-3 SYR2K trailing update.
--
-- Within the panel, column k of T is corrected for deferred updates:
--   T[:,k] -= V_panel * W_panel[k,:] + W_panel * V_panel[k,:]
-- before computing the Householder reflector.
--
-- Returns the list of beta values for the panel columns.
panelTridiagP :: MutableByteArray s -> Int -> Int  -- T matrix, offset, n
              -> MutableByteArray s -> MutableByteArray s -> MutableByteArray s  -- v, p, w temps
              -> MutableByteArray s -> MutableByteArray s -> MutableByteArray s  -- Vp, Wp, temp
              -> Int -> Int  -- k0, bs (panel start, panel size)
              -> ST s [Double]
panelTridiagP mbaT offT nn mbaV mbaP mbaW mbaVp mbaWp _mbaTemp k0 bs = do
  -- DLATRD-style: NO rank-2 updates to T within the panel.
  -- All corrections computed from V_panel, W_panel.
  -- After the panel, apply SYR2K to the full remaining submatrix.
  betas <- go 0 []

  -- After the panel: apply accumulated rank-2 update to the full remaining
  -- submatrix T[k0+1:nn, k0+1:nn]. This includes both within-panel diagonal
  -- entries and the trailing submatrix.
  --
  -- We must save/restore Householder vectors in columns k0..k0+bs-1
  -- because the SYR2K will overwrite them.
  let !remStart = k0 + 1
      !remSize = nn - remStart
  when (remSize > 0 && bs > 0) $ do
    -- Save Householder vectors from T columns k0..k0+bs-1
    -- These are T[i, k] for i > k+1, k in [k0..k0+bs-1]
    -- Also save subdiagonal entries T[k+1, k] = mu
    let !hvSize = bs * nn  -- upper bound for storage
    mbaHvSave <- newByteArray (hvSize * 8)
    forM_ [0..bs-1] $ \l -> do
      let !k = k0 + l
          !startRow = k + 1
      forM_ [startRow..nn-1] $ \i -> do
        val <- readRawD mbaT offT (i * nn + k)
        writeRawD mbaHvSave 0 (l * nn + i) val
      -- Also save T[k, k+1] (the upper subdiagonal)
      when (k + 1 < nn) $ do
        val <- readRawD mbaT offT (k * nn + (k + 1))
        writeRawD mbaHvSave 0 (l * nn + k) val  -- reuse slot k < startRow

    -- Build contiguous V_rem (remSize × bs) and W_rem (remSize × bs)
    mbaVr <- newByteArray (remSize * bs * 8)
    mbaWr <- newByteArray (remSize * bs * 8)
    forM_ [0..remSize-1] $ \i ->
      forM_ [0..bs-1] $ \j -> do
        readRawD mbaVp 0 ((remStart + i) * bs + j) >>= writeRawD mbaVr 0 (i * bs + j)
        readRawD mbaWp 0 ((remStart + i) * bs + j) >>= writeRawD mbaWr 0 (i * bs + j)

    -- Build -W_rem^T and -V_rem^T (bs × remSize)
    mbaNWrT <- newByteArray (bs * remSize * 8)
    mbaNVrT <- newByteArray (bs * remSize * 8)
    forM_ [0..remSize-1] $ \i ->
      forM_ [0..bs-1] $ \j -> do
        readRawD mbaWr 0 (i * bs + j) >>= \w -> writeRawD mbaNWrT 0 (j * remSize + i) (negate w)
        readRawD mbaVr 0 (i * bs + j) >>= \v -> writeRawD mbaNVrT 0 (j * remSize + i) (negate v)

    -- Copy T_rem to contiguous temp
    mbaRem <- newByteArray (remSize * remSize * 8)
    forM_ [0..remSize-1] $ \i ->
      forM_ [0..remSize-1] $ \j ->
        readRawD mbaT offT ((remStart + i) * nn + (remStart + j))
          >>= writeRawD mbaRem 0 (i * remSize + j)

    -- GEMM: rem += V_rem * (-W_rem^T) + W_rem * (-V_rem^T)
    baVr <- unsafeFreezeByteArray mbaVr
    baNWrT <- unsafeFreezeByteArray mbaNWrT
    rawGemmKernel baVr 0 baNWrT 0 mbaRem 0 remSize bs remSize
    baWr <- unsafeFreezeByteArray mbaWr
    baNVrT <- unsafeFreezeByteArray mbaNVrT
    rawGemmKernel baWr 0 baNVrT 0 mbaRem 0 remSize bs remSize

    -- Copy back to T
    forM_ [0..remSize-1] $ \i ->
      forM_ [0..remSize-1] $ \j ->
        readRawD mbaRem 0 (i * remSize + j)
          >>= writeRawD mbaT offT ((remStart + i) * nn + (remStart + j))

    -- Restore saved Householder vectors and subdiagonal entries
    forM_ [0..bs-1] $ \l -> do
      let !k = k0 + l
          !startRow = k + 1
      forM_ [startRow..nn-1] $ \i -> do
        val <- readRawD mbaHvSave 0 (l * nn + i)
        writeRawD mbaT offT (i * nn + k) val
      when (k + 1 < nn) $ do
        val <- readRawD mbaHvSave 0 (l * nn + k)
        writeRawD mbaT offT (k * nn + (k + 1)) val

  pure betas
  where
    go !j !acc
      | j >= bs = pure (reverse acc)
      | otherwise = do
          let !k = k0 + j
              !subSize = nn - k - 1

          -- Step 1: Read corrected column. T is ORIGINAL (no rank-2 updates applied).
          -- corrected_col[i] = T[i+k+1, k] - Σ_l (V[i+k+1,l]*W[k,l] + W[i+k+1,l]*V[k,l])
          forM_ [0..subSize-1] $ \i -> do
            tik <- readRawD mbaT offT ((k+1+i)*nn + k)
            if j == 0
              then writeRawD mbaP 0 i tik
              else do
                let corrLoop !l !accC
                      | l >= j = pure accC
                      | otherwise = do
                          vp_il <- readRawD mbaVp 0 ((k+1+i) * bs + l)
                          wp_kl <- readRawD mbaWp 0 (k * bs + l)
                          wp_il <- readRawD mbaWp 0 ((k+1+i) * bs + l)
                          vp_kl <- readRawD mbaVp 0 (k * bs + l)
                          corrLoop (l+1) (accC + vp_il * wp_kl + wp_il * vp_kl)
                corr <- corrLoop 0 0
                writeRawD mbaP 0 i (tik - corr)

          -- Step 2: Householder from corrected column
          x0 <- readRawD mbaP 0 0
          sigma <- do
            let sigLoop !i !acc_
                  | i >= subSize = pure acc_
                  | otherwise = do
                      ci <- readRawD mbaP 0 i
                      sigLoop (i+1) (acc_ + ci * ci)
            sigLoop 1 0
          if sigma == 0 && x0 >= 0
            then do
              forM_ [0..nn-1] $ \i -> do
                writeRawD mbaVp 0 (i * bs + j) 0
                writeRawD mbaWp 0 (i * bs + j) 0
              go (j+1) (0 : acc)
            else do
              let mu = sqrt (x0 * x0 + sigma)
                  v0 = if x0 <= 0 then x0 - mu else -sigma / (x0 + mu)
                  beta = 2 * v0 * v0 / (sigma + v0 * v0)

              writeRawD mbaV 0 0 1.0
              forM_ [1..subSize-1] $ \i -> do
                ci <- readRawD mbaP 0 i
                writeRawD mbaV 0 i (ci / v0)

              -- Step 3: p = beta * T * v (using ORIGINAL T, then correct via V,W)
              rawMutSymMatvecSub mbaT offT nn mbaV 0 mbaP 0 (k+1) nn
              forM_ [0..subSize-1] $ \i -> do
                pi_ <- readRawD mbaP 0 i
                writeRawD mbaP 0 i (beta * pi_)

              -- Step 4: Full V,W correction.
              -- p -= beta * (V_sub*(W_sub^T*v) + W_sub*(V_sub^T*v))
              -- where V_sub = V_panel[k+1:nn-1, 0:j-1], W_sub = W_panel[k+1:nn-1, 0:j-1]
              when (j > 0) $ do
                forM_ [0..j-1] $ \l -> do
                  let dotW !idx !accW
                        | idx >= subSize = pure accW
                        | otherwise = do
                            wp <- readRawD mbaWp 0 ((k+1+idx) * bs + l)
                            vi <- readRawD mbaV 0 idx
                            dotW (idx+1) (accW + wp * vi)
                  z1 <- dotW 0 0

                  let dotV !idx !accV
                        | idx >= subSize = pure accV
                        | otherwise = do
                            vp <- readRawD mbaVp 0 ((k+1+idx) * bs + l)
                            vi <- readRawD mbaV 0 idx
                            dotV (idx+1) (accV + vp * vi)
                  z2 <- dotV 0 0

                  forM_ [0..subSize-1] $ \i -> do
                    vp_il <- readRawD mbaVp 0 ((k+1+i) * bs + l)
                    wp_il <- readRawD mbaWp 0 ((k+1+i) * bs + l)
                    pi_ <- readRawD mbaP 0 i
                    writeRawD mbaP 0 i (pi_ - beta * (vp_il * z1 + wp_il * z2))

              -- Step 5: w = p - alpha*v
              ptv <- mutDotVec mbaP 0 mbaV 0 subSize
              let alpha_ = beta * ptv / 2
              forM_ [0..subSize-1] $ \i -> do
                pi_ <- readRawD mbaP 0 i
                vi  <- readRawD mbaV 0 i
                writeRawD mbaW 0 i (pi_ - alpha_ * vi)

              -- Step 6: Store v,w in panels
              forM_ [0..k] $ \i -> do
                writeRawD mbaVp 0 (i * bs + j) 0
                writeRawD mbaWp 0 (i * bs + j) 0
              forM_ [0..subSize-1] $ \i -> do
                readRawD mbaV 0 i >>= writeRawD mbaVp 0 ((k+1+i) * bs + j)
                readRawD mbaW 0 i >>= writeRawD mbaWp 0 ((k+1+i) * bs + j)

              -- Step 7: Store Householder vector in T (for Q accumulation)
              forM_ [1..subSize-1] $ \i ->
                readRawD mbaV 0 i >>= writeRawD mbaT offT ((k+1+i)*nn + k)

              -- Step 8: Set subdiagonal
              writeRawD mbaT offT ((k+1)*nn + k) mu
              writeRawD mbaT offT (k*nn + (k+1)) mu

              go (j+1) (beta : acc)

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
symmetricEigenP a maxIter tol
  | dimVal @n >= 100000 = symmetricEigenPDC a tol  -- disabled: D&C overhead too high
  | otherwise =
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

-- | Divide-and-conquer specialised symmetric eigenvalue decomposition for @P Double@.
-- Uses raw-primop tridiagonalisation then D&C eigensolver (GEMM-based merge)
-- instead of QR iteration. Faster than 'symmetricEigenP' at larger sizes (n ≥ 50).
symmetricEigenPDC :: forall n. KnownNat n
                  => Matrix n n M.P Double -> Double
                  -> (Vector n M.P Double, Matrix n n M.P Double)
symmetricEigenPDC a tol =
  let nn = dimVal @n
      (q0, diag_, subdiag) = tridiagonalizeP a
      (dArr, qArr) = M.withMArrayST (unMatrix q0) $ \mq -> do
        md  <- M.thawS (unVector diag_)
        msd <- M.thawS (unVector subdiag)
        let !mbaD  = unwrapMutableByteArray md
            !offD  = unwrapMutableByteArrayOffset md
            !mbaSD = unwrapMutableByteArray msd
            !offSD = unwrapMutableByteArrayOffset msd
            !mbaQ  = unwrapMutableByteArray mq
            !offQ  = unwrapMutableByteArrayOffset mq
        dcEigenTridiagOpt mbaD offD mbaSD offSD mbaQ offQ nn 0 (nn - 1) tol
        dFrozen <- M.freezeS md
        pure (MkVector dFrozen)
  in (dArr, MkMatrix qArr)
{-# NOINLINE symmetricEigenPDC #-}

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

-- | Read a Double from an immutable ByteArray at element index.
indexRawD :: ByteArray -> Int -> Int -> Double
indexRawD (ByteArray ba) (I# off) (I# i) =
  case indexDoubleArray# ba (off +# i) of
    v -> D# v
{-# INLINE indexRawD #-}

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
      rawMutApplyGivensColumns mbaQ offQ nn c (negate s) k (k+1) nn
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
--
-- The algorithm maintains a LOCAL eigenvector matrix (maxN × maxN, starting as
-- identity) throughout the D&C recursion. The z-vector for each merge step is
-- extracted from this local matrix (not the global Q), ensuring correctness.
-- At the end, the global Q is updated via a single GEMM: Q_out = Q_in * Q_local.
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
      wsQsub   <- newByteArray (maxN * maxN * 8)
      wsResult <- newByteArray (maxN * maxN * 8)
      wsQtemp  <- newByteArray (maxN * maxN * 8)

      -- Local eigenvector accumulator (maxN × maxN), initialised to identity.
      -- Indexed with LOCAL coordinates: row/col in [0..maxN-1].
      -- Local index i corresponds to global index (lo0 + i).
      wsQlocal <- newByteArray (maxN * maxN * 8)
      forM_ [0..maxN*maxN-1] $ \i -> writeRawD wsQlocal 0 i 0
      forM_ [0..maxN-1] $ \i -> writeRawD wsQlocal 0 (i * maxN + i) 1

      let !dcThreshold = 25

          -- Convenience: convert global index to local index
          toLocal !g = g - lo0

          -- Main recursive function (captures workspace via closure)
          -- All operations affect wsQlocal (local eigenvector accumulator),
          -- NOT the global Q matrix.
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
                -- Apply Givens to LOCAL eigenvector matrix columns
                rawMutApplyGivensColumns wsQlocal 0 maxN c (negate s) (toLocal lo) (toLocal hi) maxN
            | hi - lo + 1 <= dcThreshold = do
                -- QR fallback for small subproblems
                let !k = hi - lo + 1
                -- Initialise wsQtemp as k×k identity
                forM_ [0..k*k-1] $ \i -> writeRawD wsQtemp 0 i 0
                forM_ [0..k-1] $ \i -> writeRawD wsQtemp 0 (i * k + i) 1
                -- Run QR iteration on d[lo..hi], e[lo..hi-1]
                rawTridiagQRLoop mbaD (offD + lo) mbaE (offE + lo) wsQtemp 0 k (30 * k) tol
                -- Apply rotation to LOCAL eigenvector matrix:
                -- wsQlocal[:, toLocal(lo)..toLocal(lo)+k-1] *= wsQtemp
                applyRotToQlocal (toLocal lo) k wsQtemp
            | otherwise = do
                -- D&C merge for larger subproblems
                let !k  = (lo + hi) `div` 2
                    !n1 = k - lo + 1
                    !n2 = hi - k
                    !nn = hi - lo + 1
                    -- Local coordinates for the split
                    !kL  = toLocal k
                    !loL = toLocal lo

                -- Read and modify the coupling element
                beta <- readRawD mbaE offE k
                dk   <- readRawD mbaD offD k
                dk1  <- readRawD mbaD offD (k + 1)
                let !absBeta = abs beta
                    !rho = absBeta
                writeRawD mbaD offD k     (dk - absBeta)
                writeRawD mbaD offD (k+1) (dk1 - absBeta)
                writeRawD mbaE offE k 0

                -- Recurse on T1 [lo..k] and T2 [k+1..hi]
                dcGo lo k
                dcGo (k + 1) hi

                -- === Merge phase ===
                -- Extract z vector from LOCAL eigenvector matrix rows.
                -- z[0..n1-1] = last row of Q1 = row kL of wsQlocal, columns loL..loL+n1-1
                -- z[n1..nn-1] = first row of Q2 = row (kL+1) of wsQlocal, columns loL+n1..loL+nn-1
                forM_ [0..n1-1] $ \i -> do
                  qv <- readRawD wsQlocal 0 (kL * maxN + (loL + i))
                  writeRawD wsZ 0 i qv
                forM_ [0..n2-1] $ \i -> do
                  qv <- readRawD wsQlocal 0 ((kL + 1) * maxN + (loL + n1 + i))
                  let !zv = if beta < 0 then negate qv else qv
                  writeRawD wsZ 0 (n1 + i) zv

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
                -- Use a practical deflation tolerance: max of the user tolerance
                -- and a machine-epsilon-based criterion.  Entries with tiny z[i]
                -- contribute perturbations below machine precision and should be
                -- deflated to avoid ill-conditioned secular equation roots.
                dMaxAbs <- readRawD wsDSort 0 (nn - 1)
                dMinAbs <- readRawD wsDSort 0 0
                let !matNorm = max (abs dMaxAbs) (abs dMinAbs) + rho * zn2
                    !deflTol = max (tol * sqrt zn2) (8 * 2.220446049250313e-16 * matNorm)
                kND <- deflatePartition wsZSort 0 wsPerm 0 nn deflTol

                -- Extract Qlocal_sub with permuted columns into wsQsub (maxN × nn)
                -- Only rows [loL..loL+nn-1] are relevant, but we copy all maxN rows
                -- to maintain the accumulator's full row structure for the GEMM.
                forM_ [0..nn-1] $ \sortedJ -> do
                  origIdx <- readRawD wsIdx 0 sortedJ
                  let !origJ = round origIdx :: Int
                      !srcCol = loL + origJ
                  forM_ [0..maxN-1] $ \row -> do
                    qv <- readRawD wsQlocal 0 (row * maxN + srcCol)
                    writeRawD wsQsub 0 (row * nn + sortedJ) qv

                if kND == 0
                  then do
                    -- All deflated: eigenvalues = sorted d, eigenvectors = sorted Qlocal cols
                    forM_ [0..nn-1] $ \i ->
                      forM_ [0..maxN-1] $ \row -> do
                        qv <- readRawD wsQsub 0 (row * nn + i)
                        writeRawD wsQlocal 0 (row * maxN + (loL + i)) qv
                    forM_ [0..nn-1] $ \i -> do
                      di <- readRawD wsDSort 0 i
                      writeRawD mbaD offD (lo + i) di

                  else if kND == nn
                    then do
                      -- No deflation: full secular solve + full GEMM
                      secularSolve wsLam 0 wsDSort 0 wsZSort 0 rho nn deflTol
                      dcEigenvectors wsW 0 wsDSort 0 wsZSort 0 wsLam 0 rho nn
                      baQsub <- unsafeFreezeByteArray wsQsub
                      baW    <- unsafeFreezeByteArray wsW
                      forM_ [0..maxN * nn - 1] $ \i -> writeRawD wsResult 0 i 0
                      rawGemmKernel baQsub 0 baW 0 wsResult 0 maxN nn nn
                      forM_ [0..nn-1] $ \i ->
                        forM_ [0..maxN-1] $ \row -> do
                          rv <- readRawD wsResult 0 (row * nn + i)
                          writeRawD wsQlocal 0 (row * maxN + (loL + i)) rv
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

                      -- Copy deflated columns from wsQsub to wsQlocal
                      forM_ [kND..nn-1] $ \j -> do
                        pi_ <- readRawI wsPerm 0 j
                        forM_ [0..maxN-1] $ \row -> do
                          qv <- readRawD wsQsub 0 (row * nn + pi_)
                          writeRawD wsQlocal 0 (row * maxN + (loL + j)) qv

                      -- Extract Q_nd (maxN×kND) from non-deflated columns
                      forM_ [0..kND-1] $ \j -> do
                        pi_ <- readRawI wsPerm 0 j
                        forM_ [0..maxN-1] $ \row -> do
                          qv <- readRawD wsQsub 0 (row * nn + pi_)
                          writeRawD wsResult 0 (row * kND + j) qv

                      -- GEMM: wsQsub(maxN×kND) = Q_nd(maxN×kND) × W_nd(kND×kND)
                      baQnd <- unsafeFreezeByteArray wsResult
                      baW   <- unsafeFreezeByteArray wsW
                      forM_ [0..maxN * kND - 1] $ \i -> writeRawD wsQsub 0 i 0
                      rawGemmKernel baQnd 0 baW 0 wsQsub 0 maxN kND kND

                      -- Copy GEMM result (non-deflated columns) to wsQlocal
                      forM_ [0..kND-1] $ \j ->
                        forM_ [0..maxN-1] $ \row -> do
                          rv <- readRawD wsQsub 0 (row * kND + j)
                          writeRawD wsQlocal 0 (row * maxN + (loL + j)) rv

                      -- Write eigenvalues: non-deflated from wsLam, deflated from wsDSort
                      forM_ [0..kND-1] $ \i -> do
                        lam <- readRawD wsLam 0 i
                        writeRawD mbaD offD (lo + i) lam
                      forM_ [kND..nn-1] $ \j -> do
                        pi_ <- readRawI wsPerm 0 j
                        di <- readRawD wsDSort 0 pi_
                        writeRawD mbaD offD (lo + j) di

          -- Apply a k×k rotation matrix to wsQlocal columns [colOff..colOff+k-1] via GEMM
          -- colOff is in LOCAL coordinates.
          applyRotToQlocal !colOff !k rotMat = do
            -- Extract wsQlocal[:, colOff..colOff+k-1] into wsQsub (maxN × k)
            forM_ [0..k-1] $ \j ->
              forM_ [0..maxN-1] $ \row -> do
                qv <- readRawD wsQlocal 0 (row * maxN + (colOff + j))
                writeRawD wsQsub 0 (row * k + j) qv
            -- O(1) freeze for GEMM inputs
            baQsub <- unsafeFreezeByteArray wsQsub
            baRot  <- unsafeFreezeByteArray rotMat
            -- Zero result
            forM_ [0..maxN * k - 1] $ \i -> writeRawD wsResult 0 i 0
            -- GEMM: result(maxN×k) = Qsub(maxN×k) * Rot(k×k)
            rawGemmKernel baQsub 0 baRot 0 wsResult 0 maxN k k
            -- Copy result back to wsQlocal
            forM_ [0..k-1] $ \j ->
              forM_ [0..maxN-1] $ \row -> do
                rv <- readRawD wsResult 0 (row * k + j)
                writeRawD wsQlocal 0 (row * maxN + (colOff + j)) rv

      -- Run the D&C recursion (operates on wsQlocal and mbaD/mbaE)
      dcGo lo0 hi0

      -- Final step: apply wsQlocal to global Q via GEMM
      -- Q[:, lo0..hi0] = Q[:, lo0..hi0] * wsQlocal
      -- wsQlocal is maxN × maxN; Q is fullN × fullN
      -- We need: for each output column j in [lo0..hi0],
      --   Q_new[row][lo0+j'] = Σ_i Q[row][lo0+i] * wsQlocal[i][j']
      -- where j' = j - lo0, i ranges over [0..maxN-1]
      --
      -- Extract Q[:, lo0..hi0] into wsQsub (fullN × maxN)
      forM_ [0..maxN-1] $ \j ->
        forM_ [0..fullN-1] $ \row -> do
          qv <- readRawD mbaQ 0 (offQ + row * fullN + (lo0 + j))
          writeRawD wsQsub 0 (row * maxN + j) qv
      -- Freeze both for GEMM
      baQsub   <- unsafeFreezeByteArray wsQsub
      baQlocal <- unsafeFreezeByteArray wsQlocal
      -- Zero result
      forM_ [0..fullN * maxN - 1] $ \i -> writeRawD wsResult 0 i 0
      -- GEMM: result(fullN × maxN) = Qsub(fullN × maxN) * Qlocal(maxN × maxN)
      rawGemmKernel baQsub 0 baQlocal 0 wsResult 0 fullN maxN maxN
      -- Copy result back to global Q
      forM_ [0..maxN-1] $ \j ->
        forM_ [0..fullN-1] $ \row -> do
          rv <- readRawD wsResult 0 (row * maxN + j)
          writeRawD mbaQ 0 (offQ + row * fullN + (lo0 + j)) rv

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
    di <- readRawD mbaD offD i
    if abs zi <= deflTol
      then do
        -- Deflated: eigenvalue = d[i]
        writeRawD mbaLam offLam i di
      else do
        -- For small z[i], use first-order perturbation formula directly.
        -- This avoids the iterative solver's difficulty with nearly-flat
        -- secular functions near d[i].
        let !zi2 = zi * zi
            !pertTol = sqrt (2.220446049250313e-16) * (1 + abs di)
        if abs zi < pertTol
          then do
            -- Perturbation: lambda ≈ d[i] + rho * z[i]² / (1 + rho * Σ_{j≠i} z[j]²/(d[j]-d[i]))
            farSum <- farPoleSumSkip mbaD offD mbaZ offZ nn i di
            let !denom = 1 + rho * farSum
                !delta = rho * zi2 / denom
            writeRawD mbaLam offLam i (di + delta)
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
        lam0 <- fixedWeightLoop 0 mid dj dj1 gap zj2 (zj1 * zj1)
        -- Polish with Newton iterations for higher accuracy
        newtonPolish 0 lam0 dj dj1
      else do
        -- Last root when rho > 0: between d[n-1] and d[n-1] + rho*||z||²
        zn2 <- sumZSq mbaZ offZ nn
        -- Compute better initial guess via perturbation theory:
        -- f(d[nn-1]+δ) = 0 ⟹ δ ≈ rho * z[nn-1]² / (1 + rho * Σ_{i<nn-1} z[i]²/(d[i]-d[nn-1]))
        farSum <- farPoleSum mbaD offD mbaZ offZ (nn - 1) dj
        let !denominator = 1 + rho * farSum
            !delta0 = if abs denominator > 1e-300
                      then rho * zj2 / denominator
                      else rho * zn2
            !hi_ = dj + max (rho * zn2) (2 * delta0)
            !mid = dj + max delta0 (1e-15 * (1 + abs dj))
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
        farSum <- farPoleSum mbaD offD mbaZ offZ nn dj
        let !denominator = 1 + rho * farSum
            !delta0 = if abs denominator > 1e-300
                      then abs rho * zj2 / abs denominator
                      else abs rho * zn2
            !lo_ = dj - max (abs rho * zn2) (2 * delta0)
            !mid = dj - max delta0 (1e-15 * (1 + abs dj))
        newtonLoop 0 mid lo_ dj
  where
    !maxIter_ = 100 :: Int

    -- Newton polishing: 3 Newton steps to refine eigenvalue to machine precision.
    -- Uses the secular function and its derivative for rapid convergence.
    newtonPolish !iter !lam !lb !ub
      | iter >= 3 = pure lam
      | otherwise = do
          (f, fp) <- secularFuncAndDeriv mbaD offD mbaZ offZ rho nn lam
          if abs f < 1e-15 * (1 + abs lam) || abs fp < 1e-300
            then pure lam
            else do
              let !step = f / fp
                  !lamNew = lam - step
                  !clamped = max lb (min ub lamNew)
              if abs (clamped - lam) < 1e-16 * (1 + abs lam)
                then pure clamped
                else newtonPolish (iter + 1) clamped lb ub

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
                  let !(lb', ub') = if f * rho > 0 then (lb, lam) else (lam, ub)
                      !lamNew = (lb' + ub') * 0.5
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
              let !(lb', ub') = if f > 0 then (lb, lam) else (lam, ub)
                  !step = f / fp
                  !lamNew0 = lam - step
                  !lamNew = if lamNew0 <= lb' || lamNew0 >= ub'
                            then (lb' + ub') * 0.5
                            else lamNew0
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

-- | Sum of z[i]^2 / (d[i] - dj) for i in [0..skip-1], skipping near-zero denominators.
farPoleSum :: MutableByteArray s -> Int -> MutableByteArray s -> Int
           -> Int -> Double -> ST s Double
farPoleSum mbaD offD mbaZ offZ skip dj = go 0 0
  where
    go i !acc
      | i >= skip = pure acc
      | otherwise = do
          di <- readRawD mbaD offD i
          zi <- readRawD mbaZ offZ i
          let !diff = di - dj
          if abs diff < 1e-300
            then go (i+1) acc
            else go (i+1) (acc + zi * zi / diff)

-- | Sum of z[k]^2 / (d[k] - dj) for all k in [0..nn-1] except k == skip.
farPoleSumSkip :: MutableByteArray s -> Int -> MutableByteArray s -> Int
               -> Int -> Int -> Double -> ST s Double
farPoleSumSkip mbaD offD mbaZ offZ nn skip dj = go 0 0
  where
    go i !acc
      | i >= nn = pure acc
      | i == skip = go (i+1) acc
      | otherwise = do
          di <- readRawD mbaD offD i
          zi <- readRawD mbaZ offZ i
          let !diff = di - dj
          if abs diff < 1e-300
            then go (i+1) acc
            else go (i+1) (acc + zi * zi / diff)

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
