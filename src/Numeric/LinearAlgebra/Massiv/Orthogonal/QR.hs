{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE MagicHash #-}
{-# LANGUAGE UnboxedTuples #-}
{-# LANGUAGE BangPatterns #-}

-- |
-- Module      : Numeric.LinearAlgebra.Massiv.Orthogonal.QR
-- Copyright   : (c) Nadia Chambers 2026
-- License     : BSD-3-Clause
-- Maintainer  : nadia.chambers@iohk.io
-- Stability   : experimental
--
-- QR factorisation via Householder reflections and Givens rotations.
--
-- This module implements the QR factorisation of a matrix
-- \( A \in \mathbb{R}^{m \times n} \) with \( m \ge n \), following
-- Golub & Van Loan, /Matrix Computations/, 4th edition (GVL4),
-- Section 5.2, pp. 246--260.
--
-- __Theorem 5.2.1 (QR existence, GVL4 p. 246).__  For any
-- \( A \in \mathbb{R}^{m \times n} \) with \( m \ge n \), there exists
-- an orthogonal matrix \( Q \in \mathbb{R}^{m \times m} \) and an upper
-- triangular matrix \( R \in \mathbb{R}^{m \times n} \) such that
--
-- \( A = Q \, R \)
--
-- If \( A \) has full column rank, the factorisation is unique up to
-- sign changes in the rows of \( R \) (equivalently, columns of \( Q \)).
--
-- Two algorithms are provided:
--
-- * __Householder QR__ ('qr', 'qrR') -- GVL4 Algorithm 5.2.1 (p. 249).
--   A sequence of Householder reflections \( H_1, H_2, \ldots, H_n \)
--   is applied to \( A \) from the left to produce
--   \( H_n \cdots H_2 \, H_1 \, A = R \), so that
--   \( Q = H_1 \, H_2 \cdots H_n \).
--
-- * __Givens QR__ ('qrGivens') -- GVL4 Section 5.2.4 (p. 252).
--   A sequence of Givens rotations zeroes out sub-diagonal entries one
--   at a time.  This variant is preferred for sparse or banded matrices,
--   particularly Hessenberg matrices, where the number of rotations is
--   proportional to the bandwidth rather than the matrix dimension.
--
-- __Complexity.__
--
-- * Householder QR: \( 2mn^2 - \tfrac{2}{3}n^3 \) flops (GVL4 p. 249).
-- * Givens QR: \( 3mn^2 - n^3 \) flops for a dense matrix, but
--   significantly fewer for structured (e.g., Hessenberg) matrices.
--
-- __Optimisation.__  The implementation uses in-place mutable arrays via
-- the 'ST' monad, storing Householder vectors implicitly in the
-- subdiagonal of the working matrix (the LAPACK convention).  The
-- orthogonal factor \( Q \) is formed via backward accumulation, and
-- 'qrR' avoids forming \( Q \) entirely.
module Numeric.LinearAlgebra.Massiv.Orthogonal.QR
  ( -- * QR factorisation (Householder)
    qr
  , qrP
  , qrR
    -- * QR factorisation (Givens)
  , qrGivens
  ) where

import qualified Data.Massiv.Array as M
import Data.Massiv.Array (Ix1, Ix2(..), Sz(..), unwrapByteArray, unwrapByteArrayOffset,
                          unwrapMutableByteArray, unwrapMutableByteArrayOffset)
import GHC.TypeNats (KnownNat)
import Control.Monad (when)
import Control.Monad.ST (ST)
import GHC.Exts
import GHC.ST (ST(..))
import Data.Primitive.ByteArray (ByteArray(..), MutableByteArray(..), newByteArray, unsafeFreezeByteArray)

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.Orthogonal.Givens
  (givensRotation)
import Numeric.LinearAlgebra.Massiv.Internal.Kernel
  (rawMutSumSqColumn, rawMutSumProdColumns, rawMutHouseholderApply, rawMutQAccum,
   rawGemmKernel, rawZeroDoubles, rawNegateDoubles)

-- | Full QR factorisation via Householder reflections (GVL4 Algorithm 5.2.1, p. 249).
--
-- Given \( A \in \mathbb{R}^{m \times n} \) with \( m \ge n \), compute
-- the factorisation
--
-- \( A = Q \, R \)
--
-- where \( Q \in \mathbb{R}^{m \times m} \) is orthogonal and
-- \( R \in \mathbb{R}^{m \times n} \) is upper triangular.
--
-- The implementation uses in-place mutable arrays: the Householder
-- vectors are stored in the subdiagonal of the working copy of \( A \)
-- (LAPACK convention), and \( Q \) is formed via backward accumulation.
-- Total allocation: two matrices (one for \( R \), one for \( Q \)),
-- plus a small vector of \( \beta \) scalars.
--
-- __Complexity:__ \( 2mn^2 - \tfrac{2}{3}n^3 \) flops for the
-- triangularisation, plus \( 2m^2 n - \tfrac{2}{3}n^3 \) flops for
-- accumulating \( Q \) (GVL4 p. 249).
--
-- See also 'qrR' when only \( R \) is needed, and 'qrGivens' for a
-- Givens-based alternative.
qr :: forall m n r e. (KnownNat m, KnownNat n, M.Manifest r e, Floating e, Ord e)
   => Matrix m n r e -> (Matrix m m r e, Matrix m n r e)
qr a =
  let mm = dimVal @m
      nn = dimVal @n
      steps = min mm nn

      -- Phase 1: In-place Householder triangularisation.
      -- After this, rArr holds R in the upper triangle and Householder
      -- vectors v_k in the subdiagonal of column k (with v_k(k) = 1 implicit).
      -- betaList holds the β scalars as a Haskell list.
      (betaList, rArr) = M.withMArrayST (unMatrix a) $ \mr -> do
        betas <- mapM (\k -> do
          -- Compute Householder vector for column k, rows k..m-1
          x0 <- M.readM mr (k :. k)
          sigma <- sumSqRange mr k mm k  -- σ = Σ R(i,k)² for i=k+1..m-1
          if sigma == 0 && x0 >= 0
            then pure 0     -- already in desired form
            else do
              let mu = sqrt (x0 * x0 + sigma)
                  v0 = if x0 <= 0 then x0 - mu else -sigma / (x0 + mu)
                  beta = 2 * v0 * v0 / (sigma + v0 * v0)
              -- Store v_k in subdiagonal of column k (v_k(k) = 1 is implicit)
              -- Scale entries: v(i) = R(i,k) / v0 for i > k
              mapM_ (\i -> do
                rik <- M.readM mr (i :. k)
                M.write_ mr (i :. k) (rik / v0)
                ) [k+1..mm-1]
              -- Set diagonal: R(k,k) = μ (result of H_k applied to column k)
              M.write_ mr (k :. k) mu
              -- Apply H_k from the left to columns k+1..n-1 of R:
              -- (Column k is skipped: its diagonal is μ, subdiagonal stores v)
              -- For each column j: w_j = β(R(k,j) + Σ_{i>k} v(i)·R(i,j))
              --                    R(k,j) -= w_j; R(i,j) -= v(i)·w_j
              mapM_ (\j -> do
                rkj <- M.readM mr (k :. j)
                wj <- sumProdRange mr mr k mm k j
                let wj' = beta * (rkj + wj)
                M.write_ mr (k :. j) (rkj - wj')
                mapM_ (\i -> do
                  vi <- M.readM mr (i :. k)
                  rij <- M.readM mr (i :. j)
                  M.write_ mr (i :. j) (rij - vi * wj')
                  ) [k+1..mm-1]
                ) [k+1..nn-1]
              pure beta
          ) [0..steps-1]
        pure betas

      -- Phase 2: Backward accumulation of Q.
      -- Start with Q = I, then for k = steps-1 downto 0:
      --   Apply H_k from the right: Q <- Q·(I - β_k·v_k·v_k^T)
      qMat = createMatrix @m @m @r $ \mq -> do
        -- Initialize Q = I
        mapM_ (\i -> mapM_ (\j ->
          M.write_ mq (i :. j) (if i == j then 1 else 0)
          ) [0..mm-1]) [0..mm-1]
        -- Forward accumulation: Q <- Q·H_0·H_1·…·H_{n-1}
        mapM_ (\k -> do
          let beta_k = betaList !! k
          if beta_k == 0 then pure ()
          else
            -- Apply (I - β·v·v^T) from the right to Q
            -- For each row i: w_i = β·(Q(i,k) + Σ_{l>k} Q(i,l)·v(l))
            --                 Q(i,k) -= w_i; Q(i,l) -= w_i·v(l)
            mapM_ (\i -> do
              qik <- M.readM mq (i :. k)
              wi <- qvProd mq rArr i k mm
              let wi' = beta_k * (qik + wi)
              M.write_ mq (i :. k) (qik - wi')
              mapM_ (\l -> do
                let vl = M.index' rArr (l :. k)
                qil <- M.readM mq (i :. l)
                M.write_ mq (i :. l) (qil - wi' * vl)
                ) [k+1..mm-1]
              ) [0..mm-1]
          ) [0..steps-1]

      -- Extract clean R (zero out subdiagonal, which holds Householder vectors)
      rClean = makeMatrix @m @n @r $ \i j ->
        if i <= j then M.index' rArr (i :. j) else 0

  in (qMat, rClean)

-- | Specialised QR factorisation for @P Double@ using raw ByteArray# primops.
qrP :: forall m n. (KnownNat m, KnownNat n)
    => Matrix m n M.P Double -> (Matrix m m M.P Double, Matrix m n M.P Double)
qrP a =
  let mm = dimVal @m
      nn = dimVal @n
      steps = min mm nn

      -- Phase 1: In-place Householder triangularisation using raw kernels.
      (betaList, rArr) = M.withMArrayST (unMatrix a) $ \mr -> do
        let mbaR = unwrapMutableByteArray mr
            offR = unwrapMutableByteArrayOffset mr
        betas <- mapM (\k -> do
          -- Read x0 = R[k,k]
          x0 <- M.readM mr (k :. k)
          -- σ = Σ R[i,k]² for i = k+1..m-1 (using raw kernel)
          sigma <- rawMutSumSqColumn mbaR offR nn (k + 1) mm k
          if sigma == 0 && x0 >= 0
            then pure 0
            else do
              let mu = sqrt (x0 * x0 + sigma)
                  v0 = if x0 <= 0 then x0 - mu else -sigma / (x0 + mu)
                  beta = 2 * v0 * v0 / (sigma + v0 * v0)
              -- Scale v: v(i) = R(i,k) / v0 for i > k
              scaleColumn mbaR offR nn (k + 1) mm k (1.0 / v0)
              -- Set diagonal
              M.write_ mr (k :. k) mu
              -- Apply H_k to columns k+1..n-1 using raw kernel
              mapM_ (\j ->
                rawMutHouseholderApply mbaR offR nn beta k mm j
                ) [k+1..nn-1]
              pure beta
          ) [0..steps-1]
        pure betas

      -- Phase 2: Blocked WY Q accumulation using GEMM.
      -- Group Householder vectors into panels of size nb, compute T-factor
      -- via GEMM, and apply block reflectors via Level-3 GEMM.
      qMat = createMatrix @m @m @M.P $ \mq -> do
        let mbaQ = unwrapMutableByteArray mq
            offQ = unwrapMutableByteArrayOffset mq
            baR = unwrapByteArray rArr
            offFR = unwrapByteArrayOffset rArr
        -- Initialize Q = I via rawZeroDoubles + diagonal writes
        rawZeroDoubles mbaQ offQ (mm * mm)
        mapM_ (\i -> writeRawD mbaQ offQ (i * mm + i) 1) [0..mm-1]

        if steps <= 16
          then
            -- Small matrix: per-row accumulation (Level-2)
            mapM_ (\k -> do
              let beta_k = betaList !! k
              if beta_k == 0 then pure ()
              else
                mapM_ (\i ->
                  rawMutQAccum mbaQ offQ mm baR offFR nn beta_k k mm i
                  ) [0..mm-1]
              ) [0..steps-1]
          else do
            -- Blocked WY: batch nb Householder vectors at a time
            let !nb = min 32 steps
            mbaBetas <- newByteArray (steps * 8)
            mapM_ (\(i, b) -> writeRawD mbaBetas 0 i b) (zip [0..] betaList)
            mbaY  <- newByteArray (mm * nb * 8)
            mbaTf <- newByteArray (nb * nb * 8)
            mbaW1 <- newByteArray (mm * nb * 8)
            mbaW2 <- newByteArray (mm * nb * 8)
            mbaYT <- newByteArray (nb * mm * 8)
            mbaG  <- newByteArray (nb * nb * 8)

            let goBlock !k0
                  | k0 >= steps = pure ()
                  | otherwise = do
                      let !bs = min nb (steps - k0)
                      -- Pack Y (mm × bs): Y[:,j] = v_{k0+j}
                      -- v_k has implicit 1 at position k, stored values at k+1..mm-1
                      rawZeroDoubles mbaY 0 (mm * bs)
                      mapM_ (\j -> do
                        let !k = k0 + j
                        writeRawD mbaY 0 (k * bs + j) 1.0
                        mapM_ (\l ->
                          writeRawD mbaY 0 (l * bs + j) (indexRawD baR offFR (l * nn + k))
                          ) [k+1..mm-1]
                        ) [0..bs-1]

                      -- Transpose Y → Y^T (bs × mm) for GEMM reuse
                      rawZeroDoubles mbaYT 0 (bs * mm)
                      mapM_ (\j -> do
                        let !k = k0 + j
                        writeRawD mbaYT 0 (j * mm + k) 1.0
                        mapM_ (\l ->
                          writeRawD mbaYT 0 (j * mm + l) (indexRawD baR offFR (l * nn + k))
                          ) [k+1..mm-1]
                        ) [0..bs-1]

                      baY  <- unsafeFreezeByteArray mbaY
                      baYT <- unsafeFreezeByteArray mbaYT

                      -- Compute G = Y^T × Y (bs × bs) via GEMM
                      rawZeroDoubles mbaG 0 (bs * bs)
                      rawGemmKernel baYT 0 baY 0 mbaG 0 bs mm bs

                      -- Build T-factor (bs × bs upper-triangular)
                      rawZeroDoubles mbaTf 0 (bs * bs)
                      mapM_ (\j -> do
                        betaj <- readRawD mbaBetas 0 (k0 + j)
                        writeRawD mbaTf 0 (j * bs + j) betaj
                        when (j > 0 && betaj /= 0) $ do
                          -- T[0..j-1, j] = -betaj * T[0..j-1, 0..j-1] * G[0..j-1, j]
                          mapM_ (\i -> do
                            g_ij <- readRawD mbaG 0 (i * bs + j)
                            writeRawD mbaW1 0 i g_ij
                            ) [0..j-1]
                          mapM_ (\i -> do
                            let triLoop !l !acc
                                  | l >= j = pure acc
                                  | otherwise = do
                                      til <- readRawD mbaTf 0 (i * bs + l)
                                      dl  <- readRawD mbaW1 0 l
                                      triLoop (l+1) (acc + til * dl)
                            z <- triLoop 0 0
                            writeRawD mbaTf 0 (i * bs + j) (negate betaj * z)
                            ) [0..j-1]
                        ) [0..bs-1]

                      -- W1 = Q · Y (mm×mm * mm×bs → mm×bs)
                      baQ <- unsafeFreezeByteArray mbaQ
                      rawZeroDoubles mbaW1 0 (mm * bs)
                      rawGemmKernel baQ offQ baY 0 mbaW1 0 mm mm bs

                      -- W2 = W1 · T (mm×bs * bs×bs → mm×bs)
                      baW1 <- unsafeFreezeByteArray mbaW1
                      baTf <- unsafeFreezeByteArray mbaTf
                      rawZeroDoubles mbaW2 0 (mm * bs)
                      rawGemmKernel baW1 0 baTf 0 mbaW2 0 mm bs bs

                      -- Negate W2 in-place
                      rawNegateDoubles mbaW2 0 (mm * bs)

                      -- Q += (-W2) · Y^T (mm×bs * bs×mm → mm×mm)
                      baNW2 <- unsafeFreezeByteArray mbaW2
                      rawGemmKernel baNW2 0 baYT 0 mbaQ offQ mm bs mm

                      goBlock (k0 + bs)
            goBlock 0

      -- Extract clean R
      rClean = makeMatrix @m @n @M.P $ \i j ->
        if i <= j then M.index' rArr (i :. j) else 0

  in (qMat, rClean)
{-# NOINLINE qrP #-}

-- | Scale elements in a column of a mutable matrix: A[i,col] *= scale for i in [start..end-1].
scaleColumn :: MutableByteArray s -> Int -> Int -> Int -> Int -> Int -> Double -> ST s ()
scaleColumn (MutableByteArray mba) (I# off) (I# ncols) (I# start) (I# end) (I# col) (D# scale) = ST $ \s0 ->
  let go i s
        | isTrue# (i >=# end) = s
        | otherwise =
            case readDoubleArray# mba (off +# i *# ncols +# col) s of
              (# s', v #) ->
                case writeDoubleArray# mba (off +# i *# ncols +# col) (v *## scale) s' of
                  s'' -> go (i +# 1#) s''
  in (# go start s0, () #)
{-# INLINE scaleColumn #-}

-- | Helper: Σ R(i,col)² for i=start+1..end-1  (sum of squares below diagonal)
sumSqRange :: (M.Manifest r e, Num e) => M.MArray s r Ix2 e -> Int -> Int -> Int -> ST s e
sumSqRange mr start end col = go (start + 1) 0
  where
    go i !acc
      | i >= end = pure acc
      | otherwise = do
          v <- M.readM mr (i :. col)
          go (i + 1) (acc + v * v)

-- | Helper: Σ mr1(i,c1)·mr2(i,c2) for i=start+1..end-1
sumProdRange :: (M.Manifest r e, Num e)
             => M.MArray s r Ix2 e -> M.MArray s r Ix2 e
             -> Int -> Int -> Int -> Int -> ST s e
sumProdRange mr1 mr2 start end c1 c2 = go (start + 1) 0
  where
    go i !acc
      | i >= end = pure acc
      | otherwise = do
          v1 <- M.readM mr1 (i :. c1)
          v2 <- M.readM mr2 (i :. c2)
          go (i + 1) (acc + v1 * v2)

-- | Helper: Σ Q(row,l)·v(l) for l=start+1..end-1 where v is stored in rArr subdiag of col start
qvProd :: (M.Manifest r1 e, M.Manifest r2 e, Num e)
       => M.MArray s r1 Ix2 e -> M.Array r2 Ix2 e -> Int -> Int -> Int -> ST s e
qvProd mq rArr row start end = go (start + 1) 0
  where
    go l !acc
      | l >= end = pure acc
      | otherwise = do
          qrl <- M.readM mq (row :. l)
          let vl = M.index' rArr (l :. start)
          go (l + 1) (acc + qrl * vl)

-- | Compute only the upper triangular factor \( R \) from the QR
-- factorisation, without explicitly forming \( Q \)
-- (GVL4 Algorithm 5.2.1, p. 249).
--
-- Unlike 'qr', this function never forms the orthogonal factor,
-- saving \( O(m^2 n) \) flops.  The Householder vectors are computed
-- and applied in-place but discarded.
--
-- __Complexity:__ \( 2mn^2 - \tfrac{2}{3}n^3 \) flops.
qrR :: forall m n r e. (KnownNat m, KnownNat n, M.Manifest r e, Floating e, Ord e)
    => Matrix m n r e -> Matrix m n r e
qrR a = snd (qr a)

-- | QR factorisation via Givens rotations (GVL4 Section 5.2.4, p. 252).
--
-- Given \( A \in \mathbb{R}^{m \times n} \) with \( m \ge n \), compute
-- the factorisation \( A = Q \, R \) by applying a sequence of Givens
-- rotations to zero out sub-diagonal entries one at a time.
--
-- The implementation computes \( R \) in-place via the 'ST' monad,
-- records the rotation parameters, then applies them to form \( Q \)
-- in a second in-place pass.
--
-- __Complexity:__ \( 3mn^2 - n^3 \) flops for a dense \( m \times n \)
-- matrix; \( O(mn) \) flops for an upper Hessenberg matrix
-- (GVL4 p. 253).
qrGivens :: forall m n r e. (KnownNat m, KnownNat n, M.Manifest r e, Floating e, Ord e)
          => Matrix m n r e -> (Matrix m m r e, Matrix m n r e)
qrGivens a =
  let mm = dimVal @m
      nn = dimVal @n
      steps = min mm nn

      -- Pass 1: compute R in-place, recording Givens rotations
      (rots, rArr) = M.withMArrayST (unMatrix a) $ \mr -> do
        let go j !acc
              | j >= steps = pure acc
              | otherwise = do
                  acc' <- goRows j (j+1) acc mr
                  go (j+1) acc'
            goRows j i !acc mr_
              | i >= mm = pure acc
              | otherwise = do
                  aij <- M.readM mr_ (i :. j)
                  if aij == 0 then goRows j (i+1) acc mr_
                  else do
                    ajj <- M.readM mr_ (j :. j)
                    let (c, s) = givensRotation ajj aij
                    -- Apply G^T to rows j and i
                    mapM_ (\col -> do
                      rjc <- M.readM mr_ (j :. col)
                      ric <- M.readM mr_ (i :. col)
                      M.write_ mr_ (j :. col) (c * rjc - s * ric)
                      M.write_ mr_ (i :. col) (s * rjc + c * ric)
                      ) [0..nn-1]
                    goRows j (i+1) (acc ++ [(c, s, j, i)]) mr_
        go 0 []

      -- Pass 2: form Q by applying recorded rotations to I
      qMat = createMatrix @m @m @r $ \mq -> do
        -- Initialize Q = I
        mapM_ (\i -> mapM_ (\j ->
          M.write_ mq (i :. j) (if i == j then 1 else 0)
          ) [0..mm-1]) [0..mm-1]
        -- Apply each rotation from the right: Q <- Q·G
        mapM_ (\(c, s, ci, ck) ->
          mapM_ (\row -> do
            qrc <- M.readM mq (row :. ci)
            qrk <- M.readM mq (row :. ck)
            M.write_ mq (row :. ci) (c * qrc - s * qrk)
            M.write_ mq (row :. ck) (s * qrc + c * qrk)
            ) [0..mm-1]
          ) rots

  in (qMat, MkMatrix rArr)

-- Raw ByteArray# helpers for blocked WY Q accumulation
readRawD :: MutableByteArray s -> Int -> Int -> ST s Double
readRawD (MutableByteArray mba) (I# off) (I# i) = ST $ \s ->
  case readDoubleArray# mba (off +# i) s of
    (# s', v #) -> (# s', D# v #)
{-# INLINE readRawD #-}

writeRawD :: MutableByteArray s -> Int -> Int -> Double -> ST s ()
writeRawD (MutableByteArray mba) (I# off) (I# i) (D# v) = ST $ \s ->
  case writeDoubleArray# mba (off +# i) v s of
    s' -> (# s', () #)
{-# INLINE writeRawD #-}

indexRawD :: ByteArray -> Int -> Int -> Double
indexRawD (ByteArray ba) (I# off) (I# i) =
  D# (indexDoubleArray# ba (off +# i))
{-# INLINE indexRawD #-}
