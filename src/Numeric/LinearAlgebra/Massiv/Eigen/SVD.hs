{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE MagicHash #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE UnboxedTuples #-}

-- |
-- Module      : Numeric.LinearAlgebra.Massiv.Eigen.SVD
-- Copyright   : (c) Nadia Chambers 2026
-- License     : BSD-3-Clause
-- Maintainer  : nadia.chambers@iohk.io
-- Stability   : experimental
--
-- Singular Value Decomposition (SVD) of a general real matrix, following
-- Golub & Van Loan, /Matrix Computations/, 4th edition (GVL4), Section 8.6,
-- pp. 498--512.
--
-- __Theorem 8.6.1 (SVD Existence, p. 499):__ For any
-- \(A \in \mathbb{R}^{m \times n}\) with \(m \geq n\) there exist orthogonal
-- matrices \(U \in \mathbb{R}^{m \times m}\) and
-- \(V \in \mathbb{R}^{n \times n}\) such that
--
-- \[
--   A = U \, \Sigma \, V^T, \qquad
--   \Sigma = \mathrm{diag}(\sigma_1, \ldots, \sigma_n),
--   \qquad \sigma_1 \geq \sigma_2 \geq \cdots \geq \sigma_n \geq 0
-- \]
--
-- The \(\sigma_i\) are the /singular values/ of \(A\) and equal the
-- non-negative square roots of the eigenvalues of \(A^T A\).
module Numeric.LinearAlgebra.Massiv.Eigen.SVD
  ( -- * Full SVD
    svd
  , svdP
  , svdAtAP
  , svdGKP
    -- * Singular values only
  , singularValues
  , singularValuesP
  ) where

import qualified Data.Massiv.Array as M
import Data.Massiv.Array (Ix2(..), Sz(..), unwrapByteArray, unwrapByteArrayOffset,
                          unwrapMutableByteArray, unwrapMutableByteArrayOffset)
import Data.Primitive.ByteArray (ByteArray(..), MutableByteArray(..), newByteArray,
                                 unsafeFreezeByteArray)
import GHC.TypeNats (KnownNat)
import Control.Monad (forM_, when)
import Control.Monad.ST (runST)
import Data.List (sortBy)
import Data.Ord ()
import GHC.Exts
import GHC.ST (ST(..))

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.BLAS.Level3 (matMul, matMulP, transpose, matMulAtAP)
-- matvecP no longer needed: U-matrix now computed via single GEMM
import Numeric.LinearAlgebra.Massiv.Eigen.Symmetric
  ( symmetricEigen, symmetricEigenP, symmetricEigenPDC
  -- D&C secular equation infrastructure (reused for bidiagonal SVD)
  , secularSolve, deflatePartition, dcEigenvectors
  , sumZSq
  , readRawI
  )
import Numeric.LinearAlgebra.Massiv.Internal.Kernel
  ( rawMutSumSqColumn, rawMutSumSqRow
  , rawMutHouseholderApply, rawMutHouseholderApplyRow
  , rawMutQAccum
  , rawMutApplyGivensColumns
  , rawMutApplyGivensColumnsCM
  , rawTransposeToColMajor, rawTransposeFromColMajor
  , rawGemmKernel, rawZeroDoubles, rawNegateDoubles
  , rawCopyColumn )

-- | Compute the full Singular Value Decomposition (GVL4 Theorem 8.6.1,
-- p. 499).
--
-- For an \(m \times n\) matrix \(A\) with \(m \geq n\), computes
--
-- \[
--   A = U \, \Sigma \, V^T
-- \]
--
-- where
--
--   * \(U \in \mathbb{R}^{m \times m}\) is orthogonal (columns are the
--     /left singular vectors/),
--   * \(\Sigma = \mathrm{diag}(\sigma_1, \ldots, \sigma_n)\) with
--     \(\sigma_1 \geq \cdots \geq \sigma_n \geq 0\) (the /singular values/),
--   * \(V \in \mathbb{R}^{n \times n}\) is orthogonal (columns are the
--     /right singular vectors/).
--
-- __Method:__ Forms \(A^T A\) and calls
-- 'Numeric.LinearAlgebra.Massiv.Eigen.Symmetric.symmetricEigen' to obtain
-- the eigendecomposition \(A^T A = V \Lambda V^T\).  Singular values are
-- recovered as \(\sigma_i = \sqrt{\max(0, \lambda_i)}\) and left singular
-- vectors as \(u_i = A v_i / \sigma_i\).  For zero singular values the
-- corresponding column of \(U\) is set to the appropriate standard basis
-- vector.
--
-- Returns @(U, sigma, V)@.
svd :: forall m n r e. (KnownNat m, KnownNat n, M.Manifest r e, Floating e, Ord e)
    => Matrix m n r e
    -> (Matrix m m r e, Vector n r e, Matrix n n r e)
svd a =
  let nn = dimVal @n
      at = transpose a
      ata = matMul at a  -- n×n symmetric positive semidefinite
      -- Eigendecomposition of AᵀA
      (eigvalsRaw_, vRaw_) = symmetricEigen ata (30 * nn) 1e-12
      -- Sort eigenvalues descending; build O(1) permutation array
      permBA_ = buildPermArray
                  (map snd $ sortBy (\(a_,_) (b_,_) -> compare (Down a_) (Down b_))
                             [(eigvalsRaw_ !. i, i) | i <- [0..nn-1]])
                  nn
      v = makeMatrix @n @n @r $ \i j -> vRaw_ ! (i, indexPermArray permBA_ j)
      -- Singular values = sqrt of sorted eigenvalues (clamp negatives to 0)
      sigma = makeVector @n @r $ \j ->
        let ev = eigvalsRaw_ !. indexPermArray permBA_ j
        in if ev > 0 then sqrt ev else 0
      -- Compute U: u_i = A·v_i / σ_i
      -- First, build U by computing A·V column by column
      u = makeMatrix @m @m @r $ \i j ->
        if j < nn then
          let sj = sigma !. j
          in if sj > 1e-14
             then -- u_j = (1/σ_j) · Σ_k A(i,k) · V(k,j)
               let av_ = foldl' (\acc k -> acc + (a ! (i, k)) * (v ! (k, j))) 0 [0..nn-1]
               in av_ / sj
             else -- Zero singular value; use arbitrary orthogonal vector
               if i == j then 1 else 0
        else
          -- Extra columns for m > n: extend to full orthogonal basis
          if i == j then 1 else 0
  in (u, sigma, v)

-- | P-specialised full SVD using raw ByteArray# SIMD kernels throughout.
--
-- Wires 'matMulP' (SIMD GEMM), 'symmetricEigenP' (raw primop QR iteration),
-- and 'matvecP' (SIMD matrix–vector product) into the SVD pipeline.
-- | P-specialised full SVD.  Uses the A^T A eigendecomposition path by default
-- as it is currently faster than the Golub-Kahan bidiagonalisation path
-- (svdGKP) at all sizes.  The GK path will become the default once blocked
-- bidiagonalisation is implemented.
svdP :: forall m n. (KnownNat m, KnownNat n)
     => Matrix m n M.P Double
     -> (Matrix m m M.P Double, Vector n M.P Double, Matrix n n M.P Double)
svdP = svdAtAP
{-# NOINLINE svdP #-}

-- | SVD via A^T A eigendecomposition.
-- Forms A^T A, eigendecomposes via 'symmetricEigenP', recovers singular
-- values as square roots and left singular vectors via matrix-vector products.
svdAtAP :: forall m n. (KnownNat m, KnownNat n)
        => Matrix m n M.P Double
        -> (Matrix m m M.P Double, Vector n M.P Double, Matrix n n M.P Double)
svdAtAP a =
  let !mm = dimVal @m
      !nn = dimVal @n
      !ata = matMulAtAP a  -- n×n symmetric positive semidefinite, fast transpose + SIMD GEMM
      -- Eigendecomposition of AᵀA: D&C for large, QR iteration for small
      (!eigvalsRaw, !vRaw) = if nn >= 50
        then symmetricEigenPDC ata 1e-12
        else symmetricEigenP ata (max 30 (6 * nn)) 1e-12
      -- Sort eigenvalues descending; build O(1) permutation array
      !permList = map snd $ sortBy (\(a_,_) (b_,_) -> compare (Down a_) (Down b_))
                        [(eigvalsRaw !. i, i) | i <- [0..nn-1]]
      !permBA = buildPermArray permList nn
      -- Rearrange V columns via O(1) indexed permutation using raw ByteArray copy
      !v = createMatrix @n @n @M.P $ \mv -> do
        let !baV   = unwrapByteArray (unMatrix vRaw)
            !offV  = unwrapByteArrayOffset (unMatrix vRaw)
            !mbaVP = unwrapMutableByteArray mv
            !offVP = unwrapMutableByteArrayOffset mv
            !(ByteArray baV#) = baV
            !(I# offV#) = offV
            !(MutableByteArray mbaVP#) = mbaVP
            !(I# offVP#) = offVP
            !(I# nnV#) = nn
        -- Copy columns: V_new[i,j] = V_raw[i, perm[j]]
        ST $ \s0 ->
          let goRow i s
                | isTrue# (i >=# nnV#) = s
                | otherwise =
                    let goCol j s1
                          | isTrue# (j >=# nnV#) = s1
                          | otherwise =
                              let !(I# pj) = indexPermArray permBA (I# j)
                                  !val = indexDoubleArray# baV# (offV# +# i *# nnV# +# pj)
                              in case writeDoubleArray# mbaVP# (offVP# +# i *# nnV# +# j) val s1 of
                                   s2 -> goCol (j +# 1#) s2
                    in goRow (i +# 1#) (goCol 0# s)
          in (# goRow 0# s0, () #)
      -- Singular values = sqrt of sorted eigenvalues (clamp negatives to 0)
      sigma = makeVector @n @M.P $ \j ->
        let !(ByteArray baEV#) = unwrapByteArray (unVector eigvalsRaw)
            !(I# offEV#) = unwrapByteArrayOffset (unVector eigvalsRaw)
            !(I# pj#) = indexPermArray permBA j
            ev = case indexDoubleArray# baEV# (offEV# +# pj#) of v_ -> D# v_
        in if ev > 0 then sqrt ev else 0
      -- Compute U = A·V via single GEMM, then scale columns by 1/σ_j
      !av = matMulP a v  -- m×n result via SIMD GEMM
      u = createMatrix @m @m @M.P $ \mu -> do
        let !baAV  = unwrapByteArray (unMatrix av)
            !offAV = unwrapByteArrayOffset (unMatrix av)
            !mbaU  = unwrapMutableByteArray mu
            !offU  = unwrapMutableByteArrayOffset mu
            !(MutableByteArray mbaU#) = mbaU
            !(I# offU#) = offU
            !(I# mm#) = mm
            !(I# nn#) = nn
            !nn4 = nn - (nn `rem` 4)
            !(I# nn4#) = nn4
            !(I# mmxmm#) = mm * mm
        -- Zero-initialise only the padding region U[i, nn..mm-1] when mm > nn.
        -- The SIMD loop below fills all U[i, 0..nn-1], so no zeroing needed there.
        when (mm > nn) $
          ST $ \s0 ->
            let goZ i s
                  | isTrue# (i >=# mm#) = s
                  | otherwise =
                      let goCol j s1
                            | isTrue# (j >=# mm#) = s1
                            | otherwise = case writeDoubleArray# mbaU# (offU# +# i *# mm# +# j) 0.0## s1 of
                                            s2 -> goCol (j +# 1#) s2
                      in goZ (i +# 1#) (goCol nn# s)
            in (# goZ 0# s0, () #)
        -- Pre-compute invSigma vector, then freeze for immutable SIMD reads
        mbaInvS <- newByteArray (nn * 8)
        forM_ [0..nn-1] $ \j -> do
          let sj = sigma !. j
          writeRawD mbaInvS 0 j (if sj > 1e-14 then 1.0 / sj else 0.0)
        !(ByteArray baInvS#) <- unsafeFreezeByteArray mbaInvS
        -- Row-oriented SIMD column-scaling: U[i,j] = invSigma[j] * AV[i,j]
        -- Both AV[i,0..nn-1] and U[i,0..nn-1] are contiguous in row-major layout
        let !(ByteArray baAV#) = baAV
            !(I# offAV#) = offAV
        -- SIMD row loop (immutable reads for AV and invSigma, mutable writes for U)
        ST $ \s0 ->
          let goRow i s
                | isTrue# (i >=# mm#) = s
                | otherwise =
                    let !avRowOff = offAV# +# i *# nn#
                        !uRowOff  = offU# +# i *# mm#
                        -- SIMD phase: process 4 columns at a time
                        goSimd j s1
                          | isTrue# (j >=# nn4#) = s1
                          | otherwise =
                              let avV  = indexDoubleArrayAsDoubleX4# baAV# (avRowOff +# j)
                                  isV  = indexDoubleArrayAsDoubleX4# baInvS# j
                                  !prod = timesDoubleX4# avV isV
                              in case writeDoubleArrayAsDoubleX4# mbaU# (uRowOff +# j) prod s1 of
                                   s2 -> goSimd (j +# 4#) s2
                        -- Scalar cleanup for remainder columns
                        goScalar j s1
                          | isTrue# (j >=# nn#) = s1
                          | otherwise =
                              let avVal = indexDoubleArray# baAV# (avRowOff +# j)
                                  isVal = indexDoubleArray# baInvS# j
                              in case writeDoubleArray# mbaU# (uRowOff +# j) (avVal *## isVal) s1 of
                                   s2 -> goScalar (j +# 1#) s2
                    in goRow (i +# 1#) (goScalar nn4# (goSimd 0# s))
          in (# goRow 0# s0, () #)
        -- Fix zero singular values: set diagonal U[j,j] = 1.0
        forM_ [0..nn-1] $ \j -> do
          let sj = sigma !. j
          when (sj <= 1e-14) $
            writeRawD mbaU offU (j * mm + j) 1.0
        -- Extra columns for m > n: extend to full orthogonal basis
        forM_ [nn..mm-1] $ \j ->
          writeRawD mbaU offU (j * mm + j) 1.0
  in (u, sigma, v)
{-# NOINLINE svdAtAP #-}

-- | Read a Double from an immutable ByteArray at element index.
readBA :: ByteArray -> Int -> Int -> Double
readBA (ByteArray ba) (I# off) (I# i) =
  case indexDoubleArray# ba (off +# i) of v -> D# v
{-# INLINE readBA #-}

-- | Build an unboxed Int permutation array from a list for O(1) indexed access.
buildPermArray :: [Int] -> Int -> ByteArray
buildPermArray xs n = runST $ do
  mba <- newByteArray (n * 8)  -- 8 bytes per Int on 64-bit
  let go _ []     = pure ()
      go i (x:rest) = do
        let !(MutableByteArray mba#) = mba
            !(I# i#) = i
            !(I# x#) = x
        ST $ \s -> case writeIntArray# mba# i# x# s of s' -> (# s', () #)
        go (i + 1) rest
  go 0 xs
  unsafeFreezeByteArray mba
{-# INLINE buildPermArray #-}

-- | O(1) index into a permutation ByteArray.
indexPermArray :: ByteArray -> Int -> Int
indexPermArray (ByteArray ba#) (I# i#) =
  case indexIntArray# ba# i# of x# -> I# x#
{-# INLINE indexPermArray #-}

-- | Compute only the singular values of \(A\), sorted in descending order.
singularValues :: forall m n r e. (KnownNat m, KnownNat n, M.Manifest r e, Floating e, Ord e)
               => Matrix m n r e -> Vector n r e
singularValues a =
  let nn = dimVal @n
      at = transpose a
      ata = matMul at a
      (eigvals, _) = symmetricEigen ata (30 * nn) 1e-12
      -- Sort eigenvalues descending, take sqrt
      evList = map (\i -> eigvals !. i) [0..nn-1]
      sorted = sortBy (\x y -> compare (Down x) (Down y)) evList
  in makeVector @n @r $ \i ->
    let ev = sorted !! i
    in if ev > 0 then sqrt ev else 0

-- | P-specialised singular values using raw SIMD GEMM and raw primop eigenvalue solver.
singularValuesP :: forall m n. (KnownNat m, KnownNat n)
                => Matrix m n M.P Double -> Vector n M.P Double
singularValuesP a =
  let nn = dimVal @n
      !ata = matMulAtAP a
      (!eigvals, _) = symmetricEigenP ata (10 * nn) 1e-12
      evList = map (\i -> eigvals !. i) [0..nn-1]
      sorted = sortBy (\x y -> compare (Down x) (Down y)) evList
  in makeVector @n @M.P $ \i ->
    let ev = sorted !! i
    in if ev > 0 then sqrt ev else 0

-- ============================================================================
-- Golub-Kahan bidiagonalisation SVD (GVL4 Algorithm 5.4.2 + 8.6.2)
-- ============================================================================

-- | Full Golub-Kahan SVD pipeline.
-- Phase 1: Bidiagonalise A → U₀ B V₀^T
-- Phase 2: Implicit-shift QR on bidiagonal B, accumulating rotations into U, V
-- Phase 3: Assemble final U, sigma, V; ensure σᵢ ≥ 0; sort descending
svdGKP :: forall m n. (KnownNat m, KnownNat n)
       => Matrix m n M.P Double
       -> (Matrix m m M.P Double, Vector n M.P Double, Matrix n n M.P Double)
svdGKP (MkMatrix a_) = runST $ do
  let !mm = dimVal @m
      !nn = dimVal @n

  -- Copy input into mutable working storage
  mA <- M.thawS a_
  let mbaA = unwrapMutableByteArray mA
      offA = unwrapMutableByteArrayOffset mA

  -- Allocate arrays for Householder betas
  mbaBetaL <- newByteArray (nn * 8)  -- left Householder betas
  mbaBetaR <- newByteArray (nn * 8)  -- right Householder betas

  -- Phase 1: Bidiagonalise A in-place
  bidiagonalizeP mbaA offA mm nn mbaBetaL mbaBetaR

  -- Extract diagonal d and superdiagonal e from bidiagonalised A
  mbaD <- newByteArray (nn * 8)
  mbaE <- newByteArray (nn * 8)
  forM_ [0..nn-1] $ \k -> do
    dk <- readRawD mbaA offA (k * nn + k)
    writeRawD mbaD 0 k dk
  forM_ [0..nn-2] $ \k -> do
    ek <- readRawD mbaA offA (k * nn + (k+1))
    writeRawD mbaE 0 k ek

  -- Freeze A for Householder vector extraction
  frozenA <- M.freezeS mA

  -- Phase 2: Accumulate U₀ and V₀ from stored Householder vectors
  -- U₀ = H₀ H₁ ... H_{n-1} (left reflectors, stored in columns of A)
  mU <- M.newMArray @M.P (Sz (mm :. mm)) (0 :: Double)
  let mbaU = unwrapMutableByteArray mU
      offU = unwrapMutableByteArrayOffset mU
  -- Initialise U = I
  forM_ [0..mm-1] $ \i ->
    writeRawD mbaU offU (i * mm + i) 1.0

  let baA = unwrapByteArray frozenA
      offFA = unwrapByteArrayOffset frozenA

  -- Accumulate left Householder reflectors into U (forward: U = H₀ H₁ ⋯ H_{n-1})
  -- Left reflector k: v stored in column k of A, rows k+1..m-1, with v[k]=1 implicit
  if nn <= 16
    then
      -- Small matrix: per-row accumulation (Level-2)
      forM_ [0..nn-1] $ \k -> do
        betaK <- readRawD mbaBetaL 0 k
        when (betaK /= 0) $
          forM_ [0..mm-1] $ \row ->
            rawMutQAccum mbaU offU mm baA offFA nn betaK k mm row
    else do
      -- Blocked WY: batch nb Householder vectors at a time
      let !nbU = min 32 nn
      mbaYU  <- newByteArray (mm * nbU * 8)
      mbaTfU <- newByteArray (nbU * nbU * 8)
      mbaW1U <- newByteArray (mm * nbU * 8)
      mbaW2U <- newByteArray (mm * nbU * 8)
      mbaYTU <- newByteArray (nbU * mm * 8)
      mbaGU  <- newByteArray (nbU * nbU * 8)

      let goBlockU !k0
            | k0 >= nn = pure ()
            | otherwise = do
                let !bsz = min nbU (nn - k0)
                -- Pack Y (mm × bsz): Y[:,j] = left Householder vector k0+j
                rawZeroDoubles mbaYU 0 (mm * bsz)
                forM_ [0..bsz-1] $ \j -> do
                  let !k = k0 + j
                  writeRawD mbaYU 0 (k * bsz + j) 1.0
                  forM_ [k+1..mm-1] $ \l ->
                    writeRawD mbaYU 0 (l * bsz + j) (readBA baA offFA (l * nn + k))

                -- Transpose Y → Y^T (bsz × mm) for GEMM reuse
                rawZeroDoubles mbaYTU 0 (bsz * mm)
                forM_ [0..bsz-1] $ \j -> do
                  let !k = k0 + j
                  writeRawD mbaYTU 0 (j * mm + k) 1.0
                  forM_ [k+1..mm-1] $ \l ->
                    writeRawD mbaYTU 0 (j * mm + l) (readBA baA offFA (l * nn + k))

                baYU  <- unsafeFreezeByteArray mbaYU
                baYTU <- unsafeFreezeByteArray mbaYTU

                -- G = Y^T × Y (bsz × bsz)
                rawZeroDoubles mbaGU 0 (bsz * bsz)
                rawGemmKernel baYTU 0 baYU 0 mbaGU 0 bsz mm bsz

                -- Build T-factor (bsz × bsz upper-triangular)
                rawZeroDoubles mbaTfU 0 (bsz * bsz)
                forM_ [0..bsz-1] $ \j -> do
                  betaj <- readRawD mbaBetaL 0 (k0 + j)
                  writeRawD mbaTfU 0 (j * bsz + j) betaj
                  when (j > 0 && betaj /= 0) $ do
                    -- T[0..j-1, j] = -betaj * T[0..j-1, 0..j-1] * G[0..j-1, j]
                    forM_ [0..j-1] $ \i -> do
                      g_ij <- readRawD mbaGU 0 (i * bsz + j)
                      writeRawD mbaW1U 0 i g_ij
                    forM_ [0..j-1] $ \i -> do
                      let triLoop !l !acc
                            | l >= j = pure acc
                            | otherwise = do
                                til <- readRawD mbaTfU 0 (i * bsz + l)
                                dl  <- readRawD mbaW1U 0 l
                                triLoop (l+1) (acc + til * dl)
                      z <- triLoop 0 0
                      writeRawD mbaTfU 0 (i * bsz + j) (negate betaj * z)

                -- W1 = Q · Y (mm×mm * mm×bsz → mm×bsz)
                baQU <- unsafeFreezeByteArray mbaU
                rawZeroDoubles mbaW1U 0 (mm * bsz)
                rawGemmKernel baQU offU baYU 0 mbaW1U 0 mm mm bsz

                -- W2 = W1 · T (mm×bsz * bsz×bsz → mm×bsz)
                baW1U <- unsafeFreezeByteArray mbaW1U
                baTfU <- unsafeFreezeByteArray mbaTfU
                rawZeroDoubles mbaW2U 0 (mm * bsz)
                rawGemmKernel baW1U 0 baTfU 0 mbaW2U 0 mm bsz bsz

                -- Negate W2
                rawNegateDoubles mbaW2U 0 (mm * bsz)

                -- Q += (-W2) · Y^T (mm×bsz * bsz×mm → mm×mm)
                baNW2U <- unsafeFreezeByteArray mbaW2U
                rawGemmKernel baNW2U 0 baYTU 0 mbaU offU mm bsz mm

                goBlockU (k0 + bsz)
      goBlockU 0

  -- V₀ = G₁ G₂ ... G_{n-3} (right reflectors)
  mV <- M.newMArray @M.P (Sz (nn :. nn)) (0 :: Double)
  let mbaV = unwrapMutableByteArray mV
      offV = unwrapMutableByteArrayOffset mV
  -- Initialise V = I
  forM_ [0..nn-1] $ \i ->
    writeRawD mbaV offV (i * nn + i) 1.0

  -- Accumulate right Householder reflectors into V (forward: V = G₀ G₁ ⋯ G_{n-3})
  -- Right reflector k: v stored in row k of A, cols k+2..n-1, with v[k+1]=1 implicit
  if nn < 19
    then
      -- Small: per-row Level-2
      when (nn >= 3) $
        forM_ [0..nn-3] $ \k -> do
          betaK <- readRawD mbaBetaR 0 k
          when (betaK /= 0) $
            forM_ [0..nn-1] $ \row ->
              rightQAccum mbaV offV nn baA offFA nn betaK k nn row
    else when (nn >= 3) $ do
      -- Blocked WY for right reflectors
      let !nRefl = nn - 2  -- right reflectors 0..nn-3
          !nbV = min 32 nRefl
      mbaYV  <- newByteArray (nn * nbV * 8)
      mbaTfV <- newByteArray (nbV * nbV * 8)
      mbaW1V <- newByteArray (nn * nbV * 8)
      mbaW2V <- newByteArray (nn * nbV * 8)
      mbaYTV <- newByteArray (nbV * nn * 8)
      mbaGV  <- newByteArray (nbV * nbV * 8)

      let goBlockV !k0
            | k0 >= nRefl = pure ()
            | otherwise = do
                let !bsz = min nbV (nRefl - k0)
                -- Pack Y (nn × bsz): Y[:,j] = right Householder vector k0+j
                -- Right vector k has implicit 1 at position k+1, stored values at k+2..nn-1
                rawZeroDoubles mbaYV 0 (nn * bsz)
                forM_ [0..bsz-1] $ \j -> do
                  let !k = k0 + j
                  writeRawD mbaYV 0 ((k+1) * bsz + j) 1.0
                  forM_ [k+2..nn-1] $ \l ->
                    writeRawD mbaYV 0 (l * bsz + j) (readBA baA offFA (k * nn + l))

                -- Transpose Y → Y^T (bsz × nn)
                rawZeroDoubles mbaYTV 0 (bsz * nn)
                forM_ [0..bsz-1] $ \j -> do
                  let !k = k0 + j
                  writeRawD mbaYTV 0 (j * nn + (k+1)) 1.0
                  forM_ [k+2..nn-1] $ \l ->
                    writeRawD mbaYTV 0 (j * nn + l) (readBA baA offFA (k * nn + l))

                baYV  <- unsafeFreezeByteArray mbaYV
                baYTV <- unsafeFreezeByteArray mbaYTV

                -- G = Y^T × Y (bsz × bsz)
                rawZeroDoubles mbaGV 0 (bsz * bsz)
                rawGemmKernel baYTV 0 baYV 0 mbaGV 0 bsz nn bsz

                -- Build T-factor (bsz × bsz upper-triangular)
                rawZeroDoubles mbaTfV 0 (bsz * bsz)
                forM_ [0..bsz-1] $ \j -> do
                  betaj <- readRawD mbaBetaR 0 (k0 + j)
                  writeRawD mbaTfV 0 (j * bsz + j) betaj
                  when (j > 0 && betaj /= 0) $ do
                    -- T[0..j-1, j] = -betaj * T[0..j-1, 0..j-1] * G[0..j-1, j]
                    forM_ [0..j-1] $ \i -> do
                      g_ij <- readRawD mbaGV 0 (i * bsz + j)
                      writeRawD mbaW1V 0 i g_ij
                    forM_ [0..j-1] $ \i -> do
                      let triLoop !l !acc
                            | l >= j = pure acc
                            | otherwise = do
                                til <- readRawD mbaTfV 0 (i * bsz + l)
                                dl  <- readRawD mbaW1V 0 l
                                triLoop (l+1) (acc + til * dl)
                      z <- triLoop 0 0
                      writeRawD mbaTfV 0 (i * bsz + j) (negate betaj * z)

                -- W1 = V · Y (nn×nn * nn×bsz → nn×bsz)
                baQV <- unsafeFreezeByteArray mbaV
                rawZeroDoubles mbaW1V 0 (nn * bsz)
                rawGemmKernel baQV offV baYV 0 mbaW1V 0 nn nn bsz

                -- W2 = W1 · T (nn×bsz * bsz×bsz → nn×bsz)
                baW1V <- unsafeFreezeByteArray mbaW1V
                baTfV <- unsafeFreezeByteArray mbaTfV
                rawZeroDoubles mbaW2V 0 (nn * bsz)
                rawGemmKernel baW1V 0 baTfV 0 mbaW2V 0 nn bsz bsz

                -- Negate W2
                rawNegateDoubles mbaW2V 0 (nn * bsz)

                -- V += (-W2) · Y^T (nn×bsz * bsz×nn → nn×nn)
                baNW2V <- unsafeFreezeByteArray mbaW2V
                rawGemmKernel baNW2V 0 baYTV 0 mbaV offV nn bsz nn

                goBlockV (k0 + bsz)
      goBlockV 0

  -- Phase 3: Bidiagonal SVD (D&C for large, QR iteration for small)
  if nn >= dcBidiagThreshold
    then dcBidiagSVD mbaD 0 mbaE 0 mbaU offU mm mbaV offV nn nn 1e-14
    else bidiagQRIterPCM mbaD 0 mbaE 0 mbaU offU mm mbaV offV nn nn (30 * nn)

  -- Phase 4: Ensure σᵢ ≥ 0 (flip sign of U column if needed)
  forM_ [0..nn-1] $ \k -> do
    dk <- readRawD mbaD 0 k
    when (dk < 0) $ do
      writeRawD mbaD 0 k (negate dk)
      -- Flip column k of U
      forM_ [0..mm-1] $ \i -> do
        uik <- readRawD mbaU offU (i * mm + k)
        writeRawD mbaU offU (i * mm + k) (negate uik)

  -- Phase 5: Sort singular values descending and permute U, V columns
  pairs <- mapM (\k -> do dk <- readRawD mbaD 0 k; return (dk, k)) [0..nn-1]
  let !sorted = sortBy (\(a1,_) (b1,_) -> compare (Down a1) (Down b1)) pairs

  frozenU <- M.freezeS mU
  frozenV <- M.freezeS mV
  let baU = unwrapByteArray frozenU
      offFU = unwrapByteArrayOffset frozenU
      baV = unwrapByteArray frozenV
      offFV = unwrapByteArrayOffset frozenV

  let !sigmaVec = makeVector @n @M.P $ \i -> fst (sorted !! i)
      !uMat = makeMatrix @m @m @M.P $ \i j ->
        if j < nn
          then let origCol = snd (sorted !! j)
               in readBA baU offFU (i * mm + origCol)
          else if i == j then 1 else 0
      !vMat = makeMatrix @n @n @M.P $ \i j ->
        let origCol = snd (sorted !! j)
        in readBA baV offFV (i * nn + origCol)

  return (uMat, sigmaVec, vMat)
{-# NOINLINE svdGKP #-}

-- | In-place bidiagonalisation of an m×n matrix stored in a MutableByteArray.
-- GVL4 Algorithm 5.4.2, p. 284.
--
-- After this, the matrix has:
-- - Diagonal d[k] = A[k,k]
-- - Superdiagonal e[k] = A[k,k+1]
-- - Left Householder vectors stored in column k below diagonal (rows k+1..m-1)
-- - Right Householder vectors stored in row k right of superdiag (cols k+2..n-1)
-- - Householder betas stored in mbaBetaL and mbaBetaR
bidiagonalizeP :: MutableByteArray s -> Int -> Int -> Int
               -> MutableByteArray s -> MutableByteArray s -> ST s ()
bidiagonalizeP mbaA offA mm nn mbaBetaL mbaBetaR = do
  forM_ [0..nn-1] $ \k -> do
    -- Left Householder: zero out A[k+1:m, k]
    -- Compute Householder vector for column k, rows k..m-1
    if k < mm - 1
      then do
        -- sigma = Σ A[i,k]² for i in k+1..m-1
        sigma <- rawMutSumSqColumn mbaA offA nn (k+1) mm k
        x0 <- readRawD mbaA offA (k * nn + k)
        if sigma < 1e-300
          then writeRawD mbaBetaL 0 k 0
          else do
            let mu = sqrt (x0 * x0 + sigma)
                v0 = if x0 <= 0 then x0 - mu else -sigma / (x0 + mu)
                beta = 2 * v0 * v0 / (sigma + v0 * v0)
            -- Store v: normalise by v0
            -- v[k] will become 1 (implicit), v[k+1..m-1] = A[i,k]/v0
            forM_ [k+1..mm-1] $ \i -> do
              aik <- readRawD mbaA offA (i * nn + k)
              writeRawD mbaA offA (i * nn + k) (aik / v0)
            -- Set A[k,k] = mu (the diagonal value after reflection)
            writeRawD mbaA offA (k * nn + k) mu
            writeRawD mbaBetaL 0 k beta
            -- Apply left Householder to columns k+1..n-1
            -- Using rawMutHouseholderApply which reads v from column k, rows k+1..m-1
            forM_ [k+1..nn-1] $ \col ->
              rawMutHouseholderApply mbaA offA nn beta k mm col
      else
        writeRawD mbaBetaL 0 k 0

    -- Right Householder: zero out A[k, k+2:n]
    if k < nn - 2
      then do
        -- sigma = Σ A[k,j]² for j in k+2..n-1
        sigma <- rawMutSumSqRow mbaA offA nn k (k+2) nn
        x0 <- readRawD mbaA offA (k * nn + (k+1))
        if sigma < 1e-300
          then writeRawD mbaBetaR 0 k 0
          else do
            let mu = sqrt (x0 * x0 + sigma)
                v0 = if x0 <= 0 then x0 - mu else -sigma / (x0 + mu)
                beta = 2 * v0 * v0 / (sigma + v0 * v0)
            -- Store v: normalise by v0
            -- v[k+1] will become 1 (implicit), v[k+2..n-1] = A[k,j]/v0
            forM_ [k+2..nn-1] $ \j -> do
              akj <- readRawD mbaA offA (k * nn + j)
              writeRawD mbaA offA (k * nn + j) (akj / v0)
            -- Set A[k,k+1] = mu (the superdiagonal value)
            writeRawD mbaA offA (k * nn + (k+1)) mu
            writeRawD mbaBetaR 0 k beta
            -- Apply right Householder to rows k+1..m-1
            -- v is stored in row k, cols k+2..n-1, with implicit v[k+1]=1
            forM_ [k+1..mm-1] $ \row ->
              rawMutHouseholderApplyRow mbaA offA nn beta k (k+1) nn row
      else
        when (k < nn - 1) $ writeRawD mbaBetaR 0 k 0
{-# NOINLINE bidiagonalizeP #-}

-- | Implicit-shift bidiagonal QR iteration (GVL4 Algorithm 8.6.2).
-- Operates on diagonal d and superdiagonal e of an upper bidiagonal matrix.
-- Accumulates left rotations into U (m×n columns) and right rotations into V (n×n).
--
-- Each iteration: (1) find the active unreduced block [p..q] by scanning from
-- the bottom for deflation, then scanning up for split; (2) apply one QR step
-- to [p..q]; (3) repeat until fully deflated or maxIter reached.
bidiagQRIterP :: MutableByteArray s -> Int  -- d + offset
              -> MutableByteArray s -> Int  -- e + offset
              -> MutableByteArray s -> Int -> Int  -- U + offset + ucols
              -> MutableByteArray s -> Int -> Int  -- V + offset + vcols
              -> Int -> Int  -- n, maxIter
              -> ST s ()
bidiagQRIterP mbaD offD mbaE offE mbaU offU ucols mbaV offV vcols nn maxIter = go 0
  where
    go !iter
      | iter >= maxIter = return ()
      | otherwise = do
          -- Step 1: Find q — the bottom of the unreduced block.
          -- Scan from nn-1 downward, deflating negligible e[q-1].
          q <- deflateHi (nn - 1)
          if q <= 0
            then return ()  -- fully deflated
            else do
              -- Step 2: Find p — the top of the unreduced block.
              -- Scan from q-1 downward, looking for a split.
              p <- findLo (q - 1)
              -- Step 3: Apply one QR step to [p..q]
              bidiagQRStep mbaD offD mbaE offE mbaU offU ucols mbaV offV vcols p q
              go (iter + 1)

    -- Scan from hi down: deflate any trailing negligible superdiagonals.
    -- Returns the index of the bottom row of the active block (0 if fully deflated).
    deflateHi !hi
      | hi <= 0 = return 0
      | otherwise = do
          ehi <- readRawD mbaE offE (hi - 1)
          dhi <- readRawD mbaD offD hi
          dhi1 <- readRawD mbaD offD (hi - 1)
          let tol = 1e-14 * (abs dhi1 + abs dhi)
          if abs ehi <= tol
            then do
              writeRawD mbaE offE (hi - 1) 0
              deflateHi (hi - 1)
            else return hi

    -- Scan from idx downward to find the top of the unreduced block.
    -- Returns the smallest p such that B[p..q] is unreduced.
    findLo !idx
      | idx <= 0 = return 0
      | otherwise = do
          eidx <- readRawD mbaE offE (idx - 1)
          didx <- readRawD mbaD offD idx
          didx1 <- readRawD mbaD offD (idx - 1)
          let tol = 1e-14 * (abs didx1 + abs didx)
          if abs eidx <= tol
            then do
              writeRawD mbaE offE (idx - 1) 0
              return idx
            else findLo (idx - 1)
{-# NOINLINE bidiagQRIterP #-}

-- | Column-major bidiagonal QR iteration with AED and stall detection.
-- Transposes U (mm×mm) and V (nn×nn) to column-major layout for SIMD Givens,
-- runs bidiag QR with aggressive early deflation, then transposes back.
-- Falls back to row-major path for nn < 10 (transpose overhead dominates).
bidiagQRIterPCM :: MutableByteArray s -> Int  -- d + offset
                -> MutableByteArray s -> Int  -- e + offset
                -> MutableByteArray s -> Int -> Int  -- U + offset + ucols (= mm)
                -> MutableByteArray s -> Int -> Int  -- V + offset + vcols (= nn)
                -> Int -> Int  -- nn, maxIter
                -> ST s ()
bidiagQRIterPCM mbaD offD mbaE offE mbaU offU mm mbaV offV nn n maxIter
  | n < 10 = bidiagQRIterP mbaD offD mbaE offE mbaU offU mm mbaV offV nn n maxIter
  | otherwise = do
      -- Transpose U (mm×mm) and V (nn×nn) to column-major
      tmpU <- newByteArray (mm * mm * 8)
      tmpV <- newByteArray (nn * nn * 8)
      rawTransposeToColMajor mbaU offU tmpU 0 mm
      rawTransposeToColMajor mbaV offV tmpV 0 nn
      -- Run CM iteration
      goCM 0 (n - 1) 0
        mbaD offD mbaE offE tmpU 0 mm tmpV 0 nn n maxIter
      -- Transpose back to row-major
      rawTransposeFromColMajor tmpU 0 mbaU offU mm
      rawTransposeFromColMajor tmpV 0 mbaV offV nn
{-# NOINLINE bidiagQRIterPCM #-}

-- | CM iteration core with AED and stall detection.
-- Parameters: iter, lastQ, stallCount, then the usual d/e/U/V arrays + n + maxIter.
goCM :: Int -> Int -> Int
     -> MutableByteArray s -> Int -> MutableByteArray s -> Int
     -> MutableByteArray s -> Int -> Int
     -> MutableByteArray s -> Int -> Int
     -> Int -> Int -> ST s ()
goCM !iter !lastQ !stall mbaD offD mbaE offE mbaU offU mm mbaV offV nn n maxIter
  | iter >= maxIter = return ()
  | stall >= 20     = return ()  -- stall detection: bail after 20 steps without deflation
  | otherwise = do
      -- AED: scan bottom w superdiagonal entries for aggressive early deflation
      let w = min 6 ((n + 2) `div` 3)
      aedScan (n - 1) w
      -- Find q — bottom of unreduced block
      q <- defHiCM (n - 1)
      if q <= 0
        then return ()  -- fully deflated
        else do
          -- Find p — top of unreduced block
          p <- findLoCM (q - 1)
          -- Apply one CM QR step to [p..q]
          bidiagQRStepCM mbaD offD mbaE offE mbaU offU mm mbaV offV nn p q
          let !newStall = if q == lastQ then stall + 1 else 0
          goCM (iter + 1) q newStall mbaD offD mbaE offE mbaU offU mm mbaV offV nn n maxIter
  where
    -- AED: scan bottom w entries, deflating negligible superdiagonals
    aedScan _ 0 = return ()
    aedScan k remaining
      | k <= 0 = return ()
      | otherwise = do
          ek <- readRawD mbaE offE (k - 1)
          dk <- readRawD mbaD offD k
          dk1 <- readRawD mbaD offD (k - 1)
          let tol = 1e-14 * (abs dk1 + abs dk)
          if abs ek <= tol
            then do
              writeRawD mbaE offE (k - 1) 0
              aedScan (k - 1) (remaining - 1)
            else return ()  -- stop at first non-negligible entry

    defHiCM !hi
      | hi <= 0 = return 0
      | otherwise = do
          ehi <- readRawD mbaE offE (hi - 1)
          dhi <- readRawD mbaD offD hi
          dhi1 <- readRawD mbaD offD (hi - 1)
          let tol = 1e-14 * (abs dhi1 + abs dhi)
          if abs ehi <= tol
            then do
              writeRawD mbaE offE (hi - 1) 0
              defHiCM (hi - 1)
            else return hi

    findLoCM !idx
      | idx <= 0 = return 0
      | otherwise = do
          eidx <- readRawD mbaE offE (idx - 1)
          didx <- readRawD mbaD offD idx
          didx1 <- readRawD mbaD offD (idx - 1)
          let tol = 1e-14 * (abs didx1 + abs didx)
          if abs eidx <= tol
            then do
              writeRawD mbaE offE (idx - 1) 0
              return idx
            else findLoCM (idx - 1)
{-# NOINLINE goCM #-}

-- | One implicit-shift QR step on bidiagonal [lo..hi] using column-major U,V.
-- Same as bidiagQRStep but calls rawMutApplyGivensColumnsCM for SIMD.
bidiagQRStepCM :: MutableByteArray s -> Int  -- d + offset
               -> MutableByteArray s -> Int  -- e + offset
               -> MutableByteArray s -> Int -> Int  -- U_CM + offset + mm
               -> MutableByteArray s -> Int -> Int  -- V_CM + offset + nn
               -> Int -> Int  -- lo, hi
               -> ST s ()
bidiagQRStepCM mbaD offD mbaE offE mbaU offU mm mbaV offV nn lo hi = do
  -- Compute Wilkinson shift from trailing 2×2 of T = B^T B
  dhi1 <- readRawD mbaD offD (hi - 1)
  dhi  <- readRawD mbaD offD hi
  ehi1 <- readRawD mbaE offE (hi - 1)
  ehi2 <- if hi >= 2 then readRawD mbaE offE (hi - 2) else return 0

  let t11 = dhi1 * dhi1 + (if hi - 1 > lo then ehi2 * ehi2 else 0)
      t12 = dhi1 * ehi1
      t22 = dhi * dhi + ehi1 * ehi1
      delta = (t11 - t22) / 2
      signD = if delta >= 0 then 1 else -1
      mu = t22 - t12 * t12 / (delta + signD * sqrt (delta * delta + t12 * t12))

  dlo <- readRawD mbaD offD lo
  elo <- readRawD mbaE offE lo
  let y = dlo * dlo - mu
      z = dlo * elo

  goChase lo y z
  where
    goChase k y_ z_ = do
      let (cosR, sinR) = givens y_ z_
      dk  <- readRawD mbaD offD k
      ek  <- readRawD mbaE offE k
      dk1 <- readRawD mbaD offD (k + 1)

      let dk'  = cosR * dk + sinR * ek
          ek'  = -sinR * dk + cosR * ek
          bulgeL = sinR * dk1
          dk1'   = cosR * dk1

      writeRawD mbaD offD k dk'
      writeRawD mbaE offE k ek'
      writeRawD mbaD offD (k + 1) dk1'

      -- Update e[k-1]: right Givens rotates entry from row above
      when (k > lo) $
        writeRawD mbaE offE (k - 1) (cosR * y_ + sinR * z_)

      -- Accumulate right rotation into V (column-major, SIMD)
      rawMutApplyGivensColumnsCM mbaV offV nn cosR sinR k (k+1) nn

      let (cosL, sinL) = givens dk' bulgeL

      let dk''  = cosL * dk' + sinL * bulgeL
          ek''  = cosL * ek' + sinL * dk1'
          dk1'' = -sinL * ek' + cosL * dk1'

      writeRawD mbaD offD k dk''
      writeRawD mbaE offE k ek''
      writeRawD mbaD offD (k + 1) dk1''

      when (k + 1 < hi) $ do
        ek1 <- readRawD mbaE offE (k + 1)
        let bulgeR = sinL * ek1
            ek1'   = cosL * ek1
        writeRawD mbaE offE (k + 1) ek1'

        -- Accumulate left rotation into U (column-major, SIMD)
        rawMutApplyGivensColumnsCM mbaU offU mm cosL sinL k (k+1) mm

        goChase (k + 1) ek'' bulgeR

      when (k + 1 >= hi) $
        rawMutApplyGivensColumnsCM mbaU offU mm cosL sinL k (k+1) mm

-- | One implicit-shift QR step on bidiagonal [lo..hi].
-- Computes Wilkinson shift from bottom 2×2 of B^T B,
-- then chases bulge via Givens rotations.
bidiagQRStep :: MutableByteArray s -> Int  -- d + offset
             -> MutableByteArray s -> Int  -- e + offset
             -> MutableByteArray s -> Int -> Int  -- U + offset + ucols
             -> MutableByteArray s -> Int -> Int  -- V + offset + vcols
             -> Int -> Int  -- lo, hi
             -> ST s ()
bidiagQRStep mbaD offD mbaE offE mbaU offU ucols mbaV offV vcols lo hi = do
  -- Compute Wilkinson shift from trailing 2×2 of T = B^T B
  dhi1 <- readRawD mbaD offD (hi - 1)
  dhi  <- readRawD mbaD offD hi
  ehi1 <- readRawD mbaE offE (hi - 1)
  ehi2 <- if hi >= 2 then readRawD mbaE offE (hi - 2) else return 0

  -- T = B^T B trailing 2×2:
  -- t11 = d[hi-1]^2 + e[hi-2]^2  (e[hi-2] = 0 if hi-1 == lo)
  -- t12 = d[hi-1] * e[hi-1]
  -- t22 = d[hi]^2 + e[hi-1]^2
  let t11 = dhi1 * dhi1 + (if hi - 1 > lo then ehi2 * ehi2 else 0)
      t12 = dhi1 * ehi1
      t22 = dhi * dhi + ehi1 * ehi1

  -- Wilkinson shift: eigenvalue of [[t11,t12],[t12,t22]] closer to t22
  let delta = (t11 - t22) / 2
      signD = if delta >= 0 then 1 else -1
      mu = t22 - t12 * t12 / (delta + signD * sqrt (delta * delta + t12 * t12))

  -- Initial values for bulge chase
  dlo <- readRawD mbaD offD lo
  elo <- readRawD mbaE offE lo
  let y = dlo * dlo - mu
      z = dlo * elo

  -- Chase bulge from lo to hi
  go lo y z
  where
    go k y_ z_ = do
      -- Right Givens rotation G(k,k+1,θ) to zero z in [y; z]
      let (cosR, sinR) = givens y_ z_
      -- Apply to columns k, k+1 of B (affects d[k], e[k], d[k+1], and possibly e[k-1])
      dk  <- readRawD mbaD offD k
      ek  <- readRawD mbaE offE k
      dk1 <- readRawD mbaD offD (k + 1)

      -- B * G^T: columns k and k+1 get mixed
      let dk'  = cosR * dk + sinR * ek
          ek'  = -sinR * dk + cosR * ek
          -- This creates a bulge at B[k+1,k]
          bulgeL = sinR * dk1
          dk1'   = cosR * dk1

      writeRawD mbaD offD k dk'
      writeRawD mbaE offE k ek'
      writeRawD mbaD offD (k + 1) dk1'

      -- Update e[k-1]: the right Givens also rotates the entry from the row above.
      -- For k > lo: B[k-1,k] was y_, B[k-1,k+1] was z_ (the bulge).
      -- After rotation: B[k-1,k] = cosR*y_ + sinR*z_ (= r), B[k-1,k+1] = 0.
      when (k > lo) $
        writeRawD mbaE offE (k - 1) (cosR * y_ + sinR * z_)

      -- Accumulate right rotation into V (columns k, k+1)
      rawMutApplyGivensColumns mbaV offV vcols cosR sinR k (k+1) vcols

      -- Left Givens rotation G(k,k+1,θ) to zero the bulge at (k+1, k)
      let (cosL, sinL) = givens dk' bulgeL

      -- G * B: rows k and k+1 get mixed
      let dk''  = cosL * dk' + sinL * bulgeL
          ek''  = cosL * ek' + sinL * dk1'
          dk1'' = -sinL * ek' + cosL * dk1'

      writeRawD mbaD offD k dk''
      writeRawD mbaE offE k ek''
      writeRawD mbaD offD (k + 1) dk1''

      -- This may create a new bulge at position (k, k+2) if k+1 < hi
      when (k + 1 < hi) $ do
        ek1 <- readRawD mbaE offE (k + 1)
        let bulgeR = sinL * ek1
            ek1'   = cosL * ek1
        writeRawD mbaE offE (k + 1) ek1'

        -- Accumulate left rotation into U (columns k, k+1)
        rawMutApplyGivensColumns mbaU offU ucols cosL sinL k (k+1) ucols

        -- Continue chase
        go (k + 1) ek'' bulgeR

      when (k + 1 >= hi) $
        -- Accumulate final left rotation
        rawMutApplyGivensColumns mbaU offU ucols cosL sinL k (k+1) ucols

-- | Compute Givens rotation coefficients (c, s) such that
-- @[c, s; -s, c] [a; b] = [r; 0]@, i.e. @-s*a + c*b = 0@ and @r = c*a + s*b > 0@.
--
-- This convention is chosen so that the bidiag QR bulge-chase formulas
-- @dk' = c*dk + s*ek@, @ek' = -s*dk + c*ek@ etc. are directly correct
-- for both left and right Givens rotations (GVL4 Algorithm 8.6.2).
givens :: Double -> Double -> (Double, Double)
givens a b
  | b == 0    = (1, 0)
  | abs b > abs a =
      let tau = a / b
          s   = 1 / sqrt (1 + tau * tau)
          c   = s * tau
      in (c, s)
  | otherwise =
      let tau = b / a
          c   = 1 / sqrt (1 + tau * tau)
          s   = c * tau
      in (c, s)
{-# INLINE givens #-}

-- | Accumulate a right Householder reflector into V.
-- Right reflector k: v stored in row k of frozen A, cols k+2..n-1, with v[k+1]=1 (implicit).
-- V = V * (I - beta * v * v^T)
-- For each row of V: V[row, k+1..n-1] -= (beta * Σ V[row,l] * v[l]) * v
rightQAccum :: MutableByteArray s -> Int -> Int  -- V + offset + vcols
            -> ByteArray -> Int -> Int           -- frozen A + offset + acols
            -> Double -> Int -> Int -> Int       -- beta, k, n, row
            -> ST s ()
rightQAccum mbaV offV vcols (ByteArray baA) offFA acols beta k nn row = do
  -- Phase 1: wi = beta * (V[row, k+1] + Σ_{l=k+2}^{n-1} V[row, l] * A[k, l])
  qrk1 <- readRawD mbaV offV (row * vcols + (k+1))
  acc <- goSum (k+2) 0
  let wi = beta * (qrk1 + acc)
  -- Phase 2: V[row, k+1] -= wi (implicit v[k+1]=1)
  writeRawD mbaV offV (row * vcols + (k+1)) (qrk1 - wi)
  -- V[row, l] -= wi * A[k, l] for l in k+2..n-1
  goUpdate (k+2) wi
  where
    goSum l acc_
      | l >= nn = return acc_
      | otherwise = do
          let vl = readBAI baA offFA (k * acols + l)
          qrl <- readRawD mbaV offV (row * vcols + l)
          goSum (l + 1) (acc_ + qrl * vl)

    goUpdate l wi
      | l >= nn = return ()
      | otherwise = do
          let vl = readBAI baA offFA (k * acols + l)
          qrl <- readRawD mbaV offV (row * vcols + l)
          writeRawD mbaV offV (row * vcols + l) (qrl - wi * vl)
          goUpdate (l + 1) wi

    readBAI ba_ off_ i_ =
      case indexDoubleArray# ba_ (case off_ of I# o -> o +# case i_ of I# ii -> ii) of
        v -> D# v
{-# INLINE rightQAccum #-}

-- | Read a Double from a MutableByteArray at element index.
readRawD :: MutableByteArray s -> Int -> Int -> ST s Double
readRawD (MutableByteArray mba) (I# off) (I# i) = ST $ \s ->
  case readDoubleArray# mba (off +# i) s of (# s', v #) -> (# s', D# v #)
{-# INLINE readRawD #-}

-- | Write a Double to a MutableByteArray at element index.
writeRawD :: MutableByteArray s -> Int -> Int -> Double -> ST s ()
writeRawD (MutableByteArray mba) (I# off) (I# i) (D# v) = ST $ \s ->
  case writeDoubleArray# mba (off +# i) v s of s' -> (# s', () #)
{-# INLINE writeRawD #-}

-- ============================================================================
-- Divide-and-conquer bidiagonal SVD (Gu-Eisenstat 1995, cf. LAPACK DBDSDC)
-- ============================================================================

-- | Small-subproblem threshold for D&C bidiagonal SVD.
-- Subproblems at or below this size use bidiag QR iteration.
dcBidiagThreshold :: Int
dcBidiagThreshold = 25

-- | Direct 2×2 upper bidiagonal SVD.
-- Given [[d0, e0], [0, d1]], compute singular values and rotation angles.
-- Returns (sigma_large, sigma_small, c_left, s_left, c_right, s_right) where
-- the left and right Givens rotations diagonalise B.
bidiag2x2SVD :: Double -> Double -> Double -> (Double, Double, Double, Double, Double, Double)
bidiag2x2SVD d0 e0 d1 =
  -- B^T B = [[d0², d0*e0], [d0*e0, e0²+d1²]]
  -- Eigenvalues of this 2×2 symmetric matrix give σ²
  let !a11 = d0 * d0
      !a12 = d0 * e0
      !a22 = e0 * e0 + d1 * d1
      !tr  = a11 + a22
      !det = a11 * a22 - a12 * a12
      !disc = max 0 (tr * tr - 4 * det)
      !sqrtDisc = sqrt disc
      !lam1 = (tr + sqrtDisc) / 2
      !lam2 = (tr - sqrtDisc) / 2
      !sig1 = sqrt (max 0 lam1)
      !sig2 = sqrt (max 0 lam2)
      -- Right rotation angle: diagonalise B^T B
      -- (a11 - lam2) * cv + a12 * sv = 0 => cv/sv = -a12/(a11-lam2)
      -- Or equivalently: atan2(a12, lam1 - a22)
      !cv = if abs a12 < 1e-300
            then 1
            else let d = a11 - lam2
                     r = sqrt (d * d + a12 * a12)
                 in d / r
      !sv = if abs a12 < 1e-300
            then 0
            else let d = a11 - lam2
                     r = sqrt (d * d + a12 * a12)
                 in a12 / r
      -- Left rotation: B * V = U * Sigma
      -- (d0*cv + e0*sv, -d0*sv + e0*cv)   = (sig1*cu, sig2*su_neg)
      -- (d1*sv,         d1*cv)             = (sig1*(-su), sig2*cu)
      -- From first row: sig1*cu = d0*cv + e0*sv
      !bv00 = d0 * cv + e0 * sv
      !bv10 = d1 * sv
      !r_l = sqrt (bv00 * bv00 + bv10 * bv10)
      !cu = if r_l < 1e-300 then 1 else bv00 / r_l
      !su = if r_l < 1e-300 then 0 else bv10 / r_l
  in (sig1, sig2, cu, su, cv, sv)
{-# INLINE bidiag2x2SVD #-}

-- | Divide-and-conquer bidiagonal SVD.
-- Replaces bidiagQRIterPCM for computing the SVD of a bidiagonal matrix.
--
-- Input: d[0..nn-1] (diagonal), e[0..nn-2] (superdiagonal).
-- Output: d[] overwritten with singular values,
--         U (mm×mm) and V (nn×nn) updated with accumulated rotations.
dcBidiagSVD :: forall s. MutableByteArray s -> Int    -- d + offset
            -> MutableByteArray s -> Int    -- e + offset
            -> MutableByteArray s -> Int -> Int  -- U + offset + mm
            -> MutableByteArray s -> Int -> Int  -- V + offset + nn
            -> Int                          -- nn (bidiag dimension)
            -> Double                       -- tolerance
            -> ST s ()
dcBidiagSVD mbaD offD mbaE offE mbaU offU mm mbaV offV nn0 nn tol = do
  -- Pre-allocate all workspace once at maximum size
  let !maxN = nn
  wsLam    <- newByteArray (maxN * 8)       -- new eigenvalues (squared)
  wsZ      <- newByteArray (maxN * 8)       -- z-vector
  wsDSort  <- newByteArray (maxN * 8)       -- sorted d² values
  wsZSort  <- newByteArray (maxN * 8)       -- sorted z values
  wsDOrig  <- newByteArray (maxN * 8)       -- original d values (unsquared)
  wsIdx    <- newByteArray (maxN * 8)       -- sort permutation (stored as Double)
  wsPerm   <- newByteArray (maxN * 8)       -- deflation permutation (Int)
  wsW      <- newByteArray (maxN * maxN * 8)  -- V-eigenvector matrix W_V
  wsWU     <- newByteArray (maxN * maxN * 8)  -- U-eigenvector matrix W_U

  -- Local accumulators for V (nn×nn) and U (nn×nn)
  -- U-local is nn×nn because we track rotations in singular-value index space
  wsVlocal <- newByteArray (maxN * maxN * 8)
  wsUlocal <- newByteArray (maxN * maxN * 8)

  -- GEMM workspace
  wsVsub   <- newByteArray (maxN * maxN * 8)  -- V column extraction buffer
  wsVres   <- newByteArray (maxN * maxN * 8)  -- V GEMM result
  wsUsub   <- newByteArray (maxN * maxN * 8)  -- U column extraction buffer
  wsUres   <- newByteArray (maxN * maxN * 8)  -- U GEMM result
  wsQtemp  <- newByteArray (maxN * maxN * 8)  -- QR base case scratch

  -- Initialise local accumulators as identity
  rawZeroDoubles wsVlocal 0 (maxN * maxN)
  rawZeroDoubles wsUlocal 0 (maxN * maxN)
  forM_ [0..maxN-1] $ \i -> do
    writeRawD wsVlocal 0 (i * maxN + i) 1
    writeRawD wsUlocal 0 (i * maxN + i) 1

  let -- Convert global index to local
      toLocal g = g

      -- Apply a k×k rotation matrix to wsVlocal columns [colOff..colOff+k-1]
      applyRotToVlocal !colOff !k rotMat = do
        forM_ [0..k-1] $ \j ->
          rawCopyColumn wsVlocal 0 maxN (colOff + j) wsVsub 0 k j maxN
        baVsub <- unsafeFreezeByteArray wsVsub
        baRot  <- unsafeFreezeByteArray rotMat
        rawZeroDoubles wsVres 0 (maxN * k)
        rawGemmKernel baVsub 0 baRot 0 wsVres 0 maxN k k
        forM_ [0..k-1] $ \j ->
          rawCopyColumn wsVres 0 k j wsVlocal 0 maxN (colOff + j) maxN

      -- Apply a k×k rotation matrix to wsUlocal columns [colOff..colOff+k-1]
      applyRotToUlocal !colOff !k rotMat = do
        forM_ [0..k-1] $ \j ->
          rawCopyColumn wsUlocal 0 maxN (colOff + j) wsUsub 0 k j maxN
        baUsub <- unsafeFreezeByteArray wsUsub
        baRot  <- unsafeFreezeByteArray rotMat
        rawZeroDoubles wsUres 0 (maxN * k)
        rawGemmKernel baUsub 0 baRot 0 wsUres 0 maxN k k
        forM_ [0..k-1] $ \j ->
          rawCopyColumn wsUres 0 k j wsUlocal 0 maxN (colOff + j) maxN

      -- The recursive D&C function
      dcGo :: Int -> Int -> ST s ()  -- s from ScopedTypeVariables
      dcGo lo hi
        -- Trivial: single element
        | lo >= hi = return ()

        -- Base case: 2×2 direct SVD
        | hi == lo + 1 = do
            d0_ <- readRawD mbaD offD lo
            e0_ <- readRawD mbaE offE lo
            d1_ <- readRawD mbaD offD hi
            let (!sig1, !sig2, !cu, !su, !cv, !sv) = bidiag2x2SVD d0_ e0_ d1_
            writeRawD mbaD offD lo sig1
            writeRawD mbaD offD hi sig2
            writeRawD mbaE offE lo 0
            -- Apply left Givens to Ulocal columns
            let !loL = toLocal lo
                !hiL = toLocal hi
            rawMutApplyGivensColumns wsUlocal 0 maxN cu su loL hiL maxN
            -- Apply right Givens to Vlocal columns
            rawMutApplyGivensColumns wsVlocal 0 maxN cv sv loL hiL maxN

        -- Small subproblem: use bidiag QR + GEMM to local accumulators
        | hi - lo + 1 <= dcBidiagThreshold = do
            let !k = hi - lo + 1
                !loL = toLocal lo
            -- Initialise k×k identities for U and V rotations
            rawZeroDoubles wsQtemp 0 (k * k)
            rawZeroDoubles wsW 0 (k * k)
            forM_ [0..k-1] $ \i -> do
              writeRawD wsQtemp 0 (i * k + i) 1
              writeRawD wsW     0 (i * k + i) 1
            -- Run bidiag QR iteration: wsQtemp accumulates left, wsW accumulates right
            bidiagQRIterP mbaD (offD + lo) mbaE (offE + lo) wsQtemp 0 k wsW 0 k k (30 * k)
            -- Apply rotations to local accumulators via GEMM
            applyRotToUlocal loL k wsQtemp
            applyRotToVlocal loL k wsW

        -- D&C merge
        | otherwise = do
            let !k   = (lo + hi) `div` 2
                !n1  = k - lo + 1
                !n2  = hi - k
                !nn_ = hi - lo + 1
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

            -- Recurse on left [lo..k] and right [k+1..hi] subproblems
            dcGo lo k
            dcGo (k + 1) hi

            -- === Merge phase ===

            -- Extract z-vector from Vlocal accumulator rows
            -- z[0..n1-1] = last row of V₁ = row kL, columns loL..loL+n1-1
            forM_ [0..n1-1] $ \i -> do
              qv <- readRawD wsVlocal 0 (kL * maxN + (loL + i))
              writeRawD wsZ 0 i qv
            -- z[n1..nn_-1] = first row of V₂ = row (kL+1), columns loL+n1..loL+nn_-1
            forM_ [0..n2-1] $ \i -> do
              qv <- readRawD wsVlocal 0 ((kL + 1) * maxN + (loL + n1 + i))
              let !zv = if beta < 0 then negate qv else qv
              writeRawD wsZ 0 (n1 + i) zv

            -- Save original d-values (unsquared) for U-eigenvector computation
            forM_ [0..nn_-1] $ \i -> do
              di <- readRawD mbaD offD (lo + i)
              writeRawD wsDOrig 0 i di

            -- Square d-values for secular equation and copy into sort buffers
            forM_ [0..nn_-1] $ \i -> do
              di <- readRawD mbaD offD (lo + i)
              writeRawD wsDSort 0 i (di * di)
              writeRawD wsZSort 0 i =<< readRawD wsZ 0 i
              writeRawD wsIdx 0 i (fromIntegral i)

            -- Sort by d² values (insertion sort)
            forM_ [1..nn_-1] $ \i -> do
              di   <- readRawD wsDSort 0 i
              zi   <- readRawD wsZSort 0 i
              idxi <- readRawD wsIdx 0 i
              dOi  <- readRawD wsDOrig 0 i
              let insertAt !j
                    | j < 0 = do
                        writeRawD wsDSort 0 0 di
                        writeRawD wsZSort 0 0 zi
                        writeRawD wsIdx   0 0 idxi
                        writeRawD wsDOrig 0 0 dOi
                    | otherwise = do
                        dj <- readRawD wsDSort 0 j
                        if dj > di
                          then do
                            writeRawD wsDSort 0 (j+1) dj
                            writeRawD wsZSort 0 (j+1) =<< readRawD wsZSort 0 j
                            writeRawD wsIdx   0 (j+1) =<< readRawD wsIdx 0 j
                            writeRawD wsDOrig 0 (j+1) =<< readRawD wsDOrig 0 j
                            insertAt (j - 1)
                          else do
                            writeRawD wsDSort 0 (j+1) di
                            writeRawD wsZSort 0 (j+1) zi
                            writeRawD wsIdx   0 (j+1) idxi
                            writeRawD wsDOrig 0 (j+1) dOi
              insertAt (i - 1)

            -- Close-d deflation on squared values (cf. LAPACK dlaed2)
            dMaxSq <- readRawD wsDSort 0 (nn_ - 1)
            dMinSq <- readRawD wsDSort 0 0
            let !closeDTol = 8 * 2.220446049250313e-16
                          * max (abs dMaxSq) (abs dMinSq + rho)
            forM_ [0..nn_-2] $ \i -> do
              di  <- readRawD wsDSort 0 i
              di1 <- readRawD wsDSort 0 (i + 1)
              when (abs (di1 - di) <= closeDTol) $ do
                zi  <- readRawD wsZSort 0 i
                zi1 <- readRawD wsZSort 0 (i + 1)
                let !r = sqrt (zi * zi + zi1 * zi1)
                when (r > 1e-300) $ do
                  let !c = zi1 / r
                      !s = zi / r
                  writeRawD wsZSort 0 i 0
                  writeRawD wsZSort 0 (i + 1) r
                  -- Apply Givens to Vlocal columns (same as tridiagonal)
                  origI  <- readRawD wsIdx 0 i
                  origI1 <- readRawD wsIdx 0 (i + 1)
                  let !colI  = loL + (round origI  :: Int)
                      !colI1 = loL + (round origI1 :: Int)
                  rawMutApplyGivensColumns wsVlocal 0 maxN c s colI colI1 maxN
                  -- Also apply to Ulocal
                  rawMutApplyGivensColumns wsUlocal 0 maxN c s colI colI1 maxN

            -- Perturbation-based deflation
            zn2 <- sumZSq wsZSort 0 nn_
            let !eps_ = 2.220446049250313e-16
                !matNorm = max (abs dMaxSq) (abs dMinSq) + rho * zn2
                !basicDeflTol = max (tol * sqrt zn2) (8 * eps_ * matNorm)
                !pertDeflTol = sqrt (eps_ * (1 + matNorm) / max rho 1e-300)
                !deflTol = max basicDeflTol pertDeflTol
            kND <- deflatePartition wsZSort 0 wsPerm 0 nn_ deflTol

            -- Extract Vlocal columns permuted by sort order into wsVsub (maxN × nn_)
            forM_ [0..nn_-1] $ \sortedJ -> do
              origIdx <- readRawD wsIdx 0 sortedJ
              let !origJ = round origIdx :: Int
                  !srcCol = loL + origJ
              rawCopyColumn wsVlocal 0 maxN srcCol wsVsub 0 nn_ sortedJ maxN

            -- Also extract Ulocal columns
            forM_ [0..nn_-1] $ \sortedJ -> do
              origIdx <- readRawD wsIdx 0 sortedJ
              let !origJ = round origIdx :: Int
                  !srcCol = loL + origJ
              rawCopyColumn wsUlocal 0 maxN srcCol wsUsub 0 nn_ sortedJ maxN

            if kND == 0
              then do
                -- All deflated: singular values from sorted d², vectors from sorted cols
                forM_ [0..nn_-1] $ \i -> do
                  rawCopyColumn wsVsub 0 nn_ i wsVlocal 0 maxN (loL + i) maxN
                  rawCopyColumn wsUsub 0 nn_ i wsUlocal 0 maxN (loL + i) maxN
                forM_ [0..nn_-1] $ \i -> do
                  dsq <- readRawD wsDSort 0 i
                  writeRawD mbaD offD (lo + i) (sqrt (max 0 dsq))

              else if kND == nn_
                then do
                  -- No deflation: full secular solve + eigenvectors + dual GEMM
                  secularSolve wsLam 0 wsDSort 0 wsZSort 0 rho nn_ deflTol

                  -- V-eigenvectors via Gu-Eisenstat on (d², z, μ)
                  dcEigenvectors wsW 0 wsDSort 0 wsZSort 0 wsLam 0 rho nn_

                  -- U-eigenvectors: W_U[j,i] = dOrig[j] * W_V[j,i], then normalize
                  dcEigenvectorsBidiagU wsWU 0 wsDOrig 0 wsW 0 nn_

                  -- V-GEMM: wsVres = Vsub * W_V
                  do baVs <- unsafeFreezeByteArray wsVsub
                     baWv <- unsafeFreezeByteArray wsW
                     rawZeroDoubles wsVres 0 (maxN * nn_)
                     rawGemmKernel baVs 0 baWv 0 wsVres 0 maxN nn_ nn_
                  forM_ [0..nn_-1] $ \i ->
                    rawCopyColumn wsVres 0 nn_ i wsVlocal 0 maxN (loL + i) maxN

                  -- U-GEMM: wsUres = Usub * W_U
                  do baUs <- unsafeFreezeByteArray wsUsub
                     baWu <- unsafeFreezeByteArray wsWU
                     rawZeroDoubles wsUres 0 (maxN * nn_)
                     rawGemmKernel baUs 0 baWu 0 wsUres 0 maxN nn_ nn_
                  forM_ [0..nn_-1] $ \i ->
                    rawCopyColumn wsUres 0 nn_ i wsUlocal 0 maxN (loL + i) maxN

                  -- Write singular values = sqrt(|mu|)
                  forM_ [0..nn_-1] $ \i -> do
                    mu <- readRawD wsLam 0 i
                    writeRawD mbaD offD (lo + i) (sqrt (max 0 (abs mu)))

                else do
                  -- Partial deflation: reduced secular solve + reduced GEMM
                  -- Build compressed d²_nd and z_nd in wsQtemp
                  forM_ [0..kND-1] $ \j -> do
                    pi_ <- readRawI wsPerm 0 j
                    dpi <- readRawD wsDSort 0 pi_
                    zpi <- readRawD wsZSort 0 pi_
                    writeRawD wsQtemp 0 j dpi
                    writeRawD wsQtemp 0 (kND + j) zpi

                  -- Also build compressed dOrig_nd for U-eigenvectors
                  forM_ [0..kND-1] $ \j -> do
                    pi_ <- readRawI wsPerm 0 j
                    dO  <- readRawD wsDOrig 0 pi_
                    writeRawD wsQtemp 0 (2 * kND + j) dO

                  secularSolve wsLam 0 wsQtemp 0 wsQtemp kND rho kND deflTol
                  dcEigenvectors wsW 0 wsQtemp 0 wsQtemp kND wsLam 0 rho kND
                  dcEigenvectorsBidiagU wsWU 0 wsQtemp (2 * kND) wsW 0 kND

                  -- Copy deflated columns to local accumulators
                  forM_ [kND..nn_-1] $ \j -> do
                    pi_ <- readRawI wsPerm 0 j
                    rawCopyColumn wsVsub 0 nn_ pi_ wsVlocal 0 maxN (loL + j) maxN
                    rawCopyColumn wsUsub 0 nn_ pi_ wsUlocal 0 maxN (loL + j) maxN

                  -- Extract V_nd (maxN × kND) from non-deflated columns
                  forM_ [0..kND-1] $ \j -> do
                    pi_ <- readRawI wsPerm 0 j
                    rawCopyColumn wsVsub 0 nn_ pi_ wsVres 0 kND j maxN

                  -- V-GEMM: wsVsub(maxN×kND) = V_nd(maxN×kND) * W_V(kND×kND)
                  do baVnd <- unsafeFreezeByteArray wsVres
                     baWv  <- unsafeFreezeByteArray wsW
                     rawZeroDoubles wsVsub 0 (maxN * kND)
                     rawGemmKernel baVnd 0 baWv 0 wsVsub 0 maxN kND kND
                  forM_ [0..kND-1] $ \j ->
                    rawCopyColumn wsVsub 0 kND j wsVlocal 0 maxN (loL + j) maxN

                  -- Extract U_nd (maxN × kND) from non-deflated columns
                  forM_ [0..kND-1] $ \j -> do
                    pi_ <- readRawI wsPerm 0 j
                    rawCopyColumn wsUsub 0 nn_ pi_ wsUres 0 kND j maxN

                  -- U-GEMM: wsUsub(maxN×kND) = U_nd(maxN×kND) * W_U(kND×kND)
                  do baUnd <- unsafeFreezeByteArray wsUres
                     baWu  <- unsafeFreezeByteArray wsWU
                     rawZeroDoubles wsUsub 0 (maxN * kND)
                     rawGemmKernel baUnd 0 baWu 0 wsUsub 0 maxN kND kND
                  forM_ [0..kND-1] $ \j ->
                    rawCopyColumn wsUsub 0 kND j wsUlocal 0 maxN (loL + j) maxN

                  -- Write eigenvalues: non-deflated from secular, deflated from sorted d²
                  forM_ [0..kND-1] $ \i -> do
                    mu <- readRawD wsLam 0 i
                    writeRawD mbaD offD (lo + i) (sqrt (max 0 (abs mu)))
                  forM_ [kND..nn_-1] $ \j -> do
                    pi_ <- readRawI wsPerm 0 j
                    dsq <- readRawD wsDSort 0 pi_
                    writeRawD mbaD offD (lo + j) (sqrt (max 0 dsq))

  -- Run the D&C recursion
  dcGo 0 (nn - 1)

  -- Final step: apply local accumulators to global U and V via GEMM
  -- V[:, 0..nn-1] = V[:, 0..nn-1] * Vlocal
  forM_ [0..nn-1] $ \j ->
    rawCopyColumn mbaV offV nn0 j wsVsub 0 nn j nn0
  do baVs <- unsafeFreezeByteArray wsVsub
     baVl <- unsafeFreezeByteArray wsVlocal
     rawZeroDoubles wsVres 0 (nn0 * nn)
     rawGemmKernel baVs 0 baVl 0 wsVres 0 nn0 nn nn
  forM_ [0..nn-1] $ \j ->
    rawCopyColumn wsVres 0 nn j mbaV offV nn0 j nn0

  -- U[:, 0..nn-1] = U[:, 0..nn-1] * Ulocal
  forM_ [0..nn-1] $ \j ->
    rawCopyColumn mbaU offU mm j wsUsub 0 nn j mm
  do baUs <- unsafeFreezeByteArray wsUsub
     baUl <- unsafeFreezeByteArray wsUlocal
     rawZeroDoubles wsUres 0 (mm * nn)
     rawGemmKernel baUs 0 baUl 0 wsUres 0 mm nn nn
  forM_ [0..nn-1] $ \j ->
    rawCopyColumn wsUres 0 nn j mbaU offU mm j mm

-- | Compute U-eigenvectors from V-eigenvectors and original (unsquared) d-values.
-- W_U[j,i] = dOrig[j] * W_V[j,i], then normalize each column.
dcEigenvectorsBidiagU :: MutableByteArray s -> Int  -- W_U output (nn × nn)
                      -> MutableByteArray s -> Int  -- dOrig (unsquared singular values)
                      -> MutableByteArray s -> Int  -- W_V (V-eigenvectors, already computed)
                      -> Int                        -- nn
                      -> ST s ()
dcEigenvectorsBidiagU mbaWU offWU mbaDOrig offDO mbaWV offWV nn = do
  forM_ [0..nn-1] $ \i -> do
    -- W_U[:,i] = diag(dOrig) * W_V[:,i], then normalize
    norm2 <- goCol i 0 0
    let !invNorm = if norm2 > 0 then 1 / sqrt norm2 else 1
    forM_ [0..nn-1] $ \j -> do
      wuji <- readRawD mbaWU offWU (j * nn + i)
      writeRawD mbaWU offWU (j * nn + i) (wuji * invNorm)
  where
    goCol !i !j !acc
      | j >= nn = pure acc
      | otherwise = do
          dj   <- readRawD mbaDOrig offDO j
          wvji <- readRawD mbaWV offWV (j * nn + i)
          let !wuji = dj * wvji
          writeRawD mbaWU offWU (j * nn + i) wuji
          goCol i (j + 1) (acc + wuji * wuji)
{-# NOINLINE dcBidiagSVD #-}
