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
import Numeric.LinearAlgebra.Massiv.BLAS.Level3 (matMul, transpose, matMulAtAP)
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
      -- Compute U = A · V · diag(1/σ) via pre-scaled V and single GEMM.
      -- This avoids the intermediate m×n AV matrix and a separate scaling pass.
      -- V_scaled[i,j] = V[i,j] / σ_j (zero for σ_j ≤ ε).
      !vScaled = createMatrix @n @n @M.P @Double $ \mvs -> do
        let !baV   = unwrapByteArray (unMatrix v)
            !offVs = unwrapByteArrayOffset (unMatrix v)
            !mbaVS = unwrapMutableByteArray mvs
            !offVS = unwrapMutableByteArrayOffset mvs
        -- Pre-compute invSigma
        mbaInvS <- newByteArray (nn * 8)
        forM_ [0..nn-1] $ \j -> do
          let sj = sigma !. j
          writeRawD mbaInvS 0 j (if sj > 1e-14 then 1.0 / sj else 0.0)
        !(ByteArray baInvS#) <- unsafeFreezeByteArray mbaInvS
        -- Scale each column: V_scaled[i,j] = V[i,j] * invSigma[j]
        let !(ByteArray baV#) = baV
            !(I# offVs#) = offVs
            !(MutableByteArray mbaVS#) = mbaVS
            !(I# offVS#) = offVS
            !(I# nnV#) = nn
            !nn4 = nn - (nn `rem` 4)
            !(I# nn4#) = nn4
        ST $ \s0 ->
          let goRow i s
                | isTrue# (i >=# nnV#) = s
                | otherwise =
                    let !srcOff = offVs# +# i *# nnV#
                        !dstOff = offVS# +# i *# nnV#
                        goSimd j s1
                          | isTrue# (j >=# nn4#) = s1
                          | otherwise =
                              let vv = indexDoubleArrayAsDoubleX4# baV# (srcOff +# j)
                                  sv = indexDoubleArrayAsDoubleX4# baInvS# j
                                  !p  = timesDoubleX4# vv sv
                              in case writeDoubleArrayAsDoubleX4# mbaVS# (dstOff +# j) p s1 of
                                   s2 -> goSimd (j +# 4#) s2
                        goScalar j s1
                          | isTrue# (j >=# nnV#) = s1
                          | otherwise =
                              let vVal = indexDoubleArray# baV# (srcOff +# j)
                                  sVal = indexDoubleArray# baInvS# j
                              in case writeDoubleArray# mbaVS# (dstOff +# j) (vVal *## sVal) s1 of
                                   s2 -> goScalar (j +# 1#) s2
                    in goRow (i +# 1#) (goScalar nn4# (goSimd 0# s))
          in (# goRow 0# s0, () #)
      -- U = A · V_scaled: GEMM writes m×n result directly.
      -- For mm == nn (square), GEMM writes directly into U.
      -- For mm > nn (rectangular), GEMM writes into temp then copy columns.
      u = createMatrix @m @m @M.P $ \mu -> do
        let !mbaU  = unwrapMutableByteArray mu
            !offU  = unwrapMutableByteArrayOffset mu
            !(I# mm#) = mm
        -- Zero all of U
        rawZeroDoubles mbaU offU (mm * mm)
        let !baA  = unwrapByteArray (unMatrix a)
            !offA = unwrapByteArrayOffset (unMatrix a)
            !baVS = unwrapByteArray (unMatrix vScaled)
            !offVS = unwrapByteArrayOffset (unMatrix vScaled)
        if mm == nn
          then
            -- Direct GEMM into U (stride mm == nn, so layout matches)
            rawGemmKernel baA offA baVS offVS mbaU offU mm nn nn
          else do
            -- GEMM into temp (m×n), then copy columns into U (m×m)
            mbaTemp <- newByteArray (mm * nn * 8)
            rawZeroDoubles mbaTemp 0 (mm * nn)
            rawGemmKernel baA offA baVS offVS mbaTemp 0 mm nn nn
            baTemp <- unsafeFreezeByteArray mbaTemp
            -- Copy: U[i, 0..nn-1] = temp[i, 0..nn-1]
            let !(ByteArray baT#) = baTemp
                !(MutableByteArray mbaU#) = mbaU
                !(I# offU#) = offU
                !(I# nn#) = nn
            ST $ \s0 ->
              let goCopy i s
                    | isTrue# (i >=# mm#) = s
                    | otherwise =
                        let goCol j s1
                              | isTrue# (j >=# nn#) = s1
                              | otherwise =
                                  let !val = indexDoubleArray# baT# (i *# nn# +# j)
                                  in case writeDoubleArray# mbaU# (offU# +# i *# mm# +# j) val s1 of
                                       s2 -> goCol (j +# 1#) s2
                        in goCopy (i +# 1#) (goCol 0# s)
              in (# goCopy 0# s0, () #)
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

  -- Phase 1: Bidiagonalise A in-place (BLAS-3 panel for large, Level-2 for small)
  bidiagonalizePPanel mbaA offA mm nn mbaBetaL mbaBetaR

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

-- | BLAS-3 panel bidiagonalisation (DLABRD-style, GVL4 §5.4.3).
-- Processes nb columns at a time, deferring trailing updates via X, Y
-- accumulators and applying them as rank-nb GEMMs.
-- Falls back to Level-2 bidiagonalizeP for n < panelBidiagCrossover.
bidiagonalizePPanel :: MutableByteArray s -> Int -> Int -> Int
                    -> MutableByteArray s -> MutableByteArray s -> ST s ()
bidiagonalizePPanel mbaA offA mm nn mbaBetaL mbaBetaR
  | nn < panelBidiagCrossover = bidiagonalizeP mbaA offA mm nn mbaBetaL mbaBetaR
  | otherwise = do
      let !nb = min 32 (max 8 (nn `div` 6))
      -- Allocate accumulators: X (mm × nb), Y (nn × nb), row-major
      mbaX <- newByteArray (mm * nb * 8)
      mbaY <- newByteArray (nn * nb * 8)
      -- Temp vectors for dot products
      mbaZL <- newByteArray (nb * 8)  -- V_L^T * v or Y^T * u
      mbaZX <- newByteArray (nb * 8)  -- X^T * v or V_R^T * u
      -- Buffers for trailing GEMM
      mbaVLbuf <- newByteArray (mm * nb * 8)
      mbaYTbuf <- newByteArray (nb * nn * 8)
      mbaXbuf  <- newByteArray (mm * nb * 8)
      mbaVRTbuf <- newByteArray (nb * nn * 8)
      mbaTrail <- newByteArray (mm * nn * 8)

      let goPanel !k0
            | k0 >= nn - 1 = pure ()
            | otherwise = do
                let !bs = min nb (nn - 1 - k0)
                if bs < 2  -- last column: use Level-2
                  then bidiagLastCols mbaA offA mm nn mbaBetaL mbaBetaR k0
                  else do
                    rawZeroDoubles mbaX 0 (mm * bs)
                    rawZeroDoubles mbaY 0 (nn * bs)
                    -- Panel phase
                    panelBidiagStep mbaA offA mm nn mbaBetaL mbaBetaR
                                    mbaX mbaY mbaZL mbaZX k0 bs
                    -- Trailing update
                    let !remR = mm - k0 - bs
                        !remC = nn - k0 - bs
                    when (remR > 0 && remC > 0) $
                      applyTrailingUpdate mbaA offA mm nn mbaX mbaY
                                          mbaVLbuf mbaYTbuf mbaXbuf mbaVRTbuf mbaTrail
                                          k0 bs remR remC
                    goPanel (k0 + bs)
      goPanel 0
  where
    panelBidiagCrossover = 128
{-# NOINLINE bidiagonalizePPanel #-}

-- | Finish remaining columns with Level-2 bidiagonalisation.
bidiagLastCols :: MutableByteArray s -> Int -> Int -> Int
               -> MutableByteArray s -> MutableByteArray s -> Int -> ST s ()
bidiagLastCols mbaA offA mm nn mbaBetaL mbaBetaR k0 = do
  forM_ [k0..nn-1] $ \k -> do
    -- Left Householder
    if k < mm - 1
      then do
        sigma <- rawMutSumSqColumn mbaA offA nn (k+1) mm k
        x0 <- readRawD mbaA offA (k * nn + k)
        if sigma < 1e-300
          then writeRawD mbaBetaL 0 k 0
          else do
            let mu = sqrt (x0 * x0 + sigma)
                v0 = if x0 <= 0 then x0 - mu else -sigma / (x0 + mu)
                beta = 2 * v0 * v0 / (sigma + v0 * v0)
            forM_ [k+1..mm-1] $ \i -> do
              aik <- readRawD mbaA offA (i * nn + k)
              writeRawD mbaA offA (i * nn + k) (aik / v0)
            writeRawD mbaA offA (k * nn + k) mu
            writeRawD mbaBetaL 0 k beta
            forM_ [k+1..nn-1] $ \col ->
              rawMutHouseholderApply mbaA offA nn beta k mm col
      else writeRawD mbaBetaL 0 k 0
    -- Right Householder
    if k < nn - 2
      then do
        sigma <- rawMutSumSqRow mbaA offA nn k (k+2) nn
        x0 <- readRawD mbaA offA (k * nn + (k+1))
        if sigma < 1e-300
          then writeRawD mbaBetaR 0 k 0
          else do
            let mu = sqrt (x0 * x0 + sigma)
                v0 = if x0 <= 0 then x0 - mu else -sigma / (x0 + mu)
                beta = 2 * v0 * v0 / (sigma + v0 * v0)
            forM_ [k+2..nn-1] $ \j -> do
              akj <- readRawD mbaA offA (k * nn + j)
              writeRawD mbaA offA (k * nn + j) (akj / v0)
            writeRawD mbaA offA (k * nn + (k+1)) mu
            writeRawD mbaBetaR 0 k beta
            forM_ [k+1..mm-1] $ \row ->
              rawMutHouseholderApplyRow mbaA offA nn beta k (k+1) nn row
      else when (k < nn - 1) $ writeRawD mbaBetaR 0 k 0

-- | DLABRD-style panel step: compute bs left/right Householder reflectors
-- starting at column k0, maintaining X and Y accumulators.
-- After this, A_eff[i,c] = A[i,c] - V_L[i,:]*Y[c,:]^T - X[i,:]*V_R[c,:]^T
-- for all i >= k0+bs, c >= k0+bs.
panelBidiagStep :: MutableByteArray s -> Int -> Int -> Int
                -> MutableByteArray s -> MutableByteArray s
                -> MutableByteArray s -> MutableByteArray s
                -> MutableByteArray s -> MutableByteArray s
                -> Int -> Int -> ST s ()
panelBidiagStep mbaA offA mm nn mbaBetaL mbaBetaR mbaX mbaY mbaZL mbaZX k0 bs = do
  forM_ [0..bs-1] $ \j -> do
    let !k = k0 + j

    -- ================================================================
    -- PART A: Left Householder on column k
    -- ================================================================

    -- Step A1: Read corrected column k into A (in-place correction for rows k..m-1).
    -- A_corr[i,k] = A[i,k] - sum_{l<j} V_L[i,l]*Y[k,l] - sum_{l<j} X[i,l]*V_R[k,l]
    when (j > 0) $
      forM_ [k..mm-1] $ \i -> do
        aik <- readRawD mbaA offA (i * nn + k)
        -- V_L[i,l] * Y[k,l] sum
        cVLY <- panelDot_VLY mbaA offA nn mbaY bs k0 i k j
        -- X[i,l] * V_R[k,l] sum
        cXVR <- panelDot_XVR mbaA offA nn mbaX bs k0 i k j
        writeRawD mbaA offA (i * nn + k) (aik - cVLY - cXVR)

    -- Step A2: Left Householder from corrected column k, rows k..m-1
    if k < mm - 1
      then do
        sigma <- rawMutSumSqColumn mbaA offA nn (k+1) mm k
        x0 <- readRawD mbaA offA (k * nn + k)
        if sigma < 1e-300
          then do
            writeRawD mbaBetaL 0 k 0
            -- Zero Y column j
            forM_ [0..nn-1] $ \c -> writeRawD mbaY 0 (c * bs + j) 0
          else do
            let mu = sqrt (x0 * x0 + sigma)
                v0 = if x0 <= 0 then x0 - mu else -sigma / (x0 + mu)
                beta = 2 * v0 * v0 / (sigma + v0 * v0)
            -- Normalise and store HH vector in column k
            forM_ [k+1..mm-1] $ \i -> do
              aik <- readRawD mbaA offA (i * nn + k)
              writeRawD mbaA offA (i * nn + k) (aik / v0)
            writeRawD mbaA offA (k * nn + k) mu
            writeRawD mbaBetaL 0 k beta

            -- Step A3: Compute Y column j
            -- Precompute zL[l] = V_L[:,l]^T * v for l = 0..j-1
            forM_ [0..j-1] $ \l -> do
              d <- dotVL_v mbaA offA nn k0 k mm l
              writeRawD mbaZL 0 l d
            -- Precompute zX[l] = X[:,l]^T * v for l = 0..j-1
            forM_ [0..j-1] $ \l -> do
              d <- dotX_v mbaA offA nn mbaX bs k mm l
              writeRawD mbaZX 0 l d
            -- Y[c, j] = beta * (A^T*v[c] - sum_l Y[c,l]*zL[l] - sum_l V_R[c,l]*zX[l])
            forM_ [0..k] $ \c -> writeRawD mbaY 0 (c * bs + j) 0
            forM_ [k+1..nn-1] $ \c -> do
              atv <- dotAT_v mbaA offA nn k mm c
              ycorr <- dotAccum mbaY bs c mbaZL j
              vrcorr <- dotVR_zX mbaA offA nn k0 mbaZX c j
              writeRawD mbaY 0 (c * bs + j) (beta * (atv - ycorr - vrcorr))
      else do
        writeRawD mbaBetaL 0 k 0
        forM_ [0..nn-1] $ \c -> writeRawD mbaY 0 (c * bs + j) 0

    -- ================================================================
    -- PART B: Correct row k, then right Householder (if applicable)
    -- ================================================================

    -- Step B1: ALWAYS correct row k for columns k+1..n-1.
    -- This is needed both for the right HH (if k < nn-2) and for the
    -- superdiagonal entry e[k] = A[k, k+1] and trailing column values.
    -- A_eff[k,c] = A[k,c] - sum_{l<=j} V_L[k,l]*Y[c,l] - sum_{l<j} X[k,l]*V_R[c,l]
    when (k < nn - 1) $
      forM_ [k+1..nn-1] $ \c -> do
        akc <- readRawD mbaA offA (k * nn + c)
        cVLY <- panelDot_VLY mbaA offA nn mbaY bs k0 k c (j+1)
        cXVR <- panelDot_XVR mbaA offA nn mbaX bs k0 k c j
        writeRawD mbaA offA (k * nn + c) (akc - cVLY - cXVR)

    -- Step B2: Right Householder from corrected row k (only if k < nn-2)
    if k < nn - 2
      then do
        sigma <- rawMutSumSqRow mbaA offA nn k (k+2) nn
        x0 <- readRawD mbaA offA (k * nn + (k+1))
        if sigma < 1e-300
          then do
            writeRawD mbaBetaR 0 k 0
            forM_ [0..mm-1] $ \i -> writeRawD mbaX 0 (i * bs + j) 0
          else do
            let mu = sqrt (x0 * x0 + sigma)
                v0 = if x0 <= 0 then x0 - mu else -sigma / (x0 + mu)
                gamma = 2 * v0 * v0 / (sigma + v0 * v0)
            forM_ [k+2..nn-1] $ \c -> do
              akc <- readRawD mbaA offA (k * nn + c)
              writeRawD mbaA offA (k * nn + c) (akc / v0)
            writeRawD mbaA offA (k * nn + (k+1)) mu
            writeRawD mbaBetaR 0 k gamma

            -- Step B3: Compute X column j
            -- Precompute zL'[l] = Y[:,l]^T * u for l = 0..j
            forM_ [0..j] $ \l -> do
              d <- dotY_u mbaA offA nn mbaY bs k l
              writeRawD mbaZL 0 l d
            -- Precompute zX'[l] = V_R[:,l]^T * u for l = 0..j-1
            forM_ [0..j-1] $ \l -> do
              d <- dotVR_u mbaA offA nn k0 k l
              writeRawD mbaZX 0 l d
            -- X[i, j] = gamma * (A*u[i] - sum_l V_L[i,l]*zL'[l] - sum_l X[i,l]*zX'[l])
            forM_ [0..k] $ \i -> writeRawD mbaX 0 (i * bs + j) 0
            forM_ [k+1..mm-1] $ \i -> do
              au <- dotA_u mbaA offA nn k i
              vlcorr <- dotVL_zL mbaA offA nn k0 mbaZL i (j+1)
              xcorr <- dotX_zX mbaX bs mbaZX i j
              writeRawD mbaX 0 (i * bs + j) (gamma * (au - vlcorr - xcorr))
      else do
        when (k < nn - 1) $ writeRawD mbaBetaR 0 k 0
        forM_ [0..mm-1] $ \i -> writeRawD mbaX 0 (i * bs + j) 0

-- Helper: sum_l V_L[i,l]*Y[c,l] for l = 0..nL-1
panelDot_VLY :: MutableByteArray s -> Int -> Int -> MutableByteArray s -> Int
             -> Int -> Int -> Int -> Int -> ST s Double
panelDot_VLY mbaA offA nn mbaY bs k0 i c nL = go 0 0
  where
    go !l !acc
      | l >= nL = pure acc
      | otherwise = do
          let !kl = k0 + l
          vl <- if i == kl then pure 1.0
                else if i > kl then readRawD mbaA offA (i * nn + kl)
                else pure 0.0
          ycl <- readRawD mbaY 0 (c * bs + l)
          go (l+1) (acc + vl * ycl)

-- Helper: sum_l X[i,l]*V_R[c,l] for l = 0..nR-1
panelDot_XVR :: MutableByteArray s -> Int -> Int -> MutableByteArray s -> Int
             -> Int -> Int -> Int -> Int -> ST s Double
panelDot_XVR mbaA offA nn mbaX bs k0 i c nR = go 0 0
  where
    go !l !acc
      | l >= nR = pure acc
      | otherwise = do
          xil <- readRawD mbaX 0 (i * bs + l)
          let !kl = k0 + l
          vr <- if c == kl + 1 then pure 1.0
                else if c > kl + 1 then readRawD mbaA offA (kl * nn + c)
                else pure 0.0
          go (l+1) (acc + xil * vr)

-- Helper: V_L[:,l]^T * v where v = [1, A[k+1:m-1, k]]
dotVL_v :: MutableByteArray s -> Int -> Int -> Int -> Int -> Int -> Int -> ST s Double
dotVL_v mbaA offA nn k0 k mm l = do
  -- v[i-k]: v[0]=1, v[i-k]=A[i,k] for i>k
  -- V_L[i,l]: 1 if i==kl, A[i,kl] if i>kl, 0 if i<kl
  -- Since k >= k0+j and l < j, kl < k, so V_L[k,l] = A[k,kl]
  vlk <- readRawD mbaA offA (k * nn + kl)
  go (k+1) vlk
  where
    !kl = k0 + l
    go !i !acc
      | i >= mm = pure acc
      | otherwise = do
          vli <- readRawD mbaA offA (i * nn + kl)
          vi  <- readRawD mbaA offA (i * nn + k)
          go (i+1) (acc + vli * vi)

-- Helper: X[:,l]^T * v where v = [1, A[k+1:m-1, k]]
dotX_v :: MutableByteArray s -> Int -> Int -> MutableByteArray s -> Int
       -> Int -> Int -> Int -> ST s Double
dotX_v mbaA offA nn mbaX bs k mm l = do
  xkl <- readRawD mbaX 0 (k * bs + l)
  go (k+1) xkl
  where
    go !i !acc
      | i >= mm = pure acc
      | otherwise = do
          xil <- readRawD mbaX 0 (i * bs + l)
          vi  <- readRawD mbaA offA (i * nn + k)
          go (i+1) (acc + xil * vi)

-- Helper: A^T * v at column c, where v = [1, A[k+1:m-1, k]]
dotAT_v :: MutableByteArray s -> Int -> Int -> Int -> Int -> Int -> ST s Double
dotAT_v mbaA offA nn k mm c = do
  akc <- readRawD mbaA offA (k * nn + c)
  go (k+1) akc
  where
    go !i !acc
      | i >= mm = pure acc
      | otherwise = do
          aic <- readRawD mbaA offA (i * nn + c)
          vi  <- readRawD mbaA offA (i * nn + k)
          go (i+1) (acc + aic * vi)

-- Helper: sum_l Y[c,l]*zL[l] for l = 0..nL-1
dotAccum :: MutableByteArray s -> Int -> Int -> MutableByteArray s -> Int -> ST s Double
dotAccum mbaY bs c mbaZL nL = go 0 0
  where
    go !l !acc
      | l >= nL = pure acc
      | otherwise = do
          ycl <- readRawD mbaY 0 (c * bs + l)
          zl  <- readRawD mbaZL 0 l
          go (l+1) (acc + ycl * zl)

-- Helper: sum_l V_R[c,l]*zX[l] for l = 0..nR-1
dotVR_zX :: MutableByteArray s -> Int -> Int -> Int -> MutableByteArray s
         -> Int -> Int -> ST s Double
dotVR_zX mbaA offA nn k0 mbaZX c nR = go 0 0
  where
    go !l !acc
      | l >= nR = pure acc
      | otherwise = do
          let !kl = k0 + l
          vr <- if c == kl + 1 then pure 1.0
                else if c > kl + 1 then readRawD mbaA offA (kl * nn + c)
                else pure 0.0
          zx <- readRawD mbaZX 0 l
          go (l+1) (acc + vr * zx)

-- Helper: Y[:,l]^T * u where u = [1, A[k, k+2:n-1]]
dotY_u :: MutableByteArray s -> Int -> Int -> MutableByteArray s -> Int
       -> Int -> Int -> ST s Double
dotY_u mbaA offA nn mbaY bs k l = do
  yk1l <- readRawD mbaY 0 ((k+1) * bs + l)
  go (k+2) yk1l
  where
    go !c !acc
      | c >= nn = pure acc
      | otherwise = do
          ycl <- readRawD mbaY 0 (c * bs + l)
          uc  <- readRawD mbaA offA (k * nn + c)
          go (c+1) (acc + ycl * uc)

-- Helper: V_R[:,l]^T * u where u = [1, A[k, k+2:n-1]]
dotVR_u :: MutableByteArray s -> Int -> Int -> Int -> Int -> Int -> ST s Double
dotVR_u mbaA offA nn k0 k l = do
  -- u[c-k-1]: u[0]=1 at c=k+1, u[c-k-1]=A[k,c] for c>k+1
  -- V_R[c,l]: 1 if c==kl+1, A[kl,c] if c>kl+1, 0 if c<=kl
  -- We need sum_{c=k+1}^{n-1} V_R[c,l] * u[c-k-1]
  -- Since k > kl (k=k0+j, l<j), V_R[k+1,l] = A[kl, k+1] (if k+1 > kl+1, i.e., k > kl)
  vrkp1 <- if k + 1 == kl + 1 then pure 1.0
            else readRawD mbaA offA (kl * nn + (k+1))
  go (k+2) vrkp1
  where
    !kl = k0 + l
    go !c !acc
      | c >= nn = pure acc
      | otherwise = do
          vrc <- if c == kl + 1 then pure 1.0
                 else readRawD mbaA offA (kl * nn + c)
          uc  <- readRawD mbaA offA (k * nn + c)
          go (c+1) (acc + vrc * uc)

-- Helper: A * u at row i, where u = [1, A[k, k+2:n-1]]
dotA_u :: MutableByteArray s -> Int -> Int -> Int -> Int -> ST s Double
dotA_u mbaA offA nn k i = do
  aikp1 <- readRawD mbaA offA (i * nn + (k+1))
  go (k+2) aikp1
  where
    go !c !acc
      | c >= nn = pure acc
      | otherwise = do
          aic <- readRawD mbaA offA (i * nn + c)
          uc  <- readRawD mbaA offA (k * nn + c)
          go (c+1) (acc + aic * uc)

-- Helper: sum_l V_L[i,l]*zL[l] for l = 0..nL-1
dotVL_zL :: MutableByteArray s -> Int -> Int -> Int -> MutableByteArray s
         -> Int -> Int -> ST s Double
dotVL_zL mbaA offA nn k0 mbaZL i nL = go 0 0
  where
    go !l !acc
      | l >= nL = pure acc
      | otherwise = do
          let !kl = k0 + l
          vl <- if i == kl then pure 1.0
                else if i > kl then readRawD mbaA offA (i * nn + kl)
                else pure 0.0
          zl <- readRawD mbaZL 0 l
          go (l+1) (acc + vl * zl)

-- Helper: sum_l X[i,l]*zX[l] for l = 0..nR-1
dotX_zX :: MutableByteArray s -> Int -> MutableByteArray s -> Int -> Int -> ST s Double
dotX_zX mbaX bs mbaZX i nR = go 0 0
  where
    go !l !acc
      | l >= nR = pure acc
      | otherwise = do
          xil <- readRawD mbaX 0 (i * bs + l)
          zx  <- readRawD mbaZX 0 l
          go (l+1) (acc + xil * zx)

-- | Apply trailing GEMM update after a panel step.
-- A[k0+bs:m, k0+bs:n] -= V_L_trail * Y_trail^T + X_trail * V_R_trail^T
applyTrailingUpdate :: MutableByteArray s -> Int -> Int -> Int
                    -> MutableByteArray s -> MutableByteArray s
                    -> MutableByteArray s -> MutableByteArray s
                    -> MutableByteArray s -> MutableByteArray s
                    -> MutableByteArray s
                    -> Int -> Int -> Int -> Int -> ST s ()
applyTrailingUpdate mbaA offA _mm nn mbaX mbaY
                    mbaVLbuf mbaYTbuf mbaXbuf mbaVRTbuf mbaTrail
                    k0 bs remR remC = do
  let !trailRowStart = k0 + bs
      !trailColStart = k0 + bs

  -- Copy A_trail to contiguous buffer (remR × remC)
  forM_ [0..remR-1] $ \i ->
    forM_ [0..remC-1] $ \c -> do
      val <- readRawD mbaA offA ((trailRowStart + i) * nn + trailColStart + c)
      writeRawD mbaTrail 0 (i * remC + c) val

  -- Build V_L_trail (remR × bs): V_L[trailRowStart+i, l] for i=0..remR-1, l=0..bs-1
  -- For all trail rows, i >= trailRowStart > k0+l, so V_L[i,l] = A[i, k0+l]
  forM_ [0..remR-1] $ \i ->
    forM_ [0..bs-1] $ \l ->
      readRawD mbaA offA ((trailRowStart + i) * nn + (k0 + l)) >>=
        writeRawD mbaVLbuf 0 (i * bs + l)

  -- Build Y_trail^T (bs × remC): Y_trail^T[l, c] = Y[trailColStart+c, l]
  forM_ [0..bs-1] $ \l ->
    forM_ [0..remC-1] $ \c ->
      readRawD mbaY 0 ((trailColStart + c) * bs + l) >>=
        writeRawD mbaYTbuf 0 (l * remC + c)

  -- GEMM 1: trail -= V_L_trail * Y_trail^T
  -- Negate V_L_trail: nVL = -V_L_trail
  rawNegateDoubles mbaVLbuf 0 (remR * bs)
  baVL <- unsafeFreezeByteArray mbaVLbuf
  baYT <- unsafeFreezeByteArray mbaYTbuf
  rawGemmKernel baVL 0 baYT 0 mbaTrail 0 remR bs remC

  -- Build X_trail (remR × bs): X[trailRowStart+i, l]
  forM_ [0..remR-1] $ \i ->
    forM_ [0..bs-1] $ \l ->
      readRawD mbaX 0 ((trailRowStart + i) * bs + l) >>=
        writeRawD mbaXbuf 0 (i * bs + l)

  -- Build V_R_trail^T (bs × remC): V_R_trail^T[l, c] = V_R[trailColStart+c, l]
  -- V_R[c, l] = 1 if c==k0+l+1, A[k0+l, c] if c>k0+l+1, 0 if c<=k0+l
  -- Must handle implicit 1: when trailColStart+c == k0+l+1 (i.e., l=bs-1, c=0)
  forM_ [0..bs-1] $ \l ->
    forM_ [0..remC-1] $ \c -> do
      let !globalC = trailColStart + c
          !kl = k0 + l
      val <- if globalC == kl + 1 then pure 1.0
             else if globalC > kl + 1 then readRawD mbaA offA (kl * nn + globalC)
             else pure 0.0
      writeRawD mbaVRTbuf 0 (l * remC + c) val

  -- GEMM 2: trail -= X_trail * V_R_trail^T
  rawNegateDoubles mbaXbuf 0 (remR * bs)
  baX  <- unsafeFreezeByteArray mbaXbuf
  baVRT <- unsafeFreezeByteArray mbaVRTbuf
  rawGemmKernel baX 0 baVRT 0 mbaTrail 0 remR bs remC

  -- Copy trail back to A
  forM_ [0..remR-1] $ \i ->
    forM_ [0..remC-1] $ \c -> do
      val <- readRawD mbaTrail 0 (i * remC + c)
      writeRawD mbaA offA ((trailRowStart + i) * nn + trailColStart + c) val

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
