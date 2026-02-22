{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE MagicHash #-}
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
import Data.Ord (Down)
import GHC.Exts
import GHC.ST (ST(..))

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.BLAS.Level3 (matMul, matMulP, transpose, matMulAtAP)
-- matvecP no longer needed: U-matrix now computed via single GEMM
import Numeric.LinearAlgebra.Massiv.Eigen.Symmetric (symmetricEigen, symmetricEigenP)
import Numeric.LinearAlgebra.Massiv.Internal.Kernel
  ( rawMutSumSqColumn, rawMutSumSqRow
  , rawMutHouseholderApply, rawMutHouseholderApplyRow
  , rawMutQAccum
  , rawMutApplyGivensColumns )

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
      -- Eigendecomposition of AᵀA using raw primop QR iteration
      (!eigvalsRaw, !vRaw) = symmetricEigenP ata (max 30 (6 * nn)) 1e-12
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

  -- Accumulate left Householder reflectors into U (backward)
  -- Left reflector k: v stored in column k of A, rows k+1..m-1, with v[k]=1 implicit
  forM_ (reverse [0..nn-1]) $ \k -> do
    betaK <- readRawD mbaBetaL 0 k
    when (betaK /= 0) $
      forM_ [0..mm-1] $ \row ->
        rawMutQAccum mbaU offU mm baA offFA nn betaK k mm row

  -- V₀ = G₁ G₂ ... G_{n-3} (right reflectors)
  mV <- M.newMArray @M.P (Sz (nn :. nn)) (0 :: Double)
  let mbaV = unwrapMutableByteArray mV
      offV = unwrapMutableByteArrayOffset mV
  -- Initialise V = I
  forM_ [0..nn-1] $ \i ->
    writeRawD mbaV offV (i * nn + i) 1.0

  -- Accumulate right Householder reflectors into V (backward)
  -- Right reflector k: v stored in row k of A, cols k+2..n-1, with v[k+1]=1 implicit
  when (nn >= 3) $
    forM_ (reverse [0..nn-3]) $ \k -> do
      betaK <- readRawD mbaBetaR 0 k
      when (betaK /= 0) $
        forM_ [0..nn-1] $ \row ->
          rightQAccum mbaV offV nn baA offFA nn betaK k nn row

  -- Phase 3: Implicit-shift bidiagonal QR iteration
  bidiagQRIterP mbaD 0 mbaE 0 mbaU offU mm mbaV offV nn nn (30 * nn)

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

      -- Also handle e[k-1] if k > lo (it gets mixed by previous left rotation)
      -- The bulge from the right rotation is at position (k+1, k), i.e. below the diagonal

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
-- [c s; -s c]^T [a; b] = [r; 0].
givens :: Double -> Double -> (Double, Double)
givens a b
  | b == 0    = (1, 0)
  | abs b > abs a =
      let tau = -a / b
          s   = 1 / sqrt (1 + tau * tau)
          c   = s * tau
      in (c, s)
  | otherwise =
      let tau = -b / a
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
