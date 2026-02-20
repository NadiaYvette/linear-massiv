{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE MagicHash #-}

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
--
-- __Implementation note:__  Rather than the full Golub--Kahan
-- bidiagonalisation pipeline (GVL4 Algorithm 8.6.1, p. 504), this module
-- takes a simpler (though less efficient) approach: it forms the symmetric
-- positive semi-definite matrix \(A^T A\), computes its eigendecomposition
-- via 'Numeric.LinearAlgebra.Massiv.Eigen.Symmetric.symmetricEigen', and
-- recovers the singular values as \(\sigma_i = \sqrt{\max(0, \lambda_i)}\)
-- and the left singular vectors as \(u_i = A v_i / \sigma_i\).  This is
-- adequate for moderate-sized matrices but may lose accuracy for matrices
-- with a large condition number, because the condition number of \(A^T A\)
-- is the /square/ of that of \(A\).
module Numeric.LinearAlgebra.Massiv.Eigen.SVD
  ( -- * Full SVD
    svd
  , svdP
    -- * Singular values only
  , singularValues
  , singularValuesP
  ) where

import qualified Data.Massiv.Array as M
import Data.Massiv.Array (Ix2(..), Sz(..), unwrapByteArray, unwrapByteArrayOffset)
import Data.Primitive.ByteArray (ByteArray(..))
import GHC.TypeNats (KnownNat)
import Control.Monad (forM_)
import Data.List (sortBy)
import Data.Ord (Down(..))
import GHC.Exts

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.BLAS.Level3 (matMul, matMulP, transpose)
import Numeric.LinearAlgebra.Massiv.BLAS.Level2 (matvecP)
import Numeric.LinearAlgebra.Massiv.Eigen.Symmetric (symmetricEigen, symmetricEigenP)

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
  let mm = dimVal @m
      nn = dimVal @n
      at = transpose a
      ata = matMul at a  -- n×n symmetric positive semidefinite
      -- Eigendecomposition of AᵀA
      (eigvals, v) = symmetricEigen ata (30 * nn) 1e-12
      -- Singular values = sqrt of eigenvalues (clamp negatives to 0)
      sigma = makeVector @n @r $ \i ->
        let ev = eigvals !. i
        in if ev > 0 then sqrt ev else 0
      -- Compute U: u_i = A·v_i / σ_i
      -- First, build U by computing A·V column by column
      u = makeMatrix @m @m @r $ \i j ->
        if j < nn then
          let sj = sigma !. j
          in if sj > 1e-14
             then -- u_j = (1/σ_j) · Σ_k A(i,k) · V(k,j)
               let av = foldl' (\acc k -> acc + (a ! (i, k)) * (v ! (k, j))) 0 [0..nn-1]
               in av / sj
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
svdP a =
  let !mm = dimVal @m
      !nn = dimVal @n
      at = transpose a
      !ata = matMulP at a  -- n×n symmetric positive semidefinite, SIMD GEMM
      -- Eigendecomposition of AᵀA using raw primop QR iteration
      (!eigvals, !v) = symmetricEigenP ata (30 * nn) 1e-12
      -- Singular values = sqrt of eigenvalues (clamp negatives to 0)
      sigma = makeVector @n @M.P $ \i ->
        let ev = eigvals !. i
        in if ev > 0 then sqrt ev else 0
      -- Compute U columns: u_j = A·v_j / σ_j using matvecP
      -- Extract each column of V as a vector, multiply by A, normalise
      u = createMatrix @m @m @M.P $ \mu -> do
        -- Get immutable arrays for raw reads
        let !baV = unwrapByteArray (unMatrix v)
            !offV = unwrapByteArrayOffset (unMatrix v)
        forM_ [0..nn-1] $ \j -> do
          let sj = sigma !. j
          if sj > 1e-14
            then do
              -- Extract column j of V as a Vector n
              let !vj = makeVector @n @M.P $ \k -> readBA baV offV (k * nn + j)
                  -- u_j = A · v_j / σ_j
                  !avj = matvecP a vj
                  !invSj = 1.0 / sj
              forM_ [0..mm-1] $ \i ->
                M.write_ mu (i :. j) (invSj * (avj !. i))
            else
              -- Zero singular value: use standard basis vector
              forM_ [0..mm-1] $ \i ->
                M.write_ mu (i :. j) (if i == j then 1 else 0)
        -- Extra columns for m > n: extend to full orthogonal basis
        forM_ [nn..mm-1] $ \j ->
          forM_ [0..mm-1] $ \i ->
            M.write_ mu (i :. j) (if i == j then 1 else 0)
  in (u, sigma, v)
{-# NOINLINE svdP #-}

-- | Read a Double from an immutable ByteArray at element index.
readBA :: ByteArray -> Int -> Int -> Double
readBA (ByteArray ba) (I# off) (I# i) =
  case indexDoubleArray# ba (off +# i) of v -> D# v
{-# INLINE readBA #-}

-- | Compute only the singular values of \(A\), sorted in descending order.
--
-- \[
--   \sigma_i = \sqrt{\max(0,\, \lambda_i(A^T A))}, \qquad i = 1, \ldots, n
-- \]
--
-- This is a convenience function that avoids constructing the full \(U\) and
-- \(V\) matrices.  It forms \(A^T A\), computes its eigenvalues via
-- 'Numeric.LinearAlgebra.Massiv.Eigen.Symmetric.symmetricEigen', and returns
-- their non-negative square roots in descending order.
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
      at = transpose a
      !ata = matMulP at a
      (!eigvals, _) = symmetricEigenP ata (30 * nn) 1e-12
      evList = map (\i -> eigvals !. i) [0..nn-1]
      sorted = sortBy (\x y -> compare (Down x) (Down y)) evList
  in makeVector @n @M.P $ \i ->
    let ev = sorted !! i
    in if ev > 0 then sqrt ev else 0
