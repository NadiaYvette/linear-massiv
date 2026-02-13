{-# LANGUAGE AllowAmbiguousTypes #-}

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
--   This is the workhorse algorithm for dense QR factorisation.
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
module Numeric.LinearAlgebra.Massiv.Orthogonal.QR
  ( -- * QR factorisation (Householder)
    qr
  , qrR
    -- * QR factorisation (Givens)
  , qrGivens
  ) where

import qualified Data.Massiv.Array as M
import Data.Massiv.Array (Ix1, Ix2(..), Sz(..))
import GHC.TypeNats (KnownNat, natVal)
import Data.Proxy (Proxy(..))

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.Orthogonal.Householder
  (householderVector, applyHouseholderLeft, householderMatrix)
import Numeric.LinearAlgebra.Massiv.Orthogonal.Givens
  (givensRotation, applyGivensLeft)
import Numeric.LinearAlgebra.Massiv.BLAS.Level3 (matMul)

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
-- The algorithm applies \( \min(m, n) \) Householder reflections from the
-- left to zero out the sub-diagonal entries of each column in turn.  The
-- orthogonal factor \( Q \) is accumulated explicitly as the product
-- \( Q = H_1 \, H_2 \cdots H_k \), where each
-- \( H_j = I - \beta_j \, v_j \, v_j^T \).
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

      -- Iteratively apply Householder reflections
      -- We accumulate Q = H₁·H₂·...·Hₖ and update R
      go :: Int -> Matrix m m r e -> Matrix m n r e -> (Matrix m m r e, Matrix m n r e)
      go k q_ r_
        | k >= steps = (q_, r_)
        | otherwise =
          -- Extract column k of R from row k downwards
          let colLen = mm - k
              -- Build the sub-vector x = R(k:m, k)
              x_full = makeVector @m @r $ \i ->
                if i < colLen then r_ ! (i + k, k) else 0
              -- We need to compute Householder of the subvector
              -- But since our types are fixed, we work with full-size reflectors
              -- Compute sigma and mu for the subcolumn
              x0 = r_ ! (k, k)
              sigma = foldl' (\acc i -> acc + (r_ ! (i, k)) * (r_ ! (i, k))) 0 [k+1..mm-1]
          in if sigma == 0 && x0 >= 0
             then go (k + 1) q_ r_  -- Already in desired form
             else
              let mu = sqrt (x0 * x0 + sigma)
                  v0 = if x0 <= 0 then x0 - mu else -sigma / (x0 + mu)
                  beta = 2 * v0 * v0 / (sigma + v0 * v0)
                  -- Full-size Householder vector: v(i) = 0 for i < k, v(k) = 1, v(i) = x(i)/v0 for i > k
                  v = makeVector @m @r $ \i ->
                    if i < k then 0
                    else if i == k then 1
                    else (r_ ! (i, k)) / v0
                  -- Apply H = I - β·v·vᵀ to R from the left
                  r_new = applyHouseholderLeftRect @m @n v beta r_
                  -- Accumulate Q: Q ← Q·H (from the right)
                  -- H is symmetric so Q·H = Q·(I - β·v·vᵀ)
                  q_new = applyHouseholderRightQ @m v beta q_
              in go (k + 1) q_new r_new

      q0 = identityMatrix @m @r
  in go 0 q0 a

-- Helper: apply Householder from left to m×n matrix
applyHouseholderLeftRect :: forall m n r e. (KnownNat m, KnownNat n, M.Manifest r e, Num e)
                         => Vector m r e -> e -> Matrix m n r e -> Matrix m n r e
applyHouseholderLeftRect v beta a =
  let mm = dimVal @m
      nn = dimVal @n
  in makeMatrix @m @n @r $ \i j ->
    let wj = beta * foldl' (\acc k -> acc + (v !. k) * (a ! (k, j))) 0 [0..mm-1]
    in (a ! (i, j)) - (v !. i) * wj

-- Helper: apply Householder from right to m×m matrix (for Q accumulation)
applyHouseholderRightQ :: forall m r e. (KnownNat m, M.Manifest r e, Num e)
                       => Vector m r e -> e -> Matrix m m r e -> Matrix m m r e
applyHouseholderRightQ v beta q =
  let mm = dimVal @m
  in makeMatrix @m @m @r $ \i j ->
    let wi = beta * foldl' (\acc k -> acc + (q ! (i, k)) * (v !. k)) 0 [0..mm-1]
    in (q ! (i, j)) - wi * (v !. j)

-- | Compute only the upper triangular factor \( R \) from the QR
-- factorisation, without explicitly forming \( Q \)
-- (GVL4 Algorithm 5.2.1, p. 249).
--
-- This is a convenience wrapper around 'qr' that discards the
-- orthogonal factor.  In a production implementation the Householder
-- vectors would be stored in the sub-diagonal part of the result and
-- \( Q \) would never be formed, saving \( O(m^2 n) \) flops.
-- The current implementation computes the full QR and returns only
-- \( R \).
--
-- __Complexity:__ Same as 'qr'.
qrR :: forall m n r e. (KnownNat m, KnownNat n, M.Manifest r e, Floating e, Ord e)
    => Matrix m n r e -> Matrix m n r e
qrR a = snd (qr a)

-- | QR factorisation via Givens rotations (GVL4 Section 5.2.4, p. 252).
--
-- Given \( A \in \mathbb{R}^{m \times n} \) with \( m \ge n \), compute
-- the factorisation \( A = Q \, R \) by applying a sequence of Givens
-- rotations to zero out sub-diagonal entries one at a time.
--
-- For each column \( j \), the sub-diagonal entries
-- \( A(j+1, j), A(j+2, j), \ldots, A(m-1, j) \) are zeroed by
-- rotations in the \( (j, i) \) plane.  The orthogonal factor
-- \( Q \) is accumulated as the product of all applied rotations.
--
-- This variant is particularly useful for:
--
-- * Sparse matrices, where Householder reflections would destroy
--   sparsity structure.
-- * Upper Hessenberg matrices, where only a single sub-diagonal entry
--   per column needs zeroing, giving an \( O(mn) \) algorithm.
-- * Situations requiring incremental updates to an existing QR
--   factorisation.
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

      go' :: Int -> Matrix m m r e -> Matrix m n r e -> (Matrix m m r e, Matrix m n r e)
      go' j q_ r_
        | j >= steps = (q_, r_)
        | otherwise =
          let (q_final, r_final) = foldl
                (\(qq, rr) i ->
                  let aij = rr ! (i, j)
                  in if aij == 0 then (qq, rr)
                     else let (c, s) = givensRotation (rr ! (j, j)) aij
                              rr' = applyGivensLeft c s j i rr
                              qq' = applyGivensRightQ c s j i qq
                          in (qq', rr')
                ) (q_, r_) [j+1..mm-1]
          in go' (j + 1) q_final r_final

      q0 = identityMatrix @m @r
  in go' 0 q0 a

-- Helper: apply Givens from right to square matrix (for Q accumulation)
applyGivensRightQ :: forall m r e. (KnownNat m, M.Manifest r e, Num e)
                  => e -> e -> Int -> Int -> Matrix m m r e -> Matrix m m r e
applyGivensRightQ c s ci ck q =
  makeMatrix @m @m @r $ \i j ->
    if j == ci then
      c * (q ! (i, ci)) - s * (q ! (i, ck))
    else if j == ck then
      s * (q ! (i, ci)) + c * (q ! (i, ck))
    else
      q ! (i, j)
