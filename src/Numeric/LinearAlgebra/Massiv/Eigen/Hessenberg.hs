{-# LANGUAGE AllowAmbiguousTypes #-}

-- |
-- Module      : Numeric.LinearAlgebra.Massiv.Eigen.Hessenberg
-- Copyright   : (c) Nadia Chambers 2026
-- License     : BSD-3-Clause
-- Maintainer  : nadia.chambers@iohk.io
-- Stability   : experimental
--
-- Reduction of a general square matrix to upper Hessenberg form via
-- Householder similarity transformations, following Golub & Van Loan,
-- /Matrix Computations/, 4th edition (GVL4), Section 7.4, pp. 383--393.
--
-- An upper Hessenberg matrix \(H\) satisfies \(h_{ij} = 0\) for
-- \(i > j + 1\); that is, all entries below the first subdiagonal are zero.
-- Given \(A \in \mathbb{R}^{n \times n}\), the algorithm computes an
-- orthogonal matrix \(Q\) (accumulated as a product of \(n - 2\) Householder
-- reflectors) such that
--
-- \[
--   A = Q H Q^T
-- \]
--
-- This is the standard first phase in eigenvalue algorithms (e.g. the
-- implicit QR algorithm in "Numeric.LinearAlgebra.Massiv.Eigen.Schur"),
-- because QR steps preserve Hessenberg structure and are much cheaper on a
-- Hessenberg matrix than on a full matrix.
--
-- __Algorithm:__ Householder reduction to Hessenberg form (GVL4
-- Algorithm 7.4.2, p. 387).
--
-- __Complexity:__ \(\frac{10}{3} n^3\) floating-point operations.
module Numeric.LinearAlgebra.Massiv.Eigen.Hessenberg
  ( -- * Hessenberg reduction (Algorithm 7.4.2)
    hessenberg
  , hessenbergInPlace
  ) where

import qualified Data.Massiv.Array as M
import Data.Massiv.Array (Ix1, Ix2(..), Sz(..))
import GHC.TypeNats (KnownNat)

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal

-- | Reduce a general matrix to upper Hessenberg form (GVL4 Algorithm 7.4.2,
-- p. 387).
--
-- Computes the factorisation
--
-- \[
--   A = Q \, H \, Q^T
-- \]
--
-- where \(Q\) is orthogonal (a product of \(n - 2\) Householder reflectors)
-- and \(H\) is upper Hessenberg, i.e. \(h_{ij} = 0\) for \(i > j + 1\).
--
-- At step \(k\) the algorithm determines a Householder reflector
-- \(P_k = I - \beta_k v_k v_k^T\) that zeroes entries \(k+2, \ldots, n\) in
-- column \(k\) of the current matrix, and applies it as a similarity
-- transformation \(H \leftarrow P_k H P_k\).
--
-- __Complexity:__ \(\frac{10}{3} n^3\) flops (GVL4, p. 389).
--
-- Returns @(Q, H)@.
hessenberg :: forall n r e. (KnownNat n, M.Manifest r e, Floating e, Ord e)
           => Matrix n n r e -> (Matrix n n r e, Matrix n n r e)
hessenberg a =
  let nn = dimVal @n
      go :: Int -> Matrix n n r e -> Matrix n n r e -> (Matrix n n r e, Matrix n n r e)
      go k q_ h_
        | k >= nn - 2 = (q_, h_)
        | otherwise =
          -- Compute Householder to zero out h(k+2:n, k)
          let x0 = h_ ! (k+1, k)
              sigma = foldl' (\acc i -> acc + (h_ ! (i, k)) * (h_ ! (i, k))) 0 [k+2..nn-1]
          in if sigma == 0 && x0 >= 0
             then go (k + 1) q_ h_
             else
              let mu = sqrt (x0 * x0 + sigma)
                  v0 = if x0 <= 0 then x0 - mu else -sigma / (x0 + mu)
                  beta = 2 * v0 * v0 / (sigma + v0 * v0)
                  -- v is zero for indices 0..k, v(k+1)=1, v(i) = h(i,k)/v0 for i > k+1
                  v = makeVector @n @r $ \i ->
                    if i <= k then 0
                    else if i == k + 1 then 1
                    else (h_ ! (i, k)) / v0
                  -- H ← (I - β·v·vᵀ)·H·(I - β·v·vᵀ) = similarity transform
                  -- First: H ← (I - β·v·vᵀ)·H
                  h1 = applyFromLeft v beta h_
                  -- Then: H ← H·(I - β·v·vᵀ)
                  h2 = applyFromRight h1 v beta
                  -- Q ← Q·(I - β·v·vᵀ)
                  q_new = applyFromRight q_ v beta
              in go (k + 1) q_new h2

      q0 = identityMatrix @n @r
  in go 0 q0 a

-- Helper: apply Householder from left
applyFromLeft :: forall n r e. (KnownNat n, M.Manifest r e, Num e)
              => Vector n r e -> e -> Matrix n n r e -> Matrix n n r e
applyFromLeft v beta h =
  let nn = dimVal @n
  in makeMatrix @n @n @r $ \i j ->
    let wj = beta * foldl' (\acc k -> acc + (v !. k) * (h ! (k, j))) 0 [0..nn-1]
    in (h ! (i, j)) - (v !. i) * wj

-- Helper: apply Householder from right
applyFromRight :: forall n r e. (KnownNat n, M.Manifest r e, Num e)
               => Matrix n n r e -> Vector n r e -> e -> Matrix n n r e
applyFromRight h v beta =
  let nn = dimVal @n
  in makeMatrix @n @n @r $ \i j ->
    let wi = beta * foldl' (\acc k -> acc + (h ! (i, k)) * (v !. k)) 0 [0..nn-1]
    in (h ! (i, j)) - wi * (v !. j)

-- | Compute only the upper Hessenberg matrix \(H\), discarding the
-- orthogonal factor \(Q\).
--
-- This is a convenience wrapper around 'hessenberg' that returns only the
-- second component of the pair.  Use this when only the Hessenberg form is
-- needed and the transformation matrix is not required (e.g. when computing
-- eigenvalues but not eigenvectors).
hessenbergInPlace :: forall n r e. (KnownNat n, M.Manifest r e, Floating e, Ord e)
                  => Matrix n n r e -> Matrix n n r e
hessenbergInPlace a = snd (hessenberg a)
