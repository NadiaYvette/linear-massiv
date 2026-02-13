{-# LANGUAGE AllowAmbiguousTypes #-}

-- |
-- Module      : Numeric.LinearAlgebra.Massiv.Norms
-- Copyright   : (c) Nadia Chambers 2026
-- License     : BSD-3-Clause
-- Maintainer  : nadia.chambers@iohk.io
-- Stability   : experimental
--
-- = Matrix and Vector Norms
--
-- Norms measure the "size" of vectors and matrices and are fundamental to
-- error analysis, convergence criteria, and conditioning estimates in
-- numerical linear algebra.
--
-- This module implements the norms described in:
--
-- * Golub, G. H. & Van Loan, C. F. (2013). /Matrix Computations/, 4th ed.
--   __Chapter 2: Matrix Analysis__, Sections 2.3–2.7, pp. 71–95.
--
-- == Vector norms (GVL4 Section 2.3, p. 71)
--
-- For a vector \(x \in \mathbb{R}^n\):
--
-- * 1-norm: \(\|x\|_1 = \sum_{i=1}^{n} |x_i|\) — see 'vnorm1'
-- * 2-norm (Euclidean): \(\|x\|_2 = \sqrt{\sum_{i=1}^{n} x_i^2}\) — see 'vnorm2'
-- * \(\infty\)-norm: \(\|x\|_\infty = \max_i |x_i|\) — see 'vnormInf'
--
-- == Matrix norms (GVL4 Section 2.3, pp. 71–78)
--
-- For a matrix \(A \in \mathbb{R}^{m \times n}\):
--
-- * Frobenius norm: \(\|A\|_F = \sqrt{\sum_{i,j} a_{ij}^2} = \sqrt{\text{trace}(A^T A)}\)
--   — see 'normFrob'
-- * 1-norm (max column sum): \(\|A\|_1 = \max_j \sum_i |a_{ij}|\) — see 'norm1'
-- * \(\infty\)-norm (max row sum): \(\|A\|_\infty = \max_i \sum_j |a_{ij}|\) — see 'normInf'
--
-- These satisfy the norm axioms: non-negativity, homogeneity, and the
-- triangle inequality \(\|A + B\| \leq \|A\| + \|B\|\).
--
-- == Condition numbers (GVL4 Section 2.7, pp. 87–95)
--
-- The condition number \(\kappa(A) = \|A\| \cdot \|A^{-1}\|\) measures
-- how sensitive the solution of \(Ax = b\) is to perturbations in \(A\)
-- and \(b\). See 'condFrob' for a placeholder using the Frobenius norm.
module Numeric.LinearAlgebra.Massiv.Norms
  ( -- * Vector norms (GVL4 Section 2.3, p. 71)
    vnorm1
  , vnorm2
  , vnormInf
    -- * Matrix norms (GVL4 Section 2.3, pp. 71–78)
  , normFrob
  , norm1
  , normInf
    -- * Condition number estimate (GVL4 Section 2.7, p. 87)
  , condFrob
  ) where

import qualified Data.Massiv.Array as M
import Data.Massiv.Array (Ix2(..), Sz(..))
import GHC.TypeNats (KnownNat)

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal

-- | Vector 1-norm (GVL4 p. 71, eq. 2.3.1).
--
-- \[
--   \|x\|_1 = \sum_{i=1}^{n} |x_i|
-- \]
--
-- Complexity: \(O(n)\).
vnorm1 :: (KnownNat n, M.Manifest r e, Num e, Ord e)
        => Vector n r e -> e
vnorm1 (MkVector arr) = M.foldlS (\acc x -> acc + abs x) 0 arr

-- | Vector 2-norm, the Euclidean norm (GVL4 p. 71, eq. 2.3.2).
--
-- \[
--   \|x\|_2 = \sqrt{\sum_{i=1}^{n} x_i^2} = \sqrt{x^T x}
-- \]
--
-- Complexity: \(O(n)\).
vnorm2 :: (KnownNat n, M.Manifest r e, Floating e)
        => Vector n r e -> e
vnorm2 (MkVector arr) = sqrt $ M.foldlS (\acc x -> acc + x * x) 0 arr

-- | Vector \(\infty\)-norm (GVL4 p. 71, eq. 2.3.3).
--
-- \[
--   \|x\|_\infty = \max_{1 \leq i \leq n} |x_i|
-- \]
--
-- Complexity: \(O(n)\).
vnormInf :: (KnownNat n, M.Manifest r e, Num e, Ord e)
          => Vector n r e -> e
vnormInf (MkVector arr) = M.foldlS (\acc x -> max acc (abs x)) 0 arr

-- | Frobenius norm (GVL4 p. 72, eq. 2.3.7).
--
-- \[
--   \|A\|_F = \sqrt{\sum_{i=1}^{m} \sum_{j=1}^{n} a_{ij}^2}
--           = \sqrt{\text{trace}(A^T A)}
-- \]
--
-- The Frobenius norm is /not/ an operator norm (it is not subordinate
-- to any vector norm), but it is submultiplicative:
-- \(\|AB\|_F \leq \|A\|_F \|B\|_F\).
--
-- Complexity: \(O(mn)\).
normFrob :: (KnownNat m, KnownNat n, M.Manifest r e, Floating e)
          => Matrix m n r e -> e
normFrob (MkMatrix arr) = sqrt $ M.foldlS (\acc x -> acc + x * x) 0 arr

-- | Matrix 1-norm — maximum absolute column sum (GVL4 p. 72, eq. 2.3.10).
--
-- \[
--   \|A\|_1 = \max_{1 \leq j \leq n} \sum_{i=1}^{m} |a_{ij}|
-- \]
--
-- This is the operator norm subordinate to the vector 1-norm:
-- \(\|A\|_1 = \max_{\|x\|_1 = 1} \|Ax\|_1\).
--
-- Complexity: \(O(mn)\).
norm1 :: forall m n r e. (KnownNat m, KnownNat n, M.Manifest r e, Num e, Ord e)
      => Matrix m n r e -> e
norm1 mat =
  let c = dimVal @n
      r = dimVal @m
      colSum j = foldl (\acc i -> acc + abs (mat ! (i, j))) 0 [0..r-1]
  in maximum $ map colSum [0..c-1]

-- | Matrix \(\infty\)-norm — maximum absolute row sum (GVL4 p. 72, eq. 2.3.11).
--
-- \[
--   \|A\|_\infty = \max_{1 \leq i \leq m} \sum_{j=1}^{n} |a_{ij}|
-- \]
--
-- This is the operator norm subordinate to the vector \(\infty\)-norm.
-- Note that \(\|A\|_\infty = \|A^T\|_1\).
--
-- Complexity: \(O(mn)\).
normInf :: forall m n r e. (KnownNat m, KnownNat n, M.Manifest r e, Num e, Ord e)
        => Matrix m n r e -> e
normInf mat =
  let c = dimVal @n
      r = dimVal @m
      rowSum i = foldl (\acc j -> acc + abs (mat ! (i, j))) 0 [0..c-1]
  in maximum $ map rowSum [0..r-1]

-- | Estimate condition number using the Frobenius norm (GVL4 Section 2.7, p. 87).
--
-- The condition number is defined as
-- \(\kappa_F(A) = \|A\|_F \cdot \|A^{-1}\|_F\).
--
-- __Note__: This function currently returns only \(\|A\|_F\) as a placeholder.
-- Computing \(\|A^{-1}\|_F\) requires solving a linear system (e.g., via LU),
-- introducing a circular dependency. Users should compute the full condition
-- number by combining 'normFrob' with an explicit inverse or using SVD-based
-- estimates (\(\kappa_2 = \sigma_{\max} / \sigma_{\min}\)).
condFrob :: (KnownNat n, M.Manifest r e, Floating e)
         => Matrix n n r e -> e
condFrob = normFrob
