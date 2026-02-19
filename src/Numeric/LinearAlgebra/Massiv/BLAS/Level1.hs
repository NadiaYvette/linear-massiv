{-# LANGUAGE AllowAmbiguousTypes #-}

-- |
-- Module      : Numeric.LinearAlgebra.Massiv.BLAS.Level1
-- Copyright   : (c) Nadia Chambers 2026
-- License     : BSD-3-Clause
-- Maintainer  : nadia.chambers@iohk.io
-- Stability   : experimental
--
-- = BLAS Level 1: Vector–Vector Operations
--
-- This module provides type-safe, dimension-indexed wrappers around the
-- standard BLAS Level 1 kernels for vector–vector operations.  Every
-- function carries the vector length @n@ as a phantom type-level natural,
-- so dimension mismatches are caught at compile time rather than at run
-- time.
--
-- The algorithms implemented here correspond to the elementary building
-- blocks described in:
--
-- * Golub, G. H. & Van Loan, C. F. (2013). /Matrix Computations/,
--   4th edition (GVL4). Johns Hopkins University Press.
--   __Chapter 1, Section 1.1__, pp. 4–8.
--
-- Specifically:
--
-- * __Algorithm 1.1.1__ (p. 4) — Inner product (dot product).
--   Given vectors \(x, y \in \mathbb{R}^{n}\), compute
--   \(c = x^{T} y = \sum_{i=1}^{n} x_i y_i\).
--
-- * __Algorithm 1.1.2 (Saxpy)__ (p. 4) — Scalar \(\alpha\) times
--   vector \(x\) plus vector \(y\):
--   \(y \leftarrow \alpha x + y\).
--   This is the fundamental vector-update operation upon which the
--   higher-level BLAS Level 2 and Level 3 routines are built.
--
-- Additionally the module exposes the common vector norms
-- (\(\lVert \cdot \rVert_1\) and \(\lVert \cdot \rVert_2\)) and
-- scalar–vector scaling, which together form the complete Level 1
-- BLAS surface.
--
-- == Complexity
--
-- All operations in this module are \(O(n)\) in the vector length.
module Numeric.LinearAlgebra.Massiv.BLAS.Level1
  ( -- * Dot product (Algorithm 1.1.1, GVL4 p. 4)
    dot
  , dotP
    -- * Scalar–vector operations (Algorithm 1.1.2, GVL4 p. 4)
  , scal
  , axpy
    -- * Vector norms (GVL4 Section 1.1, pp. 4–8)
  , nrm2
  , asum
  ) where

import qualified Data.Massiv.Array as M
import Data.Massiv.Array (unwrapByteArray, unwrapByteArrayOffset)
import GHC.TypeNats (KnownNat)

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.Internal.Kernel (rawDot)

-- | Inner (dot) product of two vectors.
--
-- __GVL4 Reference:__ Algorithm 1.1.1, p. 4.
--
-- Given \(x, y \in \mathbb{R}^{n}\), computes the scalar
--
-- \[
--   c \;=\; x^{T} y \;=\; \sum_{i=1}^{n} x_i \, y_i
-- \]
--
-- ==== Type-safety guarantees
--
-- Both vectors carry the same compile-time dimension @n@, so a
-- length mismatch is a type error.
--
-- ==== Complexity
--
-- \(O(n)\) — exactly \(n\) multiplications and \(n\) additions
-- (or \(n - 1\) additions, depending on the fold seed).
dot :: (KnownNat n, M.Manifest r e, Num e)
    => Vector n r e -> Vector n r e -> e
dot (MkVector x) (MkVector y) =
  M.foldlS (+) 0 $ M.zipWith (*) x y
{-# NOINLINE [1] dot #-}

-- | Specialised raw-array dot product for P Double.
dotP :: forall n. KnownNat n => Vector n M.P Double -> Vector n M.P Double -> Double
dotP (MkVector x) (MkVector y) =
  rawDot (unwrapByteArray x) (unwrapByteArrayOffset x)
         (unwrapByteArray y) (unwrapByteArrayOffset y)
         (dimVal @n)
{-# NOINLINE dotP #-}

{-# RULES "dot/P/Double" forall (x :: Vector n M.P Double)
                                (y :: Vector n M.P Double).
    dot x y = dotP x y #-}

-- | Scale every element of a vector by a scalar.
--
-- __GVL4 Reference:__ Section 1.1, pp. 4–8 (scalar–vector operations).
--
-- Computes
--
-- \[
--   x \;\leftarrow\; \alpha \, x
-- \]
--
-- i.e., each component \(x_i\) is replaced by \(\alpha \, x_i\).
--
-- ==== Type-safety guarantees
--
-- The output vector retains the same compile-time dimension @n@ as the
-- input.
--
-- ==== Complexity
--
-- \(O(n)\) — one multiplication per element.
scal :: (KnownNat n, M.Manifest r e, Num e)
     => e -> Vector n r e -> Vector n r e
scal alpha (MkVector x) = MkVector $ M.compute $ M.map (* alpha) x

-- | Saxpy (Scalar Alpha X Plus Y) — the fundamental vector-update operation.
--
-- __GVL4 Reference:__ Algorithm 1.1.1 (Saxpy), p. 4.
--
-- Given a scalar \(\alpha\) and vectors \(x, y \in \mathbb{R}^{n}\),
-- computes
--
-- \[
--   y \;\leftarrow\; \alpha \, x + y
-- \]
--
-- The Saxpy kernel is the innermost building block of the BLAS hierarchy.
-- Every Gaxpy (Level 2) and matrix–matrix (Level 3) operation can be
-- expressed as a sequence of Saxpy calls (GVL4, Section 1.1, p. 4).
--
-- ==== Type-safety guarantees
--
-- Both input vectors and the result share the same compile-time
-- dimension @n@.
--
-- ==== Complexity
--
-- \(O(n)\) — one fused multiply-add per element.
axpy :: (KnownNat n, M.Manifest r e, Num e)
     => e -> Vector n r e -> Vector n r e -> Vector n r e
axpy alpha (MkVector x) (MkVector y) =
  MkVector $ M.compute $ M.zipWith (\xi yi -> alpha * xi + yi) x y

-- | Euclidean (2-) norm of a vector.
--
-- __GVL4 Reference:__ Section 1.1, pp. 4–8 (vector norms).
--
-- Computes
--
-- \[
--   \lVert x \rVert_2 \;=\; \sqrt{\sum_{i=1}^{n} x_i^{2}}
-- \]
--
-- ==== Type-safety guarantees
--
-- The input vector carries its length @n@ at the type level; the result
-- is a scalar of the same element type.
--
-- ==== Complexity
--
-- \(O(n)\) — one multiply-accumulate per element, plus a single square
-- root.
--
-- /Note:/ This implementation does not perform the scaling trick
-- described in GVL4 (p. 5) to avoid overflow\/underflow for
-- extreme element magnitudes.  For production use on
-- floating-point data with very large or very small entries,
-- consider a two-pass scaled variant.
nrm2 :: (KnownNat n, M.Manifest r e, Floating e)
     => Vector n r e -> e
nrm2 (MkVector x) = sqrt $ M.foldlS (\acc xi -> acc + xi * xi) 0 x

-- | Sum of absolute values — the 1-norm (Manhattan norm) of a vector.
--
-- __GVL4 Reference:__ Section 1.1, pp. 4–8 (vector norms).
--
-- Computes
--
-- \[
--   \lVert x \rVert_1 \;=\; \sum_{i=1}^{n} \lvert x_i \rvert
-- \]
--
-- ==== Type-safety guarantees
--
-- The input vector carries its length @n@ at the type level; the result
-- is a scalar of the same element type.
--
-- ==== Complexity
--
-- \(O(n)\) — one absolute value and one addition per element.
asum :: (KnownNat n, M.Manifest r e, Num e, Ord e)
     => Vector n r e -> e
asum (MkVector x) = M.foldlS (\acc xi -> acc + abs xi) 0 x
