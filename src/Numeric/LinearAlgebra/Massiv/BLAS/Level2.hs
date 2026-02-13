{-# LANGUAGE AllowAmbiguousTypes #-}

-- |
-- Module      : Numeric.LinearAlgebra.Massiv.BLAS.Level2
-- Copyright   : (c) Nadia Chambers 2026
-- License     : BSD-3-Clause
-- Maintainer  : nadia.chambers@iohk.io
-- Stability   : experimental
--
-- = BLAS Level 2: Matrix–Vector Operations
--
-- This module provides type-safe, dimension-indexed wrappers around the
-- standard BLAS Level 2 kernels.  These are /matrix–vector/ operations
-- whose arithmetic cost is \(O(m \, n)\) for an \(m \times n\) matrix,
-- one level above the \(O(n)\) vector–vector operations of Level 1.
--
-- The central operation is the /Gaxpy/ — Generalized Saxpy — which
-- computes \(y \leftarrow \alpha A x + \beta y\).  It can be viewed
-- as a sequence of Saxpy updates (one per row or one per column),
-- giving rise to two natural loop orderings:
--
-- * Golub, G. H. & Van Loan, C. F. (2013). /Matrix Computations/,
--   4th edition (GVL4). Johns Hopkins University Press.
--   __Chapter 1, Sections 1.1.3–1.1.4__, pp. 8–12.
--
-- Specifically:
--
-- * __Algorithm 1.1.3__ (Row-Oriented Gaxpy, p. 8) — The @i@-th
--   component of the result is a dot product:
--   \(y_i \leftarrow a_i^{T} x + y_i\), where \(a_i^{T}\) is the
--   @i@-th row of \(A\).
--
-- * __Algorithm 1.1.4__ (Column-Oriented Gaxpy, p. 9) — The result
--   vector is updated one column at a time via Saxpy:
--   \(y \leftarrow A(:,\!j) \, x_j + y\), for \(j = 1, \ldots, n\).
--
-- The module also provides the rank-1 outer-product update
-- \(A \leftarrow A + \alpha x y^{T}\) (BLAS @GER@), which is the
-- matrix analogue of the Saxpy at Level 1 and plays a key role in LU
-- factorisation (GVL4, Section 3.2, p. 112).
--
-- == Complexity
--
-- All operations in this module are \(O(m \, n)\) for an
-- \(m \times n\) matrix.
module Numeric.LinearAlgebra.Massiv.BLAS.Level2
  ( -- * Matrix–vector multiply — Gaxpy (Algorithms 1.1.3–1.1.4, GVL4 pp. 8–9)
    gemv
  , matvec
    -- * Rank-1 update (GVL4 Section 1.1.4, p. 10)
  , ger
  ) where

import qualified Data.Massiv.Array as M
import Data.Massiv.Array (Ix2(..), Sz(..), Comp(..))
import GHC.TypeNats (KnownNat)

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal

-- | General matrix–vector multiply (BLAS @GEMV@).
--
-- __GVL4 Reference:__ Algorithms 1.1.3 (Row Gaxpy, p. 8) and 1.1.4
-- (Column Gaxpy, p. 9).
--
-- Given \(A \in \mathbb{R}^{m \times n}\),
-- \(x \in \mathbb{R}^{n}\), \(y \in \mathbb{R}^{m}\), and scalars
-- \(\alpha, \beta\), computes
--
-- \[
--   y \;\leftarrow\; \alpha \, A \, x \;+\; \beta \, y
-- \]
--
-- The implementation uses a /row-oriented/ traversal (Algorithm 1.1.3):
-- for each row \(i\) the dot product \(a_i^{T} x\) is formed, scaled by
-- \(\alpha\), and added to \(\beta \, y_i\).
--
-- ==== Type-safety guarantees
--
-- * \(A\) is @m x n@, \(x\) is @n@, \(y\) and the result are @m@.
-- * The shared inner dimension @n@ is enforced at compile time, so a
--   dimension mismatch is a type error.
--
-- ==== Complexity
--
-- \(O(m \, n)\) — two floating-point operations per matrix entry
-- (one multiply, one add in the inner product), plus \(O(m)\) work
-- for the \(\alpha / \beta\) scaling.
gemv :: forall m n r e. (KnownNat m, KnownNat n, M.Manifest r e, Num e)
     => e -> Matrix m n r e -> Vector n r e -> e -> Vector m r e -> Vector m r e
gemv alpha mat x beta y =
  let r = dimVal @m
      c = dimVal @n
  in makeVector @m @r $ \i ->
    let axi = foldl (\acc j -> acc + (mat ! (i, j)) * (x !. j)) 0 [0..c-1]
    in alpha * axi + beta * (y !. i)

-- | Simple matrix–vector multiply (specialisation of 'gemv').
--
-- __GVL4 Reference:__ Algorithm 1.1.3 (Row Gaxpy, p. 8), with
-- \(\alpha = 1\) and \(\beta = 0\).
--
-- Given \(A \in \mathbb{R}^{m \times n}\) and
-- \(x \in \mathbb{R}^{n}\), computes
--
-- \[
--   y \;=\; A \, x
-- \]
--
-- This is a convenience wrapper equivalent to @'gemv' 1 a x 0 zero@
-- but avoids allocating or requiring an initial @y@ vector.
--
-- ==== Type-safety guarantees
--
-- * \(A\) is @m x n@, \(x\) is @n@, the result is @m@.
-- * The inner dimension @n@ is checked at compile time.
--
-- ==== Complexity
--
-- \(O(m \, n)\).
matvec :: forall m n r e. (KnownNat m, KnownNat n, M.Manifest r e, Num e)
       => Matrix m n r e -> Vector n r e -> Vector m r e
matvec mat x =
  let c = dimVal @n
  in makeVector @m @r $ \i ->
    foldl (\acc j -> acc + (mat ! (i, j)) * (x !. j)) 0 [0..c-1]

-- | Rank-1 update — outer product (BLAS @GER@).
--
-- __GVL4 Reference:__ Section 1.1.4, p. 10.  The rank-1 update is
-- the matrix-level analogue of the Saxpy and appears as the inner
-- kernel in outer-product formulations of LU factorisation
-- (GVL4 Section 3.2, Algorithm 3.2.1, p. 112).
--
-- Given \(x \in \mathbb{R}^{m}\), \(y \in \mathbb{R}^{n}\),
-- \(A \in \mathbb{R}^{m \times n}\), and a scalar \(\alpha\),
-- computes
--
-- \[
--   A \;\leftarrow\; A \;+\; \alpha \, x \, y^{T}
-- \]
--
-- Equivalently, each entry is updated as
-- \(a_{ij} \leftarrow a_{ij} + \alpha \, x_i \, y_j\).
--
-- ==== Type-safety guarantees
--
-- * \(x\) is @m@, \(y\) is @n@, \(A\) and the result are @m x n@.
-- * All dimension constraints are enforced at compile time.
--
-- ==== Complexity
--
-- \(O(m \, n)\) — one fused multiply-add per matrix entry.
ger :: forall m n r e. (KnownNat m, KnownNat n, M.Manifest r e, Num e)
    => e -> Vector m r e -> Vector n r e -> Matrix m n r e -> Matrix m n r e
ger alpha x y mat =
  makeMatrix @m @n @r $ \i j ->
    (mat ! (i, j)) + alpha * (x !. i) * (y !. j)
