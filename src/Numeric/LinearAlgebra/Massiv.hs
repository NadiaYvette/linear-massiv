-- |
-- Module      : Numeric.LinearAlgebra.Massiv
-- Copyright   : (c) Nadia Chambers 2026
-- License     : BSD-3-Clause
-- Maintainer  : nadia.chambers@iohk.io
-- Stability   : experimental
--
-- @linear-massiv@: Type-safe numerical linear algebra backed by
-- <https://hackage.haskell.org/package/massiv massiv> arrays.
--
-- This library provides native Haskell implementations of algorithms from:
--
-- * Golub, G. H., & Van Loan, C. F. (2013). /Matrix Computations/ (4th ed.).
--   Johns Hopkins University Press. ISBN 978-1-4214-0794-4.
--
-- referred to throughout as __GVL4__.
--
-- = Design Principles
--
-- 1. __Type-level dimensional safety__: Matrix dimensions are tracked in the
--    type system via GHC @DataKinds@ and @KnownNat@ constraints. Dimensionally
--    incorrect operations (e.g., multiplying an \(m \times k\) matrix by an
--    \(n \times p\) matrix where \(k \neq n\)) are rejected at compile time.
--
-- 2. __Representation polymorphism__: All operations are parametric over
--    massiv's array representation @r@ (e.g., @P@ for primitive, @U@ for
--    unboxed, @S@ for storable, @B@ for boxed), constrained by
--    @'Data.Massiv.Array.Manifest' r e@. Users choose the representation at
--    the call site.
--
-- 3. __Parallelism via massiv__: Operations that construct arrays via
--    @makeArray@ inherit massiv's computation strategies. Use 'matMulComp'
--    with @Par@ or @ParN n@ for parallel matrix multiplication.
--
-- 4. __No FFI__: All algorithms are pure Haskell, enabling portability and
--    auditability. Benchmarks compare performance across massiv representations
--    and parallelism strategies.
--
-- = Quick Start
--
-- @
-- import Numeric.LinearAlgebra.Massiv
-- import qualified Data.Massiv.Array as M
--
-- -- Create a 3x3 matrix (Primitive representation, Double elements)
-- let a = 'makeMatrix' \@3 \@3 \@M.P $ \\i j ->
--           fromIntegral (i * 3 + j + 1) :: Double
--
-- -- QR factorization
-- let (q, r) = 'qr' a
--
-- -- Solve Ax = b via LU
-- let b = 'makeVector' \@3 \@M.P $ \\i -> fromIntegral (i + 1) :: Double
-- let x = 'luSolve' a b
-- @
--
-- = Module Organisation
--
-- == Core types and construction
--
-- * "Numeric.LinearAlgebra.Massiv.Types" — 'Matrix', 'Vector' newtypes with
--   phantom dimension parameters
-- * "Numeric.LinearAlgebra.Massiv.Internal" — Unsafe constructors, dimension
--   reification, array creation helpers
--
-- == BLAS-like operations (GVL4 Ch. 1)
--
-- * "Numeric.LinearAlgebra.Massiv.BLAS.Level1" — Vector–vector: 'dot', 'axpy', 'scal', 'nrm2'
-- * "Numeric.LinearAlgebra.Massiv.BLAS.Level2" — Matrix–vector: 'gemv', 'matvec', 'ger'
-- * "Numeric.LinearAlgebra.Massiv.BLAS.Level3" — Matrix–matrix: 'gemm', 'matMul', 'transpose'
--
-- == Direct solvers (GVL4 Chs. 3–4)
--
-- * "Numeric.LinearAlgebra.Massiv.Solve.Triangular" — Forward\/back substitution
-- * "Numeric.LinearAlgebra.Massiv.Solve.LU" — LU with partial pivoting
-- * "Numeric.LinearAlgebra.Massiv.Solve.Cholesky" — Cholesky for SPD matrices
-- * "Numeric.LinearAlgebra.Massiv.Solve.Banded" — Band LU, band Cholesky, tridiagonal solver
--
-- == Orthogonal factorizations (GVL4 Ch. 5)
--
-- * "Numeric.LinearAlgebra.Massiv.Orthogonal.Householder" — Householder reflections
-- * "Numeric.LinearAlgebra.Massiv.Orthogonal.Givens" — Givens rotations
-- * "Numeric.LinearAlgebra.Massiv.Orthogonal.QR" — QR factorization (Householder and Givens)
-- * "Numeric.LinearAlgebra.Massiv.Orthogonal.LeastSquares" — Least squares via QR
--
-- == Eigenvalue problems and SVD (GVL4 Chs. 7–8)
--
-- * "Numeric.LinearAlgebra.Massiv.Eigen.Power" — Power, inverse, Rayleigh quotient iteration
-- * "Numeric.LinearAlgebra.Massiv.Eigen.Hessenberg" — Hessenberg reduction
-- * "Numeric.LinearAlgebra.Massiv.Eigen.Schur" — Schur decomposition (QR algorithm)
-- * "Numeric.LinearAlgebra.Massiv.Eigen.Symmetric" — Symmetric eigenvalue (tridiagonal QR, Jacobi)
-- * "Numeric.LinearAlgebra.Massiv.Eigen.SVD" — Singular value decomposition
--
-- == Norms and condition numbers (GVL4 Ch. 2)
--
-- * "Numeric.LinearAlgebra.Massiv.Norms" — Frobenius, 1-, \(\infty\)-, and 2-norms
--
-- == Integration with the @linear@ library
--
-- * "Numeric.LinearAlgebra.Massiv.Linear" — Conversions to\/from @linear@'s @V@, @V2@, @V3@, @V4@
module Numeric.LinearAlgebra.Massiv
  ( -- * Core types
    module Numeric.LinearAlgebra.Massiv.Types
    -- * Construction helpers
  , module Numeric.LinearAlgebra.Massiv.Internal
    -- * BLAS operations
  , module Numeric.LinearAlgebra.Massiv.BLAS.Level1
  , module Numeric.LinearAlgebra.Massiv.BLAS.Level2
  , module Numeric.LinearAlgebra.Massiv.BLAS.Level3
    -- * Direct solvers
  , module Numeric.LinearAlgebra.Massiv.Solve.Triangular
  , module Numeric.LinearAlgebra.Massiv.Solve.LU
  , module Numeric.LinearAlgebra.Massiv.Solve.Cholesky
  , module Numeric.LinearAlgebra.Massiv.Solve.Banded
    -- * Orthogonal factorizations
  , module Numeric.LinearAlgebra.Massiv.Orthogonal.Householder
  , module Numeric.LinearAlgebra.Massiv.Orthogonal.Givens
  , module Numeric.LinearAlgebra.Massiv.Orthogonal.QR
  , module Numeric.LinearAlgebra.Massiv.Orthogonal.LeastSquares
    -- * Eigenvalue problems
  , module Numeric.LinearAlgebra.Massiv.Eigen.Power
  , module Numeric.LinearAlgebra.Massiv.Eigen.Hessenberg
  , module Numeric.LinearAlgebra.Massiv.Eigen.Schur
  , module Numeric.LinearAlgebra.Massiv.Eigen.Symmetric
  , module Numeric.LinearAlgebra.Massiv.Eigen.SVD
    -- * Norms
  , module Numeric.LinearAlgebra.Massiv.Norms
    -- * Linear integration
  , module Numeric.LinearAlgebra.Massiv.Linear
  ) where

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.BLAS.Level1
import Numeric.LinearAlgebra.Massiv.BLAS.Level2
import Numeric.LinearAlgebra.Massiv.BLAS.Level3
import Numeric.LinearAlgebra.Massiv.Solve.Triangular
import Numeric.LinearAlgebra.Massiv.Solve.LU
import Numeric.LinearAlgebra.Massiv.Solve.Cholesky
import Numeric.LinearAlgebra.Massiv.Solve.Banded
import Numeric.LinearAlgebra.Massiv.Orthogonal.Householder
import Numeric.LinearAlgebra.Massiv.Orthogonal.Givens
import Numeric.LinearAlgebra.Massiv.Orthogonal.QR
import Numeric.LinearAlgebra.Massiv.Orthogonal.LeastSquares
import Numeric.LinearAlgebra.Massiv.Eigen.Power
import Numeric.LinearAlgebra.Massiv.Eigen.Hessenberg
import Numeric.LinearAlgebra.Massiv.Eigen.Schur
import Numeric.LinearAlgebra.Massiv.Eigen.Symmetric
import Numeric.LinearAlgebra.Massiv.Eigen.SVD
import Numeric.LinearAlgebra.Massiv.Norms
import Numeric.LinearAlgebra.Massiv.Linear
