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
-- = Derived Work Attribution
--
-- This code was co-authored by Claude Opus (Anthropic) and should be
-- considered a derived work of the various algorithmic examples and
-- reference implementations drawn upon during development, including but
-- not limited to:
--
-- * __LAPACK__ (Linear Algebra PACKage) — Anderson, E. et al. (1999).
--   /LAPACK Users' Guide/, 3rd ed., SIAM. The LAPACK testing methodology,
--   algorithm structures, and numerical stability techniques informed much
--   of the implementation.
--
-- * __OpenBLAS__ — Xianyi, Z., Qian, W., and Yunquan, Z. (2011--).
--   The tiled GEMM micro-kernel architecture, cache-blocking strategies,
--   and SIMD vectorisation patterns were inspired by OpenBLAS.
--
-- * __GVL4__ — The primary algorithmic reference, as noted above.
--
-- * __Higham__ — Higham, N. J. (2002). /Accuracy and Stability of Numerical
--   Algorithms/, 2nd ed., SIAM. Error analysis and numerical stability
--   frameworks.
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
-- = Internal Architecture: Two-Layer Design
--
-- @linear-massiv@ uses a two-layer architecture that separates the type-safe
-- public API from the performance-critical internal representation:
--
-- * __Public layer__: 'Matrix' and 'Vector' are @newtype@ wrappers around
--   massiv's @Array r Ix2 e@, providing compile-time dimension checking via
--   phantom @Nat@ parameters and representation polymorphism via @r@.
--
-- * __Internal layer__: Performance-critical inner loops (GEMM, QR,
--   tridiagonalisation, SVD, etc.) unwrap the massiv array to a raw
--   @ByteArray#@ \/ @MutableByteArray#@ and operate directly via GHC primops,
--   including @DoubleX4#@ AVX2 SIMD instructions compiled through the LLVM 17
--   backend.  Functions receive @(ByteArray, offset, stride)@ triples, enabling
--   zero-copy submatrix views for panel factorisations.
--
-- This separation is essential for performance.  Benchmarks (Round 3 of the
-- accompanying report) showed that massiv's per-element @M.readM@\/@M.write_@
-- abstraction layer imposed a 240–330× penalty on BLAS operations relative to
-- direct primop access, even though the underlying memory layout is identical.
-- The raw primop layer eliminates this overhead while the @newtype@ wrapper
-- preserves type safety at the API boundary.
--
-- == Why not @vector-sized@ or @linear@'s @V@?
--
-- The @<https://hackage.haskell.org/package/vector-sized vector-sized>@ package
-- provides an @Unbox (Vector n a)@ instance that stores
-- @Vector m (Vector n Double)@ as a contiguous flat @ByteArray@ of @m × n@
-- doubles.  While the __memory layout__ is correct, contiguous memory alone is
-- insufficient for high-performance numerical kernels:
--
-- * __Per-element typeclass dispatch__: Every access goes through
--   @basicUnsafeRead@ \/ @basicUnsafeWrite@ of the @Unbox@ data family.
--   Reading element @(i, j)@ requires indexing the outer vector to obtain a
--   @Vector n Double@ (constructing an intermediate slice), then indexing that.
--   @linear-massiv@ computes @off + i * stride + j@ and issues a single
--   @readDoubleArray#@ primop.
--
-- * __No SIMD access__: The 4×8 GEMM micro-kernel loads 4 consecutive doubles
--   via @indexDoubleArrayAsDoubleX4# ba# (off# +# i#)@—a direct 256-bit AVX2
--   load from a computed byte offset.  The @Unbox@ typeclass does not expose
--   the underlying @ByteArray#@, and GHC cannot optimise through the data
--   family indirection to produce equivalent code.
--
-- * __No mutable primop access__: In-place factorisations (LU, QR,
--   tridiagonalisation, bidiagonalisation) require @writeDoubleArray#@ on
--   @MutableByteArray#@ with computed offsets.  The @MVector@ abstraction
--   interposes allocation and function calls that prevent the tight unboxed
--   loops needed for competitive performance.
--
-- * __No zero-copy submatrix views__: Panel factorisations pass
--   @(ByteArray, offset, stride)@ triples to kernels, enabling zero-copy views
--   into submatrices.  @Vector n a@ does not naturally express "this row starts
--   at byte offset X in a larger backing array."
--
-- The @<https://hackage.haskell.org/package/linear linear>@ library's @V n a@
-- uses @Vector@ from the @vector@ package internally and is designed for small
-- fixed-size vectors (V2–V4) where GHC can fully unbox everything.  At
-- @n = 100–500@, @V n (V m Double)@ would be a vector-of-vectors with per-row
-- indirection—catastrophic for cache locality and SIMD vectorisation.
-- @linear-massiv@ provides conversion functions ('fromLinearV', 'fromV2', etc.)
-- in "Numeric.LinearAlgebra.Massiv.Linear" for interoperability.
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
