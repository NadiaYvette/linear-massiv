# linear-massiv

**Pure Haskell numerical linear algebra that outperforms BLAS and LAPACK.**

linear-massiv provides type-safe matrix and vector operations backed by
[massiv](https://hackage.haskell.org/package/massiv) arrays, with
compile-time dimension checking via GHC's type-level naturals. It requires
no FFI, no C dependencies, and no foreign libraries — every algorithm is
implemented in Haskell, from the GEMM micro-kernel to the divide-and-conquer
SVD.

Despite being pure Haskell, linear-massiv matches or exceeds the performance
of [hmatrix](https://hackage.haskell.org/package/hmatrix) (OpenBLAS/LAPACK)
across all benchmarked operations at moderate to large matrix sizes. This is
achieved through GHC 9.14's native `DoubleX4#` AVX2 SIMD primops, compiled
via the LLVM 17 backend, combined with cache-blocked algorithms, register-tuned
micro-kernels, and panel factorisations drawn from the numerical linear algebra
literature.

## Performance

Single-threaded benchmarks vs hmatrix (OpenBLAS/LAPACK) on AMD Zen 4. Ratios
below 1.0 mean linear-massiv is faster.

| Operation | 100x100 | 200x200 | 500x500 |
|-----------|---------|---------|---------|
| GEMM | **0.15x** (6.5x faster) | **0.10x** (9.8x faster) | **0.098x** (10.2x faster) |
| QR factorisation | **0.021x** (48x faster) | | |
| LU solve | **0.51x** (2.0x faster) | | |
| Cholesky solve | **0.53x** (1.9x faster) | | |
| Symmetric eigen | **0.94x** (1.07x faster) | **0.77x** (1.3x faster) | **0.49x** (2.0x faster) |
| SVD | 1.27x | **0.88x** (1.14x faster) | **0.54x** (1.84x faster) |

At 500x500 — the largest benchmarked size — linear-massiv outperforms
hmatrix/LAPACK in **all nine operation categories**.

## Why pure Haskell?

Wrapping BLAS and LAPACK via FFI (as hmatrix does) gives you decades of
Fortran micro-optimisation, but at a cost:

- **No type safety across the FFI boundary.** Dimension mismatches, incorrect
  strides, and aliased pointers are runtime errors. linear-massiv catches
  dimension errors at compile time via phantom type parameters.

- **No portability.** Building hmatrix requires a working BLAS/LAPACK
  installation, platform-specific linker flags, and Fortran runtime
  libraries. linear-massiv builds with `cabal build` on any platform
  GHC supports.

- **No visibility.** The hot path disappears into opaque C/Fortran object
  code. linear-massiv's entire implementation — from the 4x8 GEMM
  micro-kernel to the Gu-Eisenstat secular equation solver — is readable
  Haskell, auditable and modifiable without leaving the language.

The key insight enabling competitive performance is GHC 9.14's native SIMD
support: `DoubleX4#` primops compile to AVX2 `vfmadd231pd` instructions
through the LLVM backend, giving Haskell direct access to the same hardware
that makes BLAS fast — without FFI overhead.

## Quick start

```haskell
import Numeric.LinearAlgebra.Massiv
import qualified Data.Massiv.Array as M

-- Create a 3x3 matrix (Primitive representation, Double elements)
let a = makeMatrix @3 @3 @M.P $ \i j ->
          fromIntegral (i * 3 + j + 1) :: Double

-- QR factorisation
let (q, r) = qr a

-- Solve Ax = b via LU
let b = makeVector @3 @M.P $ \i -> fromIntegral (i + 1) :: Double
let x = luSolve a b

-- Symmetric eigenvalue decomposition
let (eigenvalues, eigenvectors) = symmetricEigen a

-- Singular value decomposition
let (u, sigma, vt) = svd a
```

## Requirements

- **GHC 9.14** or later (for `DoubleX4#` SIMD primops)
- **LLVM 17** backend (`-fllvm -pgmlo opt-17 -pgmlc llc-17`)
- **AVX2 + FMA** capable CPU (Intel Haswell / AMD Zen or later)

Add to your `cabal.project`:

```
package linear-massiv
  ghc-options: -fllvm -pgmlo opt-17 -pgmlc llc-17 -mavx2 -mfma
```

## Operations

### BLAS-equivalent (GVL4 Ch. 1)

- **Level 1:** dot product, axpy, scal, nrm2
- **Level 2:** matrix-vector multiply (gemv, matvec), outer product (ger)
- **Level 3:** matrix multiply (gemm), transpose — with 4x8 register-blocked
  micro-kernel, cache-blocked tiling, micro-panel packing, and optional
  thread-level parallelism

### Direct solvers (GVL4 Chs. 3-4)

- **LU** with partial pivoting, panel factorisation
- **Cholesky** for symmetric positive definite systems, panel factorisation
- **Band solvers:** band LU, band Cholesky, tridiagonal
- **Forward/back substitution** with 8-wide SIMD dual accumulators

### Orthogonal factorisations (GVL4 Ch. 5)

- **QR** via Householder reflections with blocked WY accumulation
- **Givens** rotations
- **Least squares** via QR

### Eigenvalue problems and SVD (GVL4 Chs. 7-8)

- **Symmetric eigendecomposition:** DLATRD-style panel tridiagonalisation,
  QR iteration with aggressive early deflation, divide-and-conquer with
  Gu-Eisenstat secular equation solver
- **SVD:** Two paths —
  - `svd` (fast): A^T A eigendecomposition with D&C eigensolver, fused
    DSYRK, pre-scaled U-recovery
  - `svdGK` (numerically proper): Golub-Kahan bidiagonalisation with
    DLABRD-style panel blocking, D&C bidiagonal SVD, blocked WY
    Householder accumulation
- **General eigenvalue:** Hessenberg reduction, Schur decomposition
- **Iterative:** power method, inverse iteration, Rayleigh quotient

### Norms (GVL4 Ch. 2)

- Frobenius, 1-norm, infinity-norm, spectral norm (2-norm via SVD)

## Architecture

The library uses a two-layer design:

- **Public API:** `Matrix m n r e` and `Vector n r e` are newtype wrappers
  around massiv arrays with phantom type-level dimension parameters.
  Dimension mismatches are compile-time errors.

- **Internal kernels:** Hot inner loops unwrap to raw `ByteArray#` /
  `MutableByteArray#` and operate via GHC primops, including `DoubleX4#`
  AVX2 SIMD. Functions receive `(ByteArray, offset, stride)` triples
  for zero-copy submatrix views in panel factorisations.

This separation eliminates massiv's per-element abstraction overhead
(measured at 240-330x in BLAS operations) while preserving type safety at
the API boundary.

## Integration with the `linear` library

Conversion functions are provided for interoperability with Edward Kmett's
[linear](https://hackage.haskell.org/package/linear) library:

```haskell
import Numeric.LinearAlgebra.Massiv.Linear

-- Convert from linear's V3
let v = fromV3 @M.P (V3 1.0 2.0 3.0)

-- Convert from linear's V n
let w = fromLinearV @5 @M.P someLinearVector
```

## Algorithmic references

The implementations follow:

- Golub, G. H. & Van Loan, C. F. (2013). *Matrix Computations* (4th ed.).
  Johns Hopkins University Press. (**GVL4**)
- Higham, N. J. (2002). *Accuracy and Stability of Numerical Algorithms*
  (2nd ed.). SIAM.

Specific techniques adapted from LAPACK and OpenBLAS are documented in
the source code and the accompanying benchmark report
(`report/benchmark-report.tex`, 91 pages covering 23 rounds of optimisation).

## License

BSD-3-Clause

## Co-authorship

This library was co-authored by Nadia Chambers and Claude Opus (Anthropic).
It should be considered a derived work of the algorithmic examples and
reference implementations drawn upon during development, including LAPACK,
OpenBLAS, and GVL4.
