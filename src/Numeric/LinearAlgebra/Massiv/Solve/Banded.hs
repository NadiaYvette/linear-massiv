{-# LANGUAGE AllowAmbiguousTypes #-}

-- |
-- Module      : Numeric.LinearAlgebra.Massiv.Solve.Banded
-- Copyright   : (c) Nadia Chambers 2026
-- License     : BSD-3-Clause
-- Maintainer  : nadia.chambers@iohk.io
-- Stability   : experimental
--
-- = Banded and Tridiagonal System Solvers
--
-- Specialised LU and Cholesky factorizations for banded matrices, plus a
-- dedicated tridiagonal solver, following Golub & Van Loan,
-- /Matrix Computations/, 4th edition (GVL4), Section 4.3, pp. 174--182.
--
-- A matrix \(A \in \mathbb{R}^{n \times n}\) has /lower bandwidth/ \(p\)
-- and /upper bandwidth/ \(q\) when \(a_{ij} = 0\) for \(i > j + p\) or
-- \(j > i + q\). Exploiting this structure reduces the factorization cost
-- from \(O(n^3)\) to \(O(npq)\) (GVL4 p. 176), and the triangular-solve
-- cost from \(O(n^2)\) to \(O(np)\) or \(O(nq)\). The important special
-- case of a tridiagonal system (\(p = q = 1\)) is solvable in \(O(n)\)
-- flops.
--
-- +-------------------+-----------------------------------+----------------------------------+
-- | Function          | Algorithm                         | Reference                        |
-- +===================+===================================+==================================+
-- | 'bandLU'          | Band Gaussian elimination          | GVL4 Algorithm 4.3.1, p. 175    |
-- +-------------------+-----------------------------------+----------------------------------+
-- | 'bandForwardSub'  | Band forward substitution          | GVL4 Algorithm 4.3.2, p. 176    |
-- +-------------------+-----------------------------------+----------------------------------+
-- | 'bandBackSub'     | Band back substitution             | GVL4 Algorithm 4.3.3, p. 176    |
-- +-------------------+-----------------------------------+----------------------------------+
-- | 'bandCholesky'    | Band Cholesky factorization        | GVL4 Algorithm 4.3.5, p. 178    |
-- +-------------------+-----------------------------------+----------------------------------+
-- | 'tridiagSolve'    | SPD tridiagonal solver             | GVL4 Algorithm 4.3.6, p. 179    |
-- +-------------------+-----------------------------------+----------------------------------+
--
-- == Complexity
--
-- * Band LU: \(O(npq)\) flops (GVL4 p. 176).
-- * Band triangular solve: \(O(np)\) or \(O(nq)\) flops (GVL4 p. 176).
-- * Band Cholesky: \(O(np^2)\) flops (GVL4 p. 178).
-- * Tridiagonal solver: \(O(n)\) flops (GVL4 p. 179).
--
-- == Type Safety
--
-- Matrix dimensions are tracked at the type level via 'KnownNat'. The
-- bandwidths \(p\) and \(q\) are passed as run-time 'Int' values because
-- they are often data-dependent. The matrix is stored in standard dense
-- format; only the band is accessed.
module Numeric.LinearAlgebra.Massiv.Solve.Banded
  ( -- * Band LU factorization
    bandLU
    -- * Band triangular solve
  , bandForwardSub
  , bandBackSub
    -- * Band Cholesky (\(A = GG^T\))
  , bandCholesky
    -- * Tridiagonal solver (\(Ax = b\), \(p = q = 1\))
  , tridiagSolve
  ) where

import qualified Data.Massiv.Array as M
import Data.Massiv.Array (Ix1, Ix2(..), Sz(..))
import GHC.TypeNats (KnownNat)

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal

-- | Band Gaussian elimination without pivoting (GVL4 Algorithm 4.3.1,
-- p. 175).
--
-- Given an \(n \times n\) matrix \(A\) with lower bandwidth \(p\) and upper
-- bandwidth \(q\), computes an in-place \(LU\) factorization where \(L\)
-- has lower bandwidth \(p\) and \(U\) has upper bandwidth \(q\).
--
-- The matrix is stored in standard dense format; only the entries within
-- the band are accessed or modified.
--
-- __Precondition.__ All leading principal submatrices of \(A\) must be
-- nonsingular (analogous to the dense case, GVL4 Theorem 3.2.1, p. 116).
--
-- ==== Type-safety guarantees
--
-- 'KnownNat' \(n\) statically ensures the matrix is square. The bandwidth
-- parameters \(p\) and \(q\) are run-time values.
--
-- ==== Complexity
--
-- \(O(npq)\) flops (GVL4 p. 176).
--
-- ==== Reference
--
-- Golub & Van Loan, /Matrix Computations/, 4th ed., Algorithm 4.3.1
-- (Band Gaussian Elimination), p. 175.
bandLU :: forall n r e. (KnownNat n, M.Manifest r e, Fractional e)
       => Int  -- ^ Lower bandwidth \(p\)
       -> Int  -- ^ Upper bandwidth \(q\)
       -> Matrix n n r e -> Matrix n n r e
bandLU p q (MkMatrix a) =
  let nn = dimVal @n
  in MkMatrix $ snd $ M.withMArrayST a $ \ma ->
    mapM_ (\k -> do
      akk <- M.readM ma (k :. k)
      let imax = min (k + p) (nn - 1)
      -- Compute multipliers
      mapM_ (\i -> do
        aik <- M.readM ma (i :. k)
        M.write_ ma (i :. k) (aik / akk)
        ) [k+1..imax]
      -- Update
      let jmax = min (k + q) (nn - 1)
      mapM_ (\j ->
        mapM_ (\i -> do
          aij <- M.readM ma (i :. j)
          aik <- M.readM ma (i :. k)
          akj <- M.readM ma (k :. j)
          M.write_ ma (i :. j) (aij - aik * akj)
          ) [k+1..imax]
        ) [k+1..jmax]
      ) [0..nn-2]

-- | Band forward substitution (GVL4 Algorithm 4.3.2, p. 176).
--
-- Solves \(Lx = b\) where \(L\) is /unit/ lower triangular with lower
-- bandwidth \(p\). Only the \(p\) subdiagonals are accessed; the unit
-- diagonal is implicit, so no division is performed and the constraint
-- relaxes to 'Num'.
--
-- ==== Type-safety guarantees
--
-- 'KnownNat' \(n\) ensures the dimensions of \(L\) and \(b\) agree at
-- compile time.
--
-- ==== Complexity
--
-- \(O(np)\) flops (GVL4 p. 176).
--
-- ==== Reference
--
-- Golub & Van Loan, /Matrix Computations/, 4th ed., Algorithm 4.3.2
-- (Band Forward Substitution), p. 176.
bandForwardSub :: forall n r e. (KnownNat n, M.Manifest r e, Num e)
               => Int  -- ^ Lower bandwidth \(p\)
               -> Matrix n n r e -> Vector n r e -> Vector n r e
bandForwardSub p l b = createVector @n $ \mx -> do
  let nn = dimVal @n
  mapM_ (\i -> M.write_ mx i (b !. i)) [0..nn-1]
  mapM_ (\j -> do
    xj <- M.readM mx j
    let imax = min (j + p) (nn - 1)
    mapM_ (\i -> do
      xi <- M.readM mx i
      M.write_ mx i (xi - (l ! (i, j)) * xj)
      ) [j+1..imax]
    ) [0..nn-1]

-- | Band back substitution (GVL4 Algorithm 4.3.3, p. 176).
--
-- Solves \(Ux = b\) where \(U\) is upper triangular with upper bandwidth
-- \(q\). Only the diagonal and the \(q\) superdiagonals are accessed.
--
-- ==== Type-safety guarantees
--
-- 'KnownNat' \(n\) ensures the dimensions of \(U\) and \(b\) agree at
-- compile time.
--
-- ==== Complexity
--
-- \(O(nq)\) flops (GVL4 p. 176).
--
-- ==== Reference
--
-- Golub & Van Loan, /Matrix Computations/, 4th ed., Algorithm 4.3.3
-- (Band Back Substitution), p. 176.
bandBackSub :: forall n r e. (KnownNat n, M.Manifest r e, Fractional e)
            => Int  -- ^ Upper bandwidth \(q\)
            -> Matrix n n r e -> Vector n r e -> Vector n r e
bandBackSub q u b = createVector @n $ \mx -> do
  let nn = dimVal @n
  mapM_ (\i -> M.write_ mx i (b !. i)) [0..nn-1]
  mapM_ (\j -> do
    xj <- M.readM mx j
    let ujj = u ! (j, j)
        xj' = xj / ujj
    M.write_ mx j xj'
    let imin = max 0 (j - q)
    mapM_ (\i -> do
      xi <- M.readM mx i
      M.write_ mx i (xi - (u ! (i, j)) * xj')
      ) [imin..j-1]
    ) [nn-1, nn-2..0]

-- | Band Cholesky factorization (GVL4 Algorithm 4.3.5, p. 178).
--
-- Given a symmetric positive definite \(n \times n\) matrix \(A\) with
-- lower bandwidth \(p\), computes the lower triangular factor \(G\) (also
-- with bandwidth \(p\)) such that
--
-- \[
--   A = G G^T
-- \]
--
-- Only the lower band of \(A\) (entries \(a_{ij}\) with
-- \(0 \le i - j \le p\)) is read.
--
-- __Precondition.__ \(A\) must be symmetric positive definite. If this
-- condition is violated, the algorithm may encounter a negative value under
-- a square root and produce NaN.
--
-- ==== Type-safety guarantees
--
-- 'KnownNat' \(n\) statically ensures the matrix is square. The bandwidth
-- \(p\) is a run-time value.
--
-- ==== Complexity
--
-- \(O(np^2)\) flops (GVL4 p. 178).
--
-- ==== Reference
--
-- Golub & Van Loan, /Matrix Computations/, 4th ed., Algorithm 4.3.5
-- (Band Cholesky), p. 178.
bandCholesky :: forall n r e. (KnownNat n, M.Manifest r e, Floating e, Ord e)
             => Int  -- ^ Bandwidth \(p\)
             -> Matrix n n r e -> Matrix n n r e
bandCholesky p (MkMatrix a) =
  let nn = dimVal @n
  in MkMatrix $ M.createArrayST_ (M.Sz2 nn nn) $ \mg -> do
    -- Initialize to zero
    mapM_ (\i -> mapM_ (\j -> M.write_ mg (i :. j) 0) [0..nn-1]) [0..nn-1]
    -- Copy lower band of A
    mapM_ (\j ->
      let imax = min (j + p) (nn - 1)
      in mapM_ (\i -> M.write_ mg (i :. j) (M.index' a (i :. j))) [j..imax]
      ) [0..nn-1]
    -- Band Cholesky
    mapM_ (\j -> do
      -- Subtract contributions
      let kmin = max 0 (j - p)
      mapM_ (\k -> do
        gjk <- M.readM mg (j :. k)
        let lam = min (j + p) (nn - 1)
        mapM_ (\i -> do
          gij <- M.readM mg (i :. j)
          gik <- M.readM mg (i :. k)
          M.write_ mg (i :. j) (gij - gik * gjk)
          ) [j..lam]
        ) [kmin..j-1]
      -- Scale
      gjj <- M.readM mg (j :. j)
      let sjj = sqrt gjj
          lam = min (j + p) (nn - 1)
      mapM_ (\i -> do
        gij <- M.readM mg (i :. j)
        M.write_ mg (i :. j) (gij / sjj)
        ) [j..lam]
      ) [0..nn-1]

-- | Symmetric positive definite tridiagonal system solver (GVL4
-- Algorithm 4.3.6, p. 179).
--
-- Solves \(Ax = b\) where \(A\) is symmetric, tridiagonal, and positive
-- definite. The matrix \(A\) is specified compactly by its diagonal
-- \(\alpha_{1:n}\) and its superdiagonal \(\beta_{1:n-1}\).
--
-- The algorithm computes the \(LDL^T\) factorization of \(A\) and folds it
-- together with forward and back substitution in a single \(O(n)\) pass:
--
-- 1. __Factor:__ Compute \(A = LDL^T\) where \(L\) is unit lower
--    bidiagonal and \(D\) is diagonal.
-- 2. __Forward substitution:__ Solve \(Lz = b\).
-- 3. __Diagonal solve:__ Solve \(Dy = z\).
-- 4. __Back substitution:__ Solve \(L^T x = y\).
--
-- __Precondition.__ \(A\) must be symmetric positive definite. This is not
-- checked at run time.
--
-- ==== Type-safety guarantees
--
-- 'KnownNat' \(n\) enforces that the diagonal, superdiagonal, right-hand
-- side, and solution vectors all have length \(n\) at compile time.
--
-- ==== Complexity
--
-- \(O(n)\) flops (GVL4 p. 179). This is optimal for tridiagonal systems.
--
-- ==== Reference
--
-- Golub & Van Loan, /Matrix Computations/, 4th ed., Algorithm 4.3.6
-- (SPD Tridiagonal System Solver), p. 179.
tridiagSolve :: forall n r e. (KnownNat n, M.Manifest r e, Fractional e)
             => Vector n r e     -- ^ Diagonal entries \(\alpha_{1:n}\)
             -> Vector n r e     -- ^ Superdiagonal entries \(\beta_{1:n-1}\) (length \(n\), only indices @0..n-2@ used)
             -> Vector n r e     -- ^ Right-hand side \(b\)
             -> Vector n r e     -- ^ Solution \(x\)
tridiagSolve diag supdiag b = createVector @n $ \mx -> do
  let nn = dimVal @n
  -- Working arrays for modified diagonal and superdiagonal
  alpha <- M.newMArray @r (M.Sz1 nn) (0 :: e)
  beta  <- M.newMArray @r (M.Sz1 nn) (0 :: e)
  -- Initialize
  mapM_ (\i -> do
    M.write_ alpha i (diag !. i)
    M.write_ beta i (if i < nn - 1 then supdiag !. i else 0)
    ) [0..nn-1]
  -- Copy b into result
  mapM_ (\i -> M.write_ mx i (b !. i)) [0..nn-1]

  -- LDLᵀ factorization and forward substitution combined
  -- for k = 2:n
  --   t = β(k-1), β(k-1) = t/α(k-1), α(k) = α(k) - t·β(k-1)
  mapM_ (\k -> do
    t <- M.readM beta (k - 1)
    ak1 <- M.readM alpha (k - 1)
    let bk1 = t / ak1
    M.write_ beta (k - 1) bk1
    ak <- M.readM alpha k
    M.write_ alpha k (ak - t * bk1)
    ) [1..nn-1]

  -- Forward substitution: b(k) = b(k) - β(k-1)·b(k-1)
  mapM_ (\k -> do
    bk <- M.readM mx k
    bk1 <- M.readM mx (k - 1)
    betaK1 <- M.readM beta (k - 1)
    M.write_ mx k (bk - betaK1 * bk1)
    ) [1..nn-1]

  -- Diagonal solve: b(n) = b(n)/α(n)
  bn <- M.readM mx (nn - 1)
  an <- M.readM alpha (nn - 1)
  M.write_ mx (nn - 1) (bn / an)

  -- Back substitution: b(k) = b(k)/α(k) - β(k)·b(k+1)
  mapM_ (\k -> do
    bk <- M.readM mx k
    ak <- M.readM alpha k
    betaK <- M.readM beta k
    bk1 <- M.readM mx (k + 1)
    M.write_ mx k (bk / ak - betaK * bk1)
    ) [nn-2, nn-3..0]
