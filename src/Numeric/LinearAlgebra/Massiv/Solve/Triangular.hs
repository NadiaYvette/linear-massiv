{-# LANGUAGE AllowAmbiguousTypes #-}

-- |
-- Module      : Numeric.LinearAlgebra.Massiv.Solve.Triangular
-- Copyright   : (c) Nadia Chambers 2026
-- License     : BSD-3-Clause
-- Maintainer  : nadia.chambers@iohk.io
-- Stability   : experimental
--
-- = Triangular System Solvers
--
-- Forward and back substitution for lower- and upper-triangular linear
-- systems, following Golub & Van Loan, /Matrix Computations/, 4th edition
-- (GVL4), Section 3.1, pp. 106--113.
--
-- The section presents four algorithmic variants organised along two axes:
--
-- * __Row-oriented vs. column-oriented__ inner loops, affecting data-access
--   patterns and cache behaviour.
-- * __General vs. unit-triangular__ coefficient matrices, where the unit
--   variants avoid division by the (known-to-be-one) diagonal.
--
-- This module exposes the following mapping:
--
-- +-------------------+---------------------------+---------------------------------+
-- | Function          | Algorithm                 | Reference                       |
-- +===================+===========================+=================================+
-- | 'forwardSub'      | Row-oriented forward sub  | GVL4 Algorithm 3.1.1, p. 106    |
-- +-------------------+---------------------------+---------------------------------+
-- | 'backSub'         | Row-oriented back sub     | GVL4 Algorithm 3.1.2, p. 107    |
-- +-------------------+---------------------------+---------------------------------+
-- | 'forwardSubUnit'  | Column-oriented fwd sub   | GVL4 Algorithm 3.1.3, p. 108    |
-- +-------------------+---------------------------+---------------------------------+
-- | 'backSubUnit'     | Column-oriented back sub  | GVL4 Algorithm 3.1.4, p. 109    |
-- +-------------------+---------------------------+---------------------------------+
--
-- == Complexity
--
-- Every solver performs \(O(n^2 / 2)\) floating-point operations (flops) for
-- an \(n \times n\) triangular system (GVL4 p. 109).
--
-- == Type Safety
--
-- The matrix dimension \(n\) is tracked at the type level via 'KnownNat',
-- so the compiler statically guarantees that the coefficient matrix is
-- square and that the right-hand side vector has a conforming length.
module Numeric.LinearAlgebra.Massiv.Solve.Triangular
  ( -- * Forward substitution (\(Lx = b\))
    forwardSub
    -- * Back substitution (\(Ux = b\))
  , backSub
    -- * Unit-triangular variants
  , forwardSubUnit
  , backSubUnit
  ) where

import qualified Data.Massiv.Array as M
import Data.Massiv.Array (Ix1, Ix2(..), Sz(..))
import GHC.TypeNats (KnownNat)

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal

-- | Row-oriented forward substitution (GVL4 Algorithm 3.1.1, p. 106).
--
-- Solves the lower-triangular system \(Lx = b\) where \(L \in \mathbb{R}^{n \times n}\)
-- is lower triangular with nonzero diagonal entries.
--
-- ==== Mathematical definition
--
-- For \(j = 1, \ldots, n\):
--
-- \[
--   x_j = \frac{1}{\ell_{jj}} \left( b_j - \sum_{k=1}^{j-1} \ell_{jk}\, x_k \right)
-- \]
--
-- ==== Type-safety guarantees
--
-- The type-level natural \(n\) ('KnownNat') ensures that \(L\) is square
-- and that \(b\) has exactly \(n\) entries. A dimension mismatch is a
-- compile-time error.
--
-- ==== Complexity
--
-- \(O(n^2 / 2)\) flops (GVL4 p. 109).
--
-- ==== Reference
--
-- Golub & Van Loan, /Matrix Computations/, 4th ed., Algorithm 3.1.1
-- (Row-Oriented Forward Substitution), p. 106.
forwardSub :: forall n r e. (KnownNat n, M.Manifest r e, Fractional e)
           => Matrix n n r e -> Vector n r e -> Vector n r e
forwardSub l b = createVector @n $ \mx -> do
  let nn = dimVal @n
  -- Copy b into the mutable result
  mapM_ (\i -> M.write_ mx i (b !. i)) [0..nn-1]
  -- Forward elimination
  mapM_ (\j -> do
    xj <- M.readM mx j
    let ldiag = l ! (j, j)
        xj' = xj / ldiag
    M.write_ mx j xj'
    -- Update remaining entries
    mapM_ (\i -> do
      xi <- M.readM mx i
      M.write_ mx i (xi - (l ! (i, j)) * xj')
      ) [j+1..nn-1]
    ) [0..nn-1]

-- | Row-oriented back substitution (GVL4 Algorithm 3.1.2, p. 107).
--
-- Solves the upper-triangular system \(Ux = b\) where \(U \in \mathbb{R}^{n \times n}\)
-- is upper triangular with nonzero diagonal entries.
--
-- ==== Mathematical definition
--
-- For \(j = n, n-1, \ldots, 1\):
--
-- \[
--   x_j = \frac{1}{u_{jj}} \left( b_j - \sum_{k=j+1}^{n} u_{jk}\, x_k \right)
-- \]
--
-- ==== Type-safety guarantees
--
-- The type-level natural \(n\) ('KnownNat') ensures that \(U\) is square
-- and that \(b\) has exactly \(n\) entries. A dimension mismatch is a
-- compile-time error.
--
-- ==== Complexity
--
-- \(O(n^2 / 2)\) flops (GVL4 p. 109).
--
-- ==== Reference
--
-- Golub & Van Loan, /Matrix Computations/, 4th ed., Algorithm 3.1.2
-- (Row-Oriented Back Substitution), p. 107.
backSub :: forall n r e. (KnownNat n, M.Manifest r e, Fractional e)
        => Matrix n n r e -> Vector n r e -> Vector n r e
backSub u b = createVector @n $ \mx -> do
  let nn = dimVal @n
  -- Copy b into the mutable result
  mapM_ (\i -> M.write_ mx i (b !. i)) [0..nn-1]
  -- Back elimination
  mapM_ (\j -> do
    xj <- M.readM mx j
    let udiag = u ! (j, j)
        xj' = xj / udiag
    M.write_ mx j xj'
    -- Update remaining entries
    mapM_ (\i -> do
      xi <- M.readM mx i
      M.write_ mx i (xi - (u ! (i, j)) * xj')
      ) [0..j-1]
    ) [nn-1, nn-2..0]

-- | Column-oriented forward substitution for unit lower triangular systems
-- (GVL4 Algorithm 3.1.3, p. 108).
--
-- Solves \(Lx = b\) where \(L \in \mathbb{R}^{n \times n}\) is /unit/ lower
-- triangular, i.e. \(\ell_{jj} = 1\) for all \(j\). Because the diagonal is
-- implicitly one, no division is needed and the constraint relaxes from
-- 'Fractional' to 'Num'.
--
-- ==== Mathematical definition
--
-- For \(j = 1, \ldots, n\):
--
-- \[
--   x_j = b_j - \sum_{k=1}^{j-1} \ell_{jk}\, x_k
-- \]
--
-- The implementation uses a /column-oriented/ (saxpy) loop: once \(x_j\) is
-- determined, rows \(i > j\) are updated by subtracting \(\ell_{ij}\, x_j\).
--
-- ==== Type-safety guarantees
--
-- Identical to 'forwardSub': the dimensions are enforced at compile time
-- via 'KnownNat'.
--
-- ==== Complexity
--
-- \(O(n^2 / 2)\) flops (GVL4 p. 109).
--
-- ==== Reference
--
-- Golub & Van Loan, /Matrix Computations/, 4th ed., Algorithm 3.1.3
-- (Column-Oriented Forward Substitution), p. 108.
forwardSubUnit :: forall n r e. (KnownNat n, M.Manifest r e, Num e)
               => Matrix n n r e -> Vector n r e -> Vector n r e
forwardSubUnit l b = createVector @n $ \mx -> do
  let nn = dimVal @n
  -- Copy b into the mutable result
  mapM_ (\i -> M.write_ mx i (b !. i)) [0..nn-1]
  -- Column-oriented forward substitution
  mapM_ (\j -> do
    xj <- M.readM mx j
    -- Subtract l(i,j)*xj from x(i) for i > j
    mapM_ (\i -> do
      xi <- M.readM mx i
      M.write_ mx i (xi - (l ! (i, j)) * xj)
      ) [j+1..nn-1]
    ) [0..nn-1]

-- | Column-oriented back substitution for unit upper triangular systems
-- (GVL4 Algorithm 3.1.4, p. 109).
--
-- Solves \(Ux = b\) where \(U \in \mathbb{R}^{n \times n}\) is /unit/ upper
-- triangular, i.e. \(u_{jj} = 1\) for all \(j\). Because the diagonal is
-- implicitly one, no division is needed and the constraint relaxes from
-- 'Fractional' to 'Num'.
--
-- ==== Mathematical definition
--
-- For \(j = n, n-1, \ldots, 1\):
--
-- \[
--   x_j = b_j - \sum_{k=j+1}^{n} u_{jk}\, x_k
-- \]
--
-- The implementation uses a /column-oriented/ loop: once \(x_j\) is known,
-- rows \(i < j\) are updated by subtracting \(u_{ij}\, x_j\).
--
-- ==== Type-safety guarantees
--
-- Identical to 'backSub': the dimensions are enforced at compile time via
-- 'KnownNat'.
--
-- ==== Complexity
--
-- \(O(n^2 / 2)\) flops (GVL4 p. 109).
--
-- ==== Reference
--
-- Golub & Van Loan, /Matrix Computations/, 4th ed., Algorithm 3.1.4
-- (Column-Oriented Back Substitution), p. 109.
backSubUnit :: forall n r e. (KnownNat n, M.Manifest r e, Num e)
            => Matrix n n r e -> Vector n r e -> Vector n r e
backSubUnit u b = createVector @n $ \mx -> do
  let nn = dimVal @n
  -- Copy b into the mutable result
  mapM_ (\i -> M.write_ mx i (b !. i)) [0..nn-1]
  -- Column-oriented back substitution
  mapM_ (\j -> do
    xj <- M.readM mx j
    -- Subtract u(i,j)*xj from x(i) for i < j
    mapM_ (\i -> do
      xi <- M.readM mx i
      M.write_ mx i (xi - (u ! (i, j)) * xj)
      ) [0..j-1]
    ) [nn-1, nn-2..0]
