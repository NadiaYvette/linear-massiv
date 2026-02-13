{-# LANGUAGE AllowAmbiguousTypes #-}

-- |
-- Module      : Numeric.LinearAlgebra.Massiv.Eigen.Power
-- Copyright   : (c) Nadia Chambers 2026
-- License     : BSD-3-Clause
-- Maintainer  : nadia.chambers@iohk.io
-- Stability   : experimental
--
-- Power iteration methods for computing individual eigenvalue\/eigenvector
-- pairs of a general square matrix.
--
-- This module implements three iterative projection techniques drawn from
-- Golub & Van Loan, /Matrix Computations/, 4th edition (GVL4), Section 7.3,
-- pp. 372--382:
--
-- * __Power method__ (Algorithm 7.3.3, p. 375) — converges to the dominant
--   eigenpair at a rate governed by the ratio \(|\lambda_2 / \lambda_1|\) per
--   iteration, where \(\lambda_1\) is the eigenvalue of largest modulus.
--
-- * __Inverse iteration__ (Section 7.3.1, p. 377) — given a shift \(\mu\),
--   converges to the eigenvalue closest to \(\mu\) by applying the power
--   method to \((A - \mu I)^{-1}\).
--
-- * __Rayleigh quotient iteration__ (Section 7.3.2, p. 379) — an adaptive
--   variant of inverse iteration in which the shift is updated at every step
--   to equal the current Rayleigh quotient.  For symmetric matrices this
--   achieves /cubic/ convergence; for general matrices the convergence is
--   /quadratic/.
--
-- All three routines return an approximate eigenvalue \(\lambda\) and its
-- associated eigenvector \(q\) satisfying \(Aq \approx \lambda q\).
module Numeric.LinearAlgebra.Massiv.Eigen.Power
  ( -- * Power method (Algorithm 7.3.3)
    powerMethod
    -- * Inverse iteration (Section 7.3.1)
  , inverseIteration
    -- * Rayleigh quotient iteration (Section 7.3.2)
  , rayleighQuotient
  ) where

import qualified Data.Massiv.Array as M
import Data.Massiv.Array (Ix2(..), Sz(..))
import GHC.TypeNats (KnownNat)

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.BLAS.Level1 (dot, scal, nrm2)
import Numeric.LinearAlgebra.Massiv.BLAS.Level2 (matvec)
import Numeric.LinearAlgebra.Massiv.Solve.LU (luSolve)

-- | Power method for computing the dominant eigenpair (GVL4 Algorithm 7.3.3,
-- p. 375).
--
-- Given a square matrix \(A \in \mathbb{R}^{n \times n}\) with eigenvalues
-- ordered \(|\lambda_1| > |\lambda_2| \geq \cdots \geq |\lambda_n|\), the
-- power method generates a sequence of vectors
--
-- \[
--   z_k = A q_{k-1}, \qquad q_k = z_k / \|z_k\|_2
-- \]
--
-- that converges to the eigenvector associated with \(\lambda_1\).  The
-- corresponding eigenvalue is estimated via the Rayleigh quotient
-- \(\lambda \approx q_k^T A q_k\).
--
-- __Convergence rate:__ the error contracts by a factor of
-- \(|\lambda_2 / \lambda_1|\) per iteration (GVL4, p. 375).  The method
-- therefore requires a dominant eigenvalue that is well-separated from the
-- rest of the spectrum.
--
-- Returns @(eigenvalue, eigenvector)@ once the eigenvalue estimate changes by
-- less than the given tolerance, or after the specified number of iterations.
powerMethod :: forall n r e. (KnownNat n, M.Manifest r e, Floating e, Ord e)
            => Matrix n n r e
            -> Vector n r e    -- ^ Initial guess \(q_0\) (should be unit norm)
            -> Int             -- ^ Maximum iterations
            -> e               -- ^ Convergence tolerance
            -> (e, Vector n r e)
powerMethod a q0 maxIter tol = go 0 q0 0
  where
    go :: Int -> Vector n r e -> e -> (e, Vector n r e)
    go iter q prevLambda
      | iter >= maxIter = (prevLambda, q)
      | otherwise =
        let z = matvec a q                        -- z = A·q
            znorm = nrm2 z                        -- ‖z‖₂
            qNew = scal (1 / znorm) z             -- q = z / ‖z‖₂
            lambda = dot qNew (matvec a qNew)     -- λ = qᵀAq (Rayleigh quotient)
        in if abs (lambda - prevLambda) < tol
           then (lambda, qNew)
           else go (iter + 1) qNew lambda

-- | Inverse iteration for computing the eigenpair closest to a given shift
-- (GVL4 Section 7.3.1, p. 377).
--
-- Given a shift \(\mu\) that approximates some eigenvalue of \(A\), inverse
-- iteration applies the power method to the matrix \((A - \mu I)^{-1}\).
-- Each step solves the linear system
--
-- \[
--   (A - \mu I)\, z_k = q_{k-1}, \qquad q_k = z_k / \|z_k\|_2
-- \]
--
-- and converges to the eigenvalue \(\lambda_j\) that minimises
-- \(|\lambda_j - \mu|\).  The convergence rate is
-- \(|\lambda_j - \mu| / |\lambda_i - \mu|\) per iteration, where
-- \(\lambda_i\) is the second-closest eigenvalue to \(\mu\).
--
-- The eigenvalue estimate is refined at each step via the Rayleigh quotient
-- \(\lambda \approx q_k^T A q_k\) (using the /original/ matrix \(A\)).
--
-- Returns @(eigenvalue, eigenvector)@.
inverseIteration :: forall n r e. (KnownNat n, M.Manifest r e, Floating e, Ord e)
                 => Matrix n n r e
                 -> e               -- ^ Shift \(\mu\)
                 -> Vector n r e    -- ^ Initial guess
                 -> Int             -- ^ Maximum iterations
                 -> e               -- ^ Convergence tolerance
                 -> (e, Vector n r e)
inverseIteration a mu q0 maxIter tol = go 0 q0 0
  where
    aShifted = makeMatrix @n @n @r $ \i j ->
      if i == j then (a ! (i, j)) - mu else a ! (i, j)

    go :: Int -> Vector n r e -> e -> (e, Vector n r e)
    go iter q prevLambda
      | iter >= maxIter = (prevLambda, q)
      | otherwise =
        let z = luSolve aShifted q                -- Solve (A - μI)z = q
            znorm = nrm2 z
            qNew = scal (1 / znorm) z
            lambda = dot qNew (matvec a qNew)     -- Rayleigh quotient with original A
        in if abs (lambda - prevLambda) < tol
           then (lambda, qNew)
           else go (iter + 1) qNew lambda

-- | Rayleigh quotient iteration (GVL4 Section 7.3.2, p. 379).
--
-- An adaptive variant of inverse iteration in which the shift \(\mu_k\) is
-- set equal to the current Rayleigh quotient at every step:
--
-- \[
--   \mu_k = q_k^T A q_k, \qquad (A - \mu_k I)\, z_{k+1} = q_k, \qquad
--   q_{k+1} = z_{k+1} / \|z_{k+1}\|_2
-- \]
--
-- __Convergence:__
--
-- * For /symmetric/ matrices \(A = A^T\), the iteration converges
--   /cubically/ — the residual \(\|Aq - \lambda q\|\) is cubed at each step
--   (GVL4, p. 379).
-- * For general (non-symmetric) matrices, convergence is /quadratic/.
--
-- Because the shift changes at each iteration, a fresh LU factorisation of
-- \(A - \mu_k I\) is computed every step.  Despite this extra cost the rapid
-- convergence usually makes Rayleigh quotient iteration the method of choice
-- when a good initial vector is available.
--
-- Returns @(eigenvalue, eigenvector)@.
rayleighQuotient :: forall n r e. (KnownNat n, M.Manifest r e, Floating e, Ord e)
                 => Matrix n n r e
                 -> Vector n r e    -- ^ Initial guess
                 -> Int             -- ^ Maximum iterations
                 -> e               -- ^ Convergence tolerance
                 -> (e, Vector n r e)
rayleighQuotient a q0 maxIter tol = go 0 q0 (dot q0 (matvec a q0))
  where
    go :: Int -> Vector n r e -> e -> (e, Vector n r e)
    go iter q lambda
      | iter >= maxIter = (lambda, q)
      | otherwise =
        let nn = dimVal @n
            aShifted = makeMatrix @n @n @r $ \i j ->
              if i == j then (a ! (i, j)) - lambda else a ! (i, j)
            z = luSolve aShifted q
            znorm = nrm2 z
            qNew = scal (1 / znorm) z
            lambdaNew = dot qNew (matvec a qNew)
        in if abs (lambdaNew - lambda) < tol
           then (lambdaNew, qNew)
           else go (iter + 1) qNew lambdaNew
