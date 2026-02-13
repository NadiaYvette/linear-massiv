{-# LANGUAGE AllowAmbiguousTypes #-}

-- |
-- Module      : Numeric.LinearAlgebra.Massiv.Orthogonal.LeastSquares
-- Copyright   : (c) Nadia Chambers 2026
-- License     : BSD-3-Clause
-- Maintainer  : nadia.chambers@iohk.io
-- Stability   : experimental
--
-- Linear least squares solvers.
--
-- This module provides two methods for solving the linear least squares
-- problem
--
-- \( \min_x \| A x - b \|_2 \)
--
-- following Golub & Van Loan, /Matrix Computations/, 4th edition (GVL4),
-- Section 5.3, pp. 260--270.
--
-- __Theorem 5.3.1 (Least squares existence, GVL4 p. 260).__  Let
-- \( A \in \mathbb{R}^{m \times n} \) with \( m \ge n \) and
-- \( \operatorname{rank}(A) = n \).  Then the least squares problem
-- \( \min_x \| A x - b \|_2 \) has a unique solution \( x^* \) given
-- by the /normal equations/
--
-- \( A^T A \, x = A^T b \)
--
-- Two solution methods are provided:
--
-- * __QR-based__ ('leastSquaresQR') -- GVL4 Algorithm 5.3.2 (p. 262).
--   Factor \( A = Q R \) via Householder QR, then solve
--   \( R_1 x = (Q^T b)_{1:n} \) by back-substitution, where \( R_1 \)
--   denotes the leading \( n \times n \) upper triangular block of
--   \( R \).  This is the recommended method: it is numerically stable
--   and does not square the condition number.
--
-- * __Normal equations__ ('leastSquaresNormal') -- GVL4 Section 5.3.2
--   (p. 261).  Form \( A^T A \) and \( A^T b \) explicitly, then solve
--   via Cholesky factorisation.  This is faster but squares the
--   condition number: \( \kappa_2(A^T A) = \kappa_2(A)^2 \), so it
--   should only be used when \( A \) is well-conditioned.
--
-- __Complexity.__
--
-- * QR-based: \( 2mn^2 \) flops (dominated by the QR factorisation).
-- * Normal equations: \( mn^2 + \tfrac{1}{3}n^3 \) flops
--   (\( mn^2 \) for forming \( A^T A \), \( \tfrac{1}{3}n^3 \) for
--   Cholesky).
module Numeric.LinearAlgebra.Massiv.Orthogonal.LeastSquares
  ( -- * Least squares via QR
    leastSquaresQR
    -- * Normal equations
  , leastSquaresNormal
  ) where

import qualified Data.Massiv.Array as M
import Data.Massiv.Array (Ix2(..), Sz(..))
import GHC.TypeNats (KnownNat)

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.Orthogonal.QR (qr)
import Numeric.LinearAlgebra.Massiv.Solve.Triangular (backSub)
import Numeric.LinearAlgebra.Massiv.BLAS.Level3 (matMul, transpose)
import Numeric.LinearAlgebra.Massiv.BLAS.Level2 (matvec)
import Numeric.LinearAlgebra.Massiv.Solve.Cholesky (choleskySolve)

-- | Solve the least squares problem \( \min_x \| A x - b \|_2 \) via QR
-- factorisation (GVL4 Algorithm 5.3.2, p. 262).
--
-- Given \( A \in \mathbb{R}^{m \times n} \) with \( m \ge n \) and full
-- column rank, and a right-hand side \( b \in \mathbb{R}^m \), compute
-- the unique least squares solution \( x^* \in \mathbb{R}^n \).
--
-- The algorithm proceeds in three steps:
--
-- 1. Factor \( A = Q R \) using Householder QR ('qr').
-- 2. Form the transformed right-hand side \( Q^T b \).
-- 3. Solve the \( n \times n \) upper triangular system
--    \( R_1 x = (Q^T b)_{1:n} \) by back-substitution, where \( R_1 \)
--    is the leading \( n \times n \) block of \( R \).
--
-- This method is numerically stable because orthogonal transformations
-- preserve the 2-norm and do not amplify rounding errors.  The condition
-- number relevant to the solution is \( \kappa_2(A) \), not
-- \( \kappa_2(A)^2 \) as with the normal equations.
--
-- __Complexity:__ \( 2mn^2 - \tfrac{2}{3}n^3 \) flops for the QR
-- factorisation plus \( 2mn \) flops for forming \( Q^T b \) and
-- \( n^2 \) flops for back-substitution, giving a total of
-- \( O(2mn^2) \) flops.
leastSquaresQR :: forall m n r e. (KnownNat m, KnownNat n, M.Manifest r e, Floating e, Ord e)
               => Matrix m n r e -> Vector m r e -> Vector n r e
leastSquaresQR a b =
  let mm = dimVal @m
      nn = dimVal @n
      (q, r) = qr a
      qt = transpose q
      -- qtb = Qᵀ·b (dimension m)
      qtb = matvec qt b
      -- Extract top n×n of R (which is upper triangular)
      r1 = makeMatrix @n @n @r $ \i j -> r ! (i, j)
      -- Extract top n entries of Qᵀb
      qtb1 = makeVector @n @r $ \i -> qtb !. i
  in backSub r1 qtb1

-- | Solve the least squares problem \( \min_x \| A x - b \|_2 \) via the
-- normal equations (GVL4 Section 5.3.2, p. 261).
--
-- Given \( A \in \mathbb{R}^{m \times n} \) with full column rank and
-- \( b \in \mathbb{R}^m \), form the \( n \times n \) symmetric positive
-- definite system
--
-- \( A^T A \, x = A^T b \)
--
-- and solve it using Cholesky factorisation (\( A^T A = L L^T \)).
--
-- __Warning:__ The normal equations square the condition number of \( A \):
-- \( \kappa_2(A^T A) = \kappa_2(A)^2 \).  For ill-conditioned problems
-- this leads to a significant loss of accuracy compared to QR-based
-- methods.  Prefer 'leastSquaresQR' unless \( A \) is known to be
-- well-conditioned (GVL4 p. 261).
--
-- __Complexity:__ \( mn^2 \) flops for forming \( A^T A \), plus
-- \( \tfrac{1}{3}n^3 \) flops for the Cholesky factorisation, giving a
-- total of \( O(mn^2 + \tfrac{1}{3}n^3) \) flops.
leastSquaresNormal :: forall m n r e. (KnownNat m, KnownNat n, M.Manifest r e, Floating e, Ord e)
                   => Matrix m n r e -> Vector m r e -> Vector n r e
leastSquaresNormal a b =
  let at = transpose a
      ata = matMul at a       -- n×n, symmetric positive definite
      atb = matvec at b       -- n×1
  in choleskySolve ata atb
