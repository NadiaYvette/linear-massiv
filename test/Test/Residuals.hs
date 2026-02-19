{-# LANGUAGE AllowAmbiguousTypes #-}

-- | LAPACK-style scaled residual functions for NLA verification.
--
-- This module provides the standard residual metrics used by LAPACK
-- (Anderson et al., 1999) and recommended by Higham (2002) for
-- verifying numerical linear algebra routines. Rather than comparing
-- element-wise against a fixed tolerance, these functions compute
-- /scaled residuals/ that account for problem size, matrix norms,
-- and machine precision.
--
-- A scaled residual less than \(O(1)\) indicates backward stability;
-- values less than ~10--100 are considered acceptable.
--
-- __References:__
--
-- * Higham, N. J. (2002). /Accuracy and Stability of Numerical
--   Algorithms/, 2nd ed., SIAM. Chapter 1.
-- * Anderson, E. et al. (1999). /LAPACK Users' Guide/, 3rd ed., SIAM.
-- * LAPACK Working Note 41: Installation Guide.
-- * Golub, G. H. & Van Loan, C. F. (2013). /Matrix Computations/,
--   4th ed., Chapter 2.
module Test.Residuals
  ( -- * Machine epsilon
    machineEps
    -- * Scaled residual functions
  , scaledResidualLinear
  , scaledResidualEigen
  , scaledResidualQR
  , scaledResidualSVD
  , scaledResidualCholesky
  , scaledResidualLU
  , orthogonalityResidual
    -- * Condition-number-aware tolerance
  , conditionTolerance
  , safeCondition2
  ) where

import qualified Data.Massiv.Array as M
import Data.Massiv.Array (Ix1, Ix2(..), Array)
import GHC.TypeNats (KnownNat)

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.BLAS.Level1 (scal)
import Numeric.LinearAlgebra.Massiv.BLAS.Level2 (matvec)
import Numeric.LinearAlgebra.Massiv.BLAS.Level3 (matMul, transpose, mSub)
import Numeric.LinearAlgebra.Massiv.Norms (normFrob, normInf, vnorm2, vnormInf)
import Numeric.LinearAlgebra.Massiv.Eigen.SVD (singularValues)

-- | IEEE 754 double precision machine epsilon (unit roundoff).
--
-- \(\varepsilon = 2^{-52} \approx 2.22 \times 10^{-16}\)
--
-- This is the smallest value such that \(1 + \varepsilon > 1\) in
-- double precision floating-point arithmetic.
machineEps :: Double
machineEps = 2.220446049250313e-16

-- | Scaled residual for a linear system solve \(Ax = b\).
--
-- \[
-- \eta(\hat{x}) = \frac{\|A\hat{x} - b\|_\infty}
--                      {\|A\|_\infty \|\hat{x}\|_\infty + \|b\|_\infty}
-- \]
--
-- Per Higham (2002), Section 7.1, this is the normwise backward error.
-- A backward-stable solver should produce \(\eta \le O(n \varepsilon)\).
scaledResidualLinear :: forall n. KnownNat n
  => Matrix n n M.P Double -> Vector n M.P Double -> Vector n M.P Double -> Double
scaledResidualLinear a x b =
  let ax = matvec a x
      nn = dimVal @n
      -- residual r = Ax - b
      r = makeVector @n @M.P $ \i -> (ax !. i) - (b !. i)
      numr = vnormInf r
      denom = normInf a * vnormInf x + vnormInf b + fromIntegral nn * machineEps
  in numr / denom

-- | Scaled residual for an eigenpair \((\\lambda, v)\).
--
-- \[
-- \frac{\|Av - \lambda v\|_2}{\|A\|_F \|v\|_2}
-- \]
--
-- Per Higham (2002), Section 14.1. A well-computed eigenpair
-- should have residual \(O(n \varepsilon)\).
scaledResidualEigen :: forall n. KnownNat n
  => Matrix n n M.P Double -> Double -> Vector n M.P Double -> Double
scaledResidualEigen a lambda v =
  let av = matvec a v
      lambdaV = scal lambda v
      r = makeVector @n @M.P $ \i -> (av !. i) - (lambdaV !. i)
      numr = vnorm2 r
      denom = normFrob a * vnorm2 v + machineEps
  in numr / denom

-- | Scaled residual for QR factorization \(A = QR\).
--
-- \[
-- \frac{\|A - QR\|_F}{\|A\|_F \cdot n \cdot \varepsilon}
-- \]
--
-- Per LAPACK testing methodology (LAWN 41). A value less than
-- \(O(1)\) (typically < 10--100) indicates the factorization
-- is backward stable.
scaledResidualQR :: forall m n. (KnownNat m, KnownNat n)
  => Matrix m n M.P Double -> Matrix m m M.P Double -> Matrix m n M.P Double -> Double
scaledResidualQR a q r =
  let qr_prod = matMul q r
      diff = mSub a qr_prod
      nn = dimVal @n
      numr = normFrob diff
      denom = normFrob a * fromIntegral nn * machineEps + machineEps
  in numr / denom

-- | Scaled residual for SVD \(A = U \Sigma V^T\).
--
-- \[
-- \frac{\|A - U \Sigma V^T\|_F}{\|A\|_F \cdot \max(m,n) \cdot \varepsilon}
-- \]
--
-- Per LAPACK testing methodology (LAWN 41).
scaledResidualSVD :: forall m n. (KnownNat m, KnownNat n)
  => Matrix m n M.P Double -> Matrix m m M.P Double -> Vector n M.P Double
  -> Matrix n n M.P Double -> Double
scaledResidualSVD a u sigma v =
  let mm = dimVal @m
      nn = dimVal @n
      -- Build U * diag(sigma) by scaling columns of U
      -- Then multiply by V^T
      -- U is m x m, sigma has n entries, V is n x n
      -- We need U(:,0:n-1) * diag(sigma) * V^T
      -- Since our SVD returns m x m U, we take the first n columns conceptually
      uSigma = makeMatrix @m @n @M.P $ \i j ->
        (u ! (i, j)) * (sigma !. j)
      uSigmaVt = matMul uSigma (transpose v)
      diff = mSub a uSigmaVt
      numr = normFrob diff
      denom = normFrob a * fromIntegral (max mm nn) * machineEps + machineEps
  in numr / denom

-- | Scaled residual for Cholesky factorization \(A = GG^T\).
--
-- \[
-- \frac{\|A - GG^T\|_F}{\|A\|_F \cdot n \cdot \varepsilon}
-- \]
scaledResidualCholesky :: forall n. KnownNat n
  => Matrix n n M.P Double -> Matrix n n M.P Double -> Double
scaledResidualCholesky a g =
  let ggt = matMul g (transpose g)
      diff = mSub a ggt
      nn = dimVal @n
      numr = normFrob diff
      denom = normFrob a * fromIntegral nn * machineEps + machineEps
  in numr / denom

-- | Scaled residual for LU factorization \(PA = LU\).
--
-- \[
-- \frac{\|PA - LU\|_F}{\|A\|_F \cdot n \cdot \varepsilon}
-- \]
--
-- Takes the original matrix, the packed LU matrix, and the pivot array.
-- Extracts L (unit lower triangular) and U (upper triangular) from the
-- packed form.
scaledResidualLU :: forall n. KnownNat n
  => Matrix n n M.P Double -> Matrix n n M.P Double
  -> Array M.P Ix1 Int -> Double
scaledResidualLU a lu_packed pivots =
  let nn = dimVal @n
      -- Extract L (unit lower triangular)
      l = makeMatrix @n @n @M.P $ \i j ->
        if i == j then 1
        else if i > j then lu_packed ! (i, j)
        else 0
      -- Extract U (upper triangular)
      u_mat = makeMatrix @n @n @M.P $ \i j ->
        if i <= j then lu_packed ! (i, j)
        else 0
      -- Construct PA by applying the permutation
      pa = makeMatrix @n @n @M.P $ \i j ->
        let pi_i = M.index' pivots i
        in a ! (pi_i, j)
      lu_prod = matMul l u_mat
      diff = mSub pa lu_prod
      numr = normFrob diff
      denom = normFrob a * fromIntegral nn * machineEps + machineEps
  in numr / denom

-- | Orthogonality residual for a matrix \(Q\).
--
-- \[
-- \frac{\|Q^T Q - I\|_F}{n \cdot \varepsilon}
-- \]
--
-- Per LAPACK testing methodology. A value less than \(O(1)\)
-- means orthogonality is preserved to within rounding error.
orthogonalityResidual :: forall n. KnownNat n
  => Matrix n n M.P Double -> Double
orthogonalityResidual q =
  let nn = dimVal @n
      qtq = matMul (transpose q) q
      diff = mSub qtq (identityMatrix @n @M.P)
      numr = normFrob diff
      denom = fromIntegral nn * machineEps
  in numr / denom

-- | Condition-number-aware tolerance: \(\kappa \cdot n \cdot \varepsilon\).
--
-- For a linear system \(Ax = b\), the forward error satisfies
-- \(\|\hat{x} - x\| / \|x\| \le \kappa(A) \cdot \eta\) where
-- \(\eta\) is the backward error. A backward-stable algorithm
-- achieves \(\eta \approx n \varepsilon\), so the expected
-- forward error is \(O(\kappa \cdot n \cdot \varepsilon)\).
--
-- See Higham (2002), Theorem 7.2.
conditionTolerance :: Double -> Int -> Double
conditionTolerance kappa n = kappa * fromIntegral n * machineEps

-- | Estimate the 2-norm condition number \(\kappa_2(A) = \sigma_{\max} / \sigma_{\min}\).
--
-- Uses SVD to compute singular values. Returns \(10^{16}\) for
-- numerically singular matrices (where \(\sigma_{\min} < \varepsilon \cdot \sigma_{\max}\)).
safeCondition2 :: forall n. KnownNat n => Matrix n n M.P Double -> Double
safeCondition2 a =
  let sv = singularValues a
      nn = dimVal @n
      sigmaMax = sv !. 0
      sigmaMin = sv !. (nn - 1)
  in if abs sigmaMin < machineEps * abs sigmaMax
     then 1e16
     else abs sigmaMax / abs sigmaMin
