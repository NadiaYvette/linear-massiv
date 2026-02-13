{-# LANGUAGE AllowAmbiguousTypes #-}

-- |
-- Module      : Numeric.LinearAlgebra.Massiv.Eigen.Schur
-- Copyright   : (c) Nadia Chambers 2026
-- License     : BSD-3-Clause
-- Maintainer  : nadia.chambers@iohk.io
-- Stability   : experimental
--
-- Real Schur decomposition via the practical QR algorithm, following
-- Golub & Van Loan, /Matrix Computations/, 4th edition (GVL4), Section 7.5,
-- pp. 393--417.
--
-- __Theorem 7.5.1 (Real Schur Decomposition, p. 393):__ For every
-- \(A \in \mathbb{R}^{n \times n}\) there exists an orthogonal matrix \(Q\)
-- such that
--
-- \[
--   A = Q \, T \, Q^T
-- \]
--
-- where \(T\) is upper /quasi/-triangular: its diagonal consists of \(1
-- \times 1\) blocks (real eigenvalues) and \(2 \times 2\) blocks whose
-- eigenvalues are complex conjugate pairs \(\alpha \pm \beta i\).
--
-- __Algorithm:__ The implementation follows GVL4 Algorithm 7.5.1 (Practical
-- QR Algorithm, p. 395):
--
--   1. Reduce \(A\) to upper Hessenberg form \(H\) via
--      "Numeric.LinearAlgebra.Massiv.Eigen.Hessenberg".
--   2. Apply implicit single-shift QR iterations with Givens rotations on
--      \(H\), using the /Wilkinson shift/ (eigenvalue of the trailing \(2
--      \times 2\) block closest to \(h_{nn}\), p. 397) to accelerate
--      convergence.
--   3. Deflate converged eigenvalues from the bottom of the active
--      Hessenberg window.
--
-- The Wilkinson shift ensures global convergence; in practice, most
-- eigenvalues converge in only one or two iterations (GVL4, p. 397).
module Numeric.LinearAlgebra.Massiv.Eigen.Schur
  ( -- * Schur decomposition (Algorithm 7.5.1)
    schur
    -- * Eigenvalues from Schur form
  , eigenvalues
  ) where

import qualified Data.Massiv.Array as M
import Data.Massiv.Array (Ix2(..), Sz(..))
import GHC.TypeNats (KnownNat)

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.Eigen.Hessenberg (hessenberg)
import Numeric.LinearAlgebra.Massiv.Orthogonal.Givens (givensRotation)

-- | Real Schur decomposition (GVL4 Theorem 7.5.1, p. 393; Algorithm 7.5.1,
-- p. 395).
--
-- Computes orthogonal \(Q\) and upper quasi-triangular \(T\) satisfying
--
-- \[
--   A = Q \, T \, Q^T
-- \]
--
-- The matrix \(T\) has the same eigenvalues as \(A\).  Its diagonal blocks
-- are either:
--
--   * \(1 \times 1\) — corresponding to a real eigenvalue, or
--   * \(2 \times 2\) — corresponding to a pair of complex conjugate
--     eigenvalues \(\alpha \pm \beta i\).
--
-- Internally the algorithm first reduces \(A\) to upper Hessenberg form via
-- 'Numeric.LinearAlgebra.Massiv.Eigen.Hessenberg.hessenberg', then applies
-- implicit single-shift QR iterations using the /Wilkinson shift/ (GVL4,
-- p. 397) and Givens rotations.
--
-- Returns @(Q, T)@.
schur :: forall n r e. (KnownNat n, M.Manifest r e, Floating e, Ord e)
      => Matrix n n r e
      -> Int   -- ^ Maximum iterations
      -> e     -- ^ Convergence tolerance
      -> (Matrix n n r e, Matrix n n r e)
schur a maxIter tol =
  let nn = dimVal @n
      -- Step 1: Reduce to Hessenberg form
      (q0, h0) = hessenberg a
      -- Step 2: QR iteration on Hessenberg matrix
      (qFinal, tFinal) = qrIteration nn q0 h0 maxIter tol
  in (qFinal, tFinal)

-- | Implicit QR iteration on an upper Hessenberg matrix.
qrIteration :: forall n r e. (KnownNat n, M.Manifest r e, Floating e, Ord e)
            => Int -> Matrix n n r e -> Matrix n n r e -> Int -> e
            -> (Matrix n n r e, Matrix n n r e)
qrIteration nn q h maxIter tol = go 0 q h (nn - 1)
  where
    go :: Int -> Matrix n n r e -> Matrix n n r e -> Int -> (Matrix n n r e, Matrix n n r e)
    go iter q_ h_ p
      | iter >= maxIter = (q_, h_)
      | p <= 0 = (q_, h_)
      | otherwise =
        -- Check for convergence of h(p, p-1)
        let subdiag = abs (h_ ! (p, p - 1))
            diagSum = abs (h_ ! (p - 1, p - 1)) + abs (h_ ! (p, p))
        in if subdiag <= tol * diagSum
           then
             -- Deflate: set subdiagonal to zero, reduce problem size
             let h_new = makeMatrix @n @n @r $ \i j ->
                   if i == p && j == p - 1 then 0 else h_ ! (i, j)
             in go iter q_ h_new (p - 1)
           else
             -- Apply one QR step with Wilkinson shift
             let shift = wilkinsonShift (h_ ! (p-1, p-1)) (h_ ! (p-1, p))
                                         (h_ ! (p, p-1))   (h_ ! (p, p))
                 -- Shifted QR step: H - σI = QR, H_new = RQ + σI
                 -- Implemented via Givens rotations on Hessenberg matrix
                 (q_new, h_new) = qrStepGivens q_ h_ shift p
             in go (iter + 1) q_new h_new p

-- | Wilkinson shift (GVL4, p. 397).
--
-- Given the trailing \(2 \times 2\) block
--
-- \[
--   \begin{bmatrix} a & b \\ c & d \end{bmatrix}
-- \]
--
-- the Wilkinson shift is the eigenvalue of this block that is closest to
-- \(d\) (the bottom-right entry).  When the eigenvalues of the block are
-- complex the shift defaults to \(d\).
wilkinsonShift :: (Floating e, Ord e) => e -> e -> e -> e -> e
wilkinsonShift a b c d =
  let trace_ = a + d
      det_ = a * d - b * c
      disc = trace_ * trace_ / 4 - det_
  in if disc < 0
     then d  -- Complex eigenvalues; use d as shift
     else
       let sqrtDisc = sqrt disc
           mu1 = trace_ / 2 + sqrtDisc
           mu2 = trace_ / 2 - sqrtDisc
       in if abs (mu1 - d) < abs (mu2 - d) then mu1 else mu2

-- | One QR step on Hessenberg matrix using Givens rotations.
-- H ← shift, QR factorize, then H = RQ + shift.
qrStepGivens :: forall n r e. (KnownNat n, M.Manifest r e, Floating e, Ord e)
             => Matrix n n r e -> Matrix n n r e -> e -> Int
             -> (Matrix n n r e, Matrix n n r e)
qrStepGivens q h shift p =
  let nn = dimVal @n
      -- Apply shift: H ← H - σI
      h_shifted = makeMatrix @n @n @r $ \i j ->
        if i == j then (h ! (i, j)) - shift else h ! (i, j)
      -- QR factorization via Givens rotations (only on the active part)
      (rotations, r) = applyGivensQR h_shifted p
      -- Form RQ + σI
      h_new = formRQ r rotations shift p
      -- Update Q
      q_new = updateQ q rotations p
  in (q_new, h_new)

-- | Apply Givens rotations to zero out subdiagonal of Hessenberg matrix.
-- Returns list of (c, s, row_index) and the resulting R.
applyGivensQR :: forall n r e. (KnownNat n, M.Manifest r e, Floating e, Ord e)
              => Matrix n n r e -> Int -> ([(e, e, Int)], Matrix n n r e)
applyGivensQR h p = foldl step ([], h) [0..p-1]
  where
    nn = dimVal @n
    step (rots, hh) k =
      let (c, s) = givensRotation (hh ! (k, k)) (hh ! (k+1, k))
          hh' = applyGivensLeftSq c s k (k+1) hh
      in (rots ++ [(c, s, k)], hh')

-- | Apply Givens from left to a square matrix.
applyGivensLeftSq :: forall n r e. (KnownNat n, M.Manifest r e, Num e)
                  => e -> e -> Int -> Int -> Matrix n n r e -> Matrix n n r e
applyGivensLeftSq c s ri rk h =
  makeMatrix @n @n @r $ \i j ->
    if i == ri then
      c * (h ! (ri, j)) - s * (h ! (rk, j))
    else if i == rk then
      s * (h ! (ri, j)) + c * (h ! (rk, j))
    else
      h ! (i, j)

-- | Form RQ + σI from R and the Givens rotations.
formRQ :: forall n r e. (KnownNat n, M.Manifest r e, Num e)
       => Matrix n n r e -> [(e, e, Int)] -> e -> Int -> Matrix n n r e
formRQ r rots shift _ =
  let -- Apply rotations from the right: R·G₁ᵀ·G₂ᵀ·...
      rq = foldl (\mat (c, s, k) ->
        applyGivensRightSq c s k (k+1) mat
        ) r rots
  in -- Add back shift
    makeMatrix @(MatDim n) @(MatDim n) $ \i j ->
      if i == j then (rq ! (i, j)) + shift else rq ! (i, j)

type MatDim n = n  -- type alias to avoid ambiguity

-- | Apply Givens from right to a square matrix.
applyGivensRightSq :: forall n r e. (KnownNat n, M.Manifest r e, Num e)
                   => e -> e -> Int -> Int -> Matrix n n r e -> Matrix n n r e
applyGivensRightSq c s ci ck h =
  makeMatrix @n @n @r $ \i j ->
    if j == ci then
      c * (h ! (i, ci)) - s * (h ! (i, ck))
    else if j == ck then
      s * (h ! (i, ci)) + c * (h ! (i, ck))
    else
      h ! (i, j)

-- | Update Q by applying Givens rotations from the right.
updateQ :: forall n r e. (KnownNat n, M.Manifest r e, Num e)
        => Matrix n n r e -> [(e, e, Int)] -> Int -> Matrix n n r e
updateQ q rots _ = foldl (\qq (c, s, k) ->
  applyGivensRightSq c s k (k+1) qq
  ) q rots

-- | Extract eigenvalues from a (quasi-)upper triangular Schur form \(T\).
--
-- The Schur matrix \(T\) produced by 'schur' has \(1 \times 1\) and
-- \(2 \times 2\) diagonal blocks.  This function walks the diagonal and
-- extracts eigenvalues:
--
--   * A \(1 \times 1\) block \([t_{ii}]\) yields the real eigenvalue
--     \(\lambda = t_{ii}\).
--   * A \(2 \times 2\) block
--     \(\bigl[\begin{smallmatrix} a & b \\ c & d \end{smallmatrix}\bigr]\)
--     yields eigenvalues \(\tfrac{a + d}{2} \pm \sqrt{\tfrac{(a+d)^2}{4} -
--     (ad - bc)}\).  When the discriminant is negative (complex conjugate
--     pair) only the real part \(\tfrac{a + d}{2}\) is returned for each
--     eigenvalue, since this module operates over real scalars.
--
-- See GVL4 Section 7.5 for the definition of the real Schur form.
eigenvalues :: forall n r e. (KnownNat n, M.Manifest r e, Floating e, Ord e)
            => Matrix n n r e -> [e]
eigenvalues t =
  let nn = dimVal @n
  in go 0
  where
    nn = dimVal @n
    go i
      | i >= nn = []
      | i == nn - 1 = [t ! (i, i)]  -- Last 1×1 block
      | abs (t ! (i+1, i)) < 1e-12 * (abs (t ! (i, i)) + abs (t ! (i+1, i+1))) =
          -- 1×1 block
          t ! (i, i) : go (i + 1)
      | otherwise =
          -- 2×2 block: eigenvalues of [[a,b],[c,d]]
          let a = t ! (i, i)
              b = t ! (i, i+1)
              c = t ! (i+1, i)
              d = t ! (i+1, i+1)
              tr = a + d
              det_ = a * d - b * c
              disc = tr * tr / 4 - det_
          in if disc >= 0
             then (tr / 2 + sqrt disc) : (tr / 2 - sqrt disc) : go (i + 2)
             else tr / 2 : tr / 2 : go (i + 2)  -- Complex pair, return real parts
