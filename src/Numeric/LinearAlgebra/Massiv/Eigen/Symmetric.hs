{-# LANGUAGE AllowAmbiguousTypes #-}

-- |
-- Module      : Numeric.LinearAlgebra.Massiv.Eigen.Symmetric
-- Copyright   : (c) Nadia Chambers 2026
-- License     : BSD-3-Clause
-- Maintainer  : nadia.chambers@iohk.io
-- Stability   : experimental
--
-- Eigenvalue algorithms specialised to real symmetric matrices, following
-- Golub & Van Loan, /Matrix Computations/, 4th edition (GVL4), Chapter 8,
-- pp. 449--512.
--
-- Every real symmetric matrix \(A = A^T\) possesses a /spectral
-- decomposition/
--
-- \[
--   A = Q \, \Lambda \, Q^T
-- \]
--
-- where \(Q\) is orthogonal and \(\Lambda = \mathrm{diag}(\lambda_1, \ldots,
-- \lambda_n)\) contains the (real) eigenvalues.  This module provides three
-- approaches for computing this decomposition:
--
-- * __Householder tridiagonalisation__ (Algorithm 8.3.1, p. 459) reduces
--   \(A\) to a symmetric tridiagonal matrix \(T\) via orthogonal similarity:
--   \(A = Q_1 T Q_1^T\).
--
-- * __Implicit symmetric QR algorithm with Wilkinson shift__
--   (Algorithm 8.3.3, p. 462) iterates on \(T\) to produce the diagonal
--   \(\Lambda\).  The /Wilkinson shift/ is chosen as the eigenvalue of the
--   trailing \(2 \times 2\) submatrix closest to \(a_{nn}\), guaranteeing
--   global (and in practice rapid) convergence.
--
-- * __Classical Jacobi method__ (Section 8.5, pp. 488--498) diagonalises
--   \(A\) directly via a sequence of plane rotations.  Although slower than
--   the QR-based approach for large matrices, the Jacobi method is highly
--   parallelisable and inherently accurate.  With /cyclic/ ordering of
--   rotations the convergence is /quadratic/ (GVL4, p. 495).
module Numeric.LinearAlgebra.Massiv.Eigen.Symmetric
  ( -- * Tridiagonal reduction (Algorithm 8.3.1)
    tridiagonalize
    -- * Symmetric QR (Algorithm 8.3.3)
  , symmetricEigen
    -- * Jacobi method (Section 8.5)
  , jacobiEigen
  ) where

import qualified Data.Massiv.Array as M
import Data.Massiv.Array (Ix2(..), Sz(..))
import GHC.TypeNats (KnownNat)

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.Orthogonal.Givens (givensRotation)

-- | Reduce a symmetric matrix to tridiagonal form via Householder similarity
-- transformations (GVL4 Algorithm 8.3.1, p. 459).
--
-- Given a symmetric \(A \in \mathbb{R}^{n \times n}\), computes an orthogonal
-- \(Q\) (product of \(n - 2\) Householder reflectors) and a symmetric
-- tridiagonal \(T\) such that
--
-- \[
--   A = Q \, T \, Q^T
-- \]
--
-- The tridiagonal matrix is returned in compact form as a pair of vectors:
-- the /diagonal/ \((t_{11}, \ldots, t_{nn})\) and the /subdiagonal/
-- \((t_{21}, \ldots, t_{n,n-1})\).  Only the lower triangle of \(A\) is
-- accessed.
--
-- __Complexity:__ \(\frac{4}{3} n^3\) flops (GVL4, p. 460).
--
-- Returns @(Q, diagonal, subdiagonal)@.
tridiagonalize :: forall n r e. (KnownNat n, M.Manifest r e, Floating e, Ord e)
               => Matrix n n r e
               -> (Matrix n n r e, Vector n r e, Vector n r e)
tridiagonalize a =
  let nn = dimVal @n
      go :: Int -> Matrix n n r e -> Matrix n n r e
          -> (Matrix n n r e, Matrix n n r e)
      go k q_ t_
        | k >= nn - 2 = (q_, t_)
        | otherwise =
          -- Householder to zero out t(k+2:n, k)
          let x0 = t_ ! (k+1, k)
              sigma = foldl' (\acc i -> acc + (t_ ! (i, k)) * (t_ ! (i, k))) 0 [k+2..nn-1]
          in if sigma == 0
             then go (k + 1) q_ t_
             else
              let mu = sqrt (x0 * x0 + sigma)
                  v0 = if x0 <= 0 then x0 - mu else -sigma / (x0 + mu)
                  beta = 2 * v0 * v0 / (sigma + v0 * v0)
                  v = makeVector @n @r $ \i ->
                    if i <= k then 0
                    else if i == k + 1 then 1
                    else (t_ ! (i, k)) / v0
                  -- Similarity transform: T ← (I - βvvᵀ)T(I - βvvᵀ)
                  t1 = applyLeft v beta t_
                  t2 = applyRight t1 v beta
                  q_new = applyRight q_ v beta
              in go (k + 1) q_new t2

      q0 = identityMatrix @n @r
      (qFinal, tFinal) = go 0 q0 a

      -- Extract diagonal and subdiagonal
      diag_ = makeVector @n @r $ \i -> tFinal ! (i, i)
      subdiag = makeVector @n @r $ \i ->
        if i < nn - 1 then tFinal ! (i + 1, i) else 0
  in (qFinal, diag_, subdiag)
  where
    nn = dimVal @n
    applyLeft v beta h =
      makeMatrix @n @n @r $ \i j ->
        let wj = beta * foldl' (\acc k -> acc + (v !. k) * (h ! (k, j))) 0 [0..nn-1]
        in (h ! (i, j)) - (v !. i) * wj
    applyRight h v beta =
      makeMatrix @n @n @r $ \i j ->
        let wi = beta * foldl' (\acc k -> acc + (h ! (i, k)) * (v !. k)) 0 [0..nn-1]
        in (h ! (i, j)) - wi * (v !. j)

-- | Full symmetric eigenvalue decomposition via the implicit symmetric QR
-- algorithm with Wilkinson shift (GVL4 Algorithm 8.3.3, p. 462).
--
-- Given a symmetric matrix \(A = A^T\), computes the spectral decomposition
--
-- \[
--   A = Q \, \Lambda \, Q^T, \qquad
--   \Lambda = \mathrm{diag}(\lambda_1, \ldots, \lambda_n)
-- \]
--
-- The algorithm proceeds in two phases:
--
--   1. Reduce \(A\) to symmetric tridiagonal form \(T\) via
--      'tridiagonalize'.
--   2. Apply the implicit symmetric QR iteration with /Wilkinson shift/ on
--      \(T\).  At each step the shift is chosen as the eigenvalue of the
--      trailing \(2 \times 2\) block
--      \(\bigl[\begin{smallmatrix} t_{p-1,p-1} & t_{p-1,p} \\ t_{p,p-1} &
--      t_{pp} \end{smallmatrix}\bigr]\) that is closest to \(t_{pp}\) (GVL4,
--      p. 462).  Converged eigenvalues are deflated from the bottom of the
--      active window.
--
-- Returns @(eigenvalues, Q)@ where @eigenvalues@ is a vector of eigenvalues
-- and the columns of @Q@ are the corresponding orthonormal eigenvectors.
symmetricEigen :: forall n r e. (KnownNat n, M.Manifest r e, Floating e, Ord e)
               => Matrix n n r e
               -> Int   -- ^ Max iterations
               -> e     -- ^ Tolerance
               -> (Vector n r e, Matrix n n r e)
symmetricEigen a maxIter tol =
  let nn = dimVal @n
      (q0, diag_, subdiag) = tridiagonalize a
      -- Apply QR iteration on the tridiagonal matrix
      (eigvals, qFinal) = tridiagQR nn q0 diag_ subdiag maxIter tol
  in (eigvals, qFinal)

-- | QR iteration on tridiagonal matrix.
tridiagQR :: forall n r e. (KnownNat n, M.Manifest r e, Floating e, Ord e)
          => Int
          -> Matrix n n r e     -- ^ Accumulated Q
          -> Vector n r e       -- ^ Diagonal
          -> Vector n r e       -- ^ Subdiagonal
          -> Int -> e
          -> (Vector n r e, Matrix n n r e)
tridiagQR nn q diag_ subdiag maxIter tol = go 0 q diag_ subdiag (nn - 1)
  where
    go :: Int -> Matrix n n r e -> Vector n r e -> Vector n r e -> Int
       -> (Vector n r e, Matrix n n r e)
    go iter q_ d sd p
      | iter >= maxIter = (d, q_)
      | p <= 0 = (d, q_)
      | otherwise =
        let sdp = abs (sd !. (p - 1))
            diagSum = abs (d !. (p - 1)) + abs (d !. p)
        in if sdp <= tol * diagSum
           then -- Converged: set subdiagonal to zero
             let sd' = makeVector @n @r $ \i ->
                   if i == p - 1 then 0 else sd !. i
             in go iter q_ d sd' (p - 1)
           else
             -- Wilkinson shift from trailing 2×2
             let dp1 = d !. (p - 1)
                 dp  = d !. p
                 sp1 = sd !. (p - 1)
                 delta = (dp1 - dp) / 2
                 sgn = if delta >= 0 then 1 else -1
                 shift = dp - sp1 * sp1 / (delta + sgn * sqrt (delta * delta + sp1 * sp1))
                 -- Apply implicit QR step via Givens rotations
                 (d', sd', q_new) = implicitQRStep d sd q_ shift p
             in go (iter + 1) q_new d' sd' p

-- | One implicit symmetric QR step via Givens rotations.
implicitQRStep :: forall n r e. (KnownNat n, M.Manifest r e, Floating e, Ord e)
               => Vector n r e -> Vector n r e -> Matrix n n r e -> e -> Int
               -> (Vector n r e, Vector n r e, Matrix n n r e)
implicitQRStep d sd q shift p =
  let nn = dimVal @n
      -- Build the tridiagonal as a full matrix, do the step, extract back
      t = makeMatrix @n @n @r $ \i j ->
        if i == j then d !. i
        else if i == j + 1 && j < nn - 1 then sd !. j
        else if j == i + 1 && i < nn - 1 then sd !. i
        else 0
      -- Shift
      t_shifted = makeMatrix @n @n @r $ \i j ->
        if i == j then (t ! (i, j)) - shift else t ! (i, j)
      -- QR step via Givens
      (rots, r) = foldl (\(rs, hh) k ->
        if k >= p then (rs, hh)
        else
          let (c, s) = givensRotation (hh ! (k, k)) (hh ! (k+1, k))
              hh' = applyGivensL c s k (k+1) hh
          in (rs ++ [(c, s, k)], hh')
        ) ([], t_shifted) [0..p-1]
      -- Form RQ + σI
      rq = foldl (\mat (c, s, k) -> applyGivensR c s k (k+1) mat) r rots
      t_new = makeMatrix @n @n @r $ \i j ->
        if i == j then (rq ! (i, j)) + shift else rq ! (i, j)
      -- Update Q
      q_new = foldl (\qq (c, s, k) -> applyGivensR c s k (k+1) qq) q rots
      -- Extract diagonal and subdiagonal
      d' = makeVector @n @r $ \i -> t_new ! (i, i)
      sd' = makeVector @n @r $ \i ->
        if i < nn - 1 then t_new ! (i+1, i) else 0
  in (d', sd', q_new)
  where
    nn = dimVal @n
    applyGivensL c s ri rk h =
      makeMatrix @n @n @r $ \i j ->
        if i == ri then c * (h ! (ri, j)) - s * (h ! (rk, j))
        else if i == rk then s * (h ! (ri, j)) + c * (h ! (rk, j))
        else h ! (i, j)
    applyGivensR c s ci ck h =
      makeMatrix @n @n @r $ \i j ->
        if j == ci then c * (h ! (i, ci)) - s * (h ! (i, ck))
        else if j == ck then s * (h ! (i, ci)) + c * (h ! (i, ck))
        else h ! (i, j)

-- | Classical Jacobi eigenvalue method (GVL4 Section 8.5, pp. 488--498).
--
-- Diagonalises a symmetric matrix \(A\) by iteratively applying plane
-- (Jacobi) rotations \(J(p, q, \theta)\) chosen to annihilate a selected
-- off-diagonal pair \((a_{pq}, a_{qp})\):
--
-- \[
--   A_{k+1} = J^T A_k \, J
-- \]
--
-- A full /sweep/ visits every off-diagonal pair \((p, q)\) with
-- \(p < q\) in cyclic (row-by-row) order.  With this ordering the
-- convergence is /quadratic/: the off-diagonal Frobenius norm
-- \(\mathrm{off}(A)\) satisfies
-- \(\mathrm{off}(A_{k+1}) \leq c \, \mathrm{off}(A_k)^2\) (GVL4, p. 495).
--
-- The Jacobi method is unconditionally convergent for every symmetric matrix
-- and is well-suited to parallel and high-accuracy settings, although it is
-- generally slower than the QR-based algorithm ('symmetricEigen') for large
-- matrices.
--
-- Returns @(eigenvalues, Q)@ where @Q@ accumulates all applied rotations so
-- that \(A = Q \Lambda Q^T\).
jacobiEigen :: forall n r e. (KnownNat n, M.Manifest r e, Floating e, Ord e)
            => Matrix n n r e
            -> Int   -- ^ Max sweeps (one sweep = all off-diagonal pairs)
            -> e     -- ^ Tolerance on off-diagonal Frobenius norm
            -> (Vector n r e, Matrix n n r e)
jacobiEigen a maxSweeps tol =
  let nn = dimVal @n
      q0 = identityMatrix @n @r
      (aFinal, qFinal) = go 0 a q0
      eigvals = makeVector @n @r $ \i -> aFinal ! (i, i)
  in (eigvals, qFinal)
  where
    nn = dimVal @n
    go :: Int -> Matrix n n r e -> Matrix n n r e -> (Matrix n n r e, Matrix n n r e)
    go sweep a_ q_
      | sweep >= maxSweeps = (a_, q_)
      | offDiagNorm a_ < tol = (a_, q_)
      | otherwise =
        -- One sweep: iterate over all off-diagonal pairs
        let (a_new, q_new) = foldl (\(aa, qq) (p, q) ->
              if abs (aa ! (p, q)) < tol * 1e-3
              then (aa, qq)
              else
                let (c, s) = jacobiRotation (aa ! (p, p)) (aa ! (p, q)) (aa ! (q, q))
                    aa' = applyJacobiSimilarity c s p q aa
                    qq' = applyJacobiRight c s p q qq
                in (aa', qq')
              ) (a_, q_) [(i, j) | i <- [0..nn-2], j <- [i+1..nn-1]]
        in go (sweep + 1) a_new q_new

    -- Off-diagonal Frobenius norm
    offDiagNorm :: Matrix n n r e -> e
    offDiagNorm m = sqrt $ foldl' (\acc (i, j) ->
      acc + (m ! (i, j)) * (m ! (i, j))
      ) 0 [(i, j) | i <- [0..nn-1], j <- [0..nn-1], i /= j]

    -- Compute Jacobi rotation angles
    jacobiRotation :: e -> e -> e -> (e, e)
    jacobiRotation app apq aqq
      | apq == 0 = (1, 0)
      | otherwise =
        let tau = (aqq - app) / (2 * apq)
            t = if tau >= 0
                then 1 / (tau + sqrt (1 + tau * tau))
                else 1 / (tau - sqrt (1 + tau * tau))
            c = 1 / sqrt (1 + t * t)
            s = t * c
        in (c, s)

    -- Apply Jacobi rotation as similarity: A ← JᵀAJ
    applyJacobiSimilarity :: e -> e -> Int -> Int -> Matrix n n r e -> Matrix n n r e
    applyJacobiSimilarity c s p q a_ =
      makeMatrix @n @n @r $ \i j ->
        let aij = a_ ! (i, j)
        in if i == p && j == p then
             c*c*(a_ ! (p,p)) - 2*s*c*(a_ ! (p,q)) + s*s*(a_ ! (q,q))
           else if i == q && j == q then
             s*s*(a_ ! (p,p)) + 2*s*c*(a_ ! (p,q)) + c*c*(a_ ! (q,q))
           else if i == p && j == q then 0
           else if i == q && j == p then 0
           else if i == p then c*(a_ ! (p,j)) - s*(a_ ! (q,j))
           else if i == q then s*(a_ ! (p,j)) + c*(a_ ! (q,j))
           else if j == p then c*(a_ ! (i,p)) - s*(a_ ! (i,q))
           else if j == q then s*(a_ ! (i,p)) + c*(a_ ! (i,q))
           else aij

    -- Apply Jacobi rotation from right: Q ← Q·J
    applyJacobiRight :: e -> e -> Int -> Int -> Matrix n n r e -> Matrix n n r e
    applyJacobiRight c s p q qq =
      makeMatrix @n @n @r $ \i j ->
        if j == p then c*(qq ! (i,p)) - s*(qq ! (i,q))
        else if j == q then s*(qq ! (i,p)) + c*(qq ! (i,q))
        else qq ! (i, j)
