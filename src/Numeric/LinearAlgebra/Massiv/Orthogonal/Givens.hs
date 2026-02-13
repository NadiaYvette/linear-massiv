{-# LANGUAGE AllowAmbiguousTypes #-}

-- |
-- Module      : Numeric.LinearAlgebra.Massiv.Orthogonal.Givens
-- Copyright   : (c) Nadia Chambers 2026
-- License     : BSD-3-Clause
-- Maintainer  : nadia.chambers@iohk.io
-- Stability   : experimental
--
-- Givens rotations for selective zeroing of matrix entries.
--
-- This module implements Givens (plane) rotations following Golub & Van
-- Loan, /Matrix Computations/, 4th edition (GVL4), Section 5.1.8,
-- pp. 240--243.
--
-- A Givens rotation is an orthogonal matrix that operates in a
-- two-dimensional subspace.  Given scalars \( a \) and \( b \), the
-- rotation matrix
--
-- \( G^T = \begin{bmatrix} c & -s \\ s & c \end{bmatrix} \)
--
-- is constructed so that
--
-- \( G^T \begin{bmatrix} a \\ b \end{bmatrix} = \begin{bmatrix} r \\ 0 \end{bmatrix} \)
--
-- where \( r = \sqrt{a^2 + b^2} \).  Our convention follows GVL4
-- Algorithm 5.1.3 (p. 240): \( c = a / r \), \( s = -b / r \).
--
-- Givens rotations are especially useful when only a small number of
-- sub-diagonal entries need to be zeroed (e.g., in Hessenberg or banded
-- matrices), whereas Householder reflections are preferred for zeroing
-- entire sub-columns at once.  Givens-based QR factorisation is the
-- method of choice for tridiagonal and Hessenberg eigenvalue problems
-- (GVL4 Section 5.2.8, p. 255).
--
-- __Complexity.__
--
-- * Computing a Givens rotation ('givensRotation'): \( O(1) \) flops
--   (one square root and a small number of divisions).
-- * Applying a Givens rotation to a row or column pair of an
--   \( m \times n \) matrix ('applyGivensLeft', 'applyGivensRight'):
--   \( O(n) \) or \( O(m) \) flops respectively (one pass over the
--   affected row or column pair).
module Numeric.LinearAlgebra.Massiv.Orthogonal.Givens
  ( -- * Givens rotation
    givensRotation
    -- * Apply Givens rotation
  , applyGivensLeft
  , applyGivensRight
  ) where

import qualified Data.Massiv.Array as M
import Data.Massiv.Array (Ix2(..), Sz(..))
import GHC.TypeNats (KnownNat)

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal

-- | Compute a Givens rotation (GVL4 Algorithm 5.1.3, p. 240).
--
-- Given scalars \( a \) and \( b \), compute cosine \( c \) and sine
-- \( s \) such that
--
-- \( \begin{bmatrix} c & s \\ -s & c \end{bmatrix}^T \begin{bmatrix} a \\ b \end{bmatrix} = \begin{bmatrix} r \\ 0 \end{bmatrix} \)
--
-- where \( r = \sqrt{a^2 + b^2} \).
--
-- The implementation avoids overflow and unnecessary computation by
-- distinguishing three cases:
--
-- * If \( b = 0 \), the rotation is the identity: \( c = 1, s = 0 \).
-- * If \( |b| > |a| \), the tangent \( \tau = -a/b \) is computed first,
--   then \( s = 1 / \sqrt{1 + \tau^2} \) and \( c = s \tau \).
-- * Otherwise, \( \tau = -b/a \), \( c = 1 / \sqrt{1 + \tau^2} \), and
--   \( s = c \tau \).
--
-- This avoids computing the potentially large quantity
-- \( r = \sqrt{a^2 + b^2} \) directly, which could overflow.
--
-- __Complexity:__ \( O(1) \) flops (one square root, a few multiplications
-- and divisions).
--
-- Returns @(c, s)@.
givensRotation :: (Floating e, Ord e) => e -> e -> (e, e)
givensRotation a b
  | b == 0    = (1, 0)
  | abs b > abs a =
      let tau = -a / b
          s = 1 / sqrt (1 + tau * tau)
          c = s * tau
      in (c, s)
  | otherwise =
      let tau = -b / a
          c = 1 / sqrt (1 + tau * tau)
          s = c * tau
      in (c, s)

-- | Apply a Givens rotation from the left to rows @i@ and @k@ of a matrix
-- (GVL4 Section 5.1.9, p. 241).
--
-- Performs the update
--
-- \( A([i,k], :) \leftarrow G^T \, A([i,k], :) \)
--
-- where \( G^T = \begin{bmatrix} c & -s \\ s & c \end{bmatrix} \).
-- Only rows @i@ and @k@ are modified; all other rows are untouched.
-- This is the standard operation used to zero out the \( (k, j) \)
-- entry of a matrix during Givens-based QR factorisation
-- (GVL4 Algorithm 5.2.3, p. 252).
--
-- __Complexity:__ \( O(n) \) flops, where \( n \) is the number of
-- columns.
applyGivensLeft :: forall m n r e. (KnownNat m, KnownNat n, M.Manifest r e, Num e)
                => e    -- ^ c
                -> e    -- ^ s
                -> Int  -- ^ row i
                -> Int  -- ^ row k
                -> Matrix m n r e -> Matrix m n r e
applyGivensLeft c s ri rk a =
  makeMatrix @m @n @r $ \i j ->
    if i == ri then
      c * (a ! (ri, j)) - s * (a ! (rk, j))
    else if i == rk then
      s * (a ! (ri, j)) + c * (a ! (rk, j))
    else
      a ! (i, j)

-- | Apply a Givens rotation from the right to columns @i@ and @k@ of a
-- matrix (GVL4 Section 5.1.9, p. 242).
--
-- Performs the update
--
-- \( A(:, [i,k]) \leftarrow A(:, [i,k]) \, G \)
--
-- where \( G = \begin{bmatrix} c & s \\ -s & c \end{bmatrix} \).
-- Only columns @i@ and @k@ are modified; all other columns are
-- untouched.  Right-multiplication by a Givens rotation is typically
-- used to accumulate the orthogonal factor \( Q \) during QR
-- factorisation (GVL4 Section 5.1.9, p. 242).
--
-- __Complexity:__ \( O(m) \) flops, where \( m \) is the number of
-- rows.
applyGivensRight :: forall m n r e. (KnownNat m, KnownNat n, M.Manifest r e, Num e)
                 => e    -- ^ c
                 -> e    -- ^ s
                 -> Int  -- ^ column i
                 -> Int  -- ^ column k
                 -> Matrix m n r e -> Matrix m n r e
applyGivensRight c s ci ck a =
  makeMatrix @m @n @r $ \i j ->
    if j == ci then
      c * (a ! (i, ci)) - s * (a ! (i, ck))
    else if j == ck then
      s * (a ! (i, ci)) + c * (a ! (i, ck))
    else
      a ! (i, j)
