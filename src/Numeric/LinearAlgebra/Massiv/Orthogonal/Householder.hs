{-# LANGUAGE AllowAmbiguousTypes #-}

-- |
-- Module      : Numeric.LinearAlgebra.Massiv.Orthogonal.Householder
-- Copyright   : (c) Nadia Chambers 2026
-- License     : BSD-3-Clause
-- Maintainer  : nadia.chambers@iohk.io
-- Stability   : experimental
--
-- Householder reflections for orthogonal triangularisation.
--
-- This module implements the Householder reflection (also known as a
-- Householder transformation), following Golub & Van Loan, /Matrix
-- Computations/, 4th edition (GVL4), Section 5.1, pp. 236--243.
--
-- As GVL4 states (p. 236): \"The Householder reflection is the most
-- important tool in matrix computations.\"  A Householder reflector is a
-- matrix of the form
--
-- \( P = I - \beta v v^T \)
--
-- where \( v \) is the /Householder vector/ and \( \beta = 2 / (v^T v) \).
-- The key property of \( P \) is that it is both symmetric and orthogonal:
--
-- \( P = P^T = P^{-1} \)
--
-- Given an input vector \( x \), the Householder vector \( v \) and scalar
-- \( \beta \) are chosen so that
--
-- \( P x = (I - \beta v v^T) x = \| x \|_2 \, e_1 \)
--
-- where \( e_1 \) is the first standard basis vector.  This is the
-- fundamental operation behind Householder QR factorisation (GVL4
-- Algorithm 5.2.1) and many other matrix decompositions.
--
-- __Complexity.__
--
-- * Computing the Householder vector ('householderVector'): \( O(n) \) flops.
-- * Applying a Householder reflection to an \( m \times n \) matrix
--   ('applyHouseholderLeft', 'applyHouseholderRight'): \( O(mn) \) flops.
-- * Forming the explicit reflector matrix ('householderMatrix'): \( O(n^2) \)
--   flops, but this should be avoided in favour of implicit application
--   whenever possible.
module Numeric.LinearAlgebra.Massiv.Orthogonal.Householder
  ( -- * Householder vector
    householderVector
    -- * Apply Householder reflection
  , applyHouseholderLeft
  , applyHouseholderRight
    -- * Construct explicit reflector
  , householderMatrix
  ) where

import qualified Data.Massiv.Array as M
import Data.Massiv.Array (Ix1, Ix2(..), Sz(..))
import GHC.TypeNats (KnownNat)

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal

-- | Compute a Householder vector (GVL4 Algorithm 5.1.1, p. 236).
--
-- Given a vector \( x \in \mathbb{R}^n \), compute the Householder vector
-- \( v \) and scalar \( \beta \) such that
--
-- \( (I - \beta \, v \, v^T) \, x = \| x \|_2 \, e_1 \)
--
-- where \( e_1 \) is the first standard basis vector.  By convention the
-- first component of \( v \) is normalised to \( v_1 = 1 \), which allows
-- it to be stored implicitly in the sub-diagonal part of a matrix during
-- QR factorisation.
--
-- The implementation follows GVL4 Algorithm 5.1.1 exactly, including the
-- careful treatment of the sign of \( x_1 \) to avoid catastrophic
-- cancellation.  When \( x \) is already a non-negative multiple of
-- \( e_1 \), the function returns \( \beta = 0 \) (i.e., the identity
-- transformation).
--
-- __Complexity:__ \( O(n) \) flops.
--
-- Returns @(v, beta)@.
householderVector :: forall n r e. (KnownNat n, M.Manifest r e, Floating e, Ord e)
                  => Vector n r e -> (Vector n r e, e)
householderVector x =
  let nn = dimVal @n
      x0 = x !. 0
      -- σ = x(2:n)ᵀ · x(2:n)
      sigma = foldl' (\acc i -> acc + (x !. i) * (x !. i)) 0 [1..nn-1]
  in if sigma == 0 && x0 >= 0
    then -- x is already a positive multiple of e1
      ( makeVector @n @r $ \i -> if i == 0 then 1 else 0
      , 0
      )
    else if sigma == 0
    then -- x = -α·e₁
      ( makeVector @n @r $ \i -> if i == 0 then 1 else 0
      , 2
      )
    else
      let mu = sqrt (x0 * x0 + sigma)
          v0 = if x0 <= 0 then x0 - mu else -sigma / (x0 + mu)
          beta = 2 * v0 * v0 / (sigma + v0 * v0)
          v = makeVector @n @r $ \i ->
            if i == 0 then 1 else (x !. i) / v0
      in (v, beta)

-- | Apply a Householder reflection from the left (GVL4 Section 5.1, p. 236).
--
-- Given a Householder vector \( v \in \mathbb{R}^m \), scalar \( \beta \),
-- and matrix \( A \in \mathbb{R}^{m \times n} \), compute
--
-- \( A \leftarrow (I - \beta \, v \, v^T) \, A = A - \beta \, v \, (A^T v)^T \)
--
-- The computation is performed without forming \( P \) explicitly.
-- Instead, the intermediate vector \( w = \beta \, A^T v \) is computed
-- first, and then the rank-1 update \( A \leftarrow A - v \, w^T \) is
-- applied.  This is the standard technique described in GVL4
-- Section 5.1.4 (p. 238).
--
-- __Complexity:__ \( O(mn) \) flops.
applyHouseholderLeft :: forall m n r e. (KnownNat m, KnownNat n, M.Manifest r e, Num e)
                     => Vector m r e -> e -> Matrix m n r e -> Matrix m n r e
applyHouseholderLeft v beta a =
  let mm = dimVal @m
      c  = dimVal @n
  in makeMatrix @m @n @r $ \i j ->
    let -- w = βAᵀv, w(j) = β · Σᵢ v(i)·A(i,j)
        wj = beta * foldl' (\acc k -> acc + (v !. k) * (a ! (k, j))) 0 [0..mm-1]
    in (a ! (i, j)) - (v !. i) * wj

-- | Apply a Householder reflection from the right (GVL4 Section 5.1, p. 236).
--
-- Given a matrix \( A \in \mathbb{R}^{m \times n} \), Householder vector
-- \( v \in \mathbb{R}^n \), and scalar \( \beta \), compute
--
-- \( A \leftarrow A \, (I - \beta \, v \, v^T) = A - \beta \, (A \, v) \, v^T \)
--
-- As with 'applyHouseholderLeft', the reflector is never formed
-- explicitly.  The intermediate vector \( w = \beta \, A \, v \) is
-- computed first, followed by the rank-1 update
-- \( A \leftarrow A - w \, v^T \).  See GVL4 Section 5.1.4 (p. 238).
--
-- __Complexity:__ \( O(mn) \) flops.
applyHouseholderRight :: forall m n r e. (KnownNat m, KnownNat n, M.Manifest r e, Num e)
                      => Matrix m n r e -> Vector n r e -> e -> Matrix m n r e
applyHouseholderRight a v beta =
  let c = dimVal @n
  in makeMatrix @m @n @r $ \i j ->
    let -- w = β·A·v, w(i) = β · Σⱼ A(i,j)·v(j)
        wi = beta * foldl' (\acc k -> acc + (a ! (i, k)) * (v !. k)) 0 [0..c-1]
    in (a ! (i, j)) - wi * (v !. j)

-- | Construct the explicit Householder reflector matrix (GVL4 Section 5.1, p. 236).
--
-- Given a Householder vector \( v \in \mathbb{R}^n \) and scalar
-- \( \beta \), form the \( n \times n \) matrix
--
-- \( H = I - \beta \, v \, v^T \)
--
-- The resulting matrix is both symmetric and orthogonal:
-- \( H = H^T = H^{-1} \).
--
-- __Note:__ In most numerical algorithms it is preferable to apply the
-- Householder transformation implicitly via 'applyHouseholderLeft' or
-- 'applyHouseholderRight' rather than forming \( H \) explicitly.
-- Forming the explicit matrix costs \( O(n^2) \) flops and storage, and
-- subsequent multiplication with it costs \( O(n^3) \) rather than the
-- \( O(mn) \) achievable by implicit application.
--
-- __Complexity:__ \( O(n^2) \) flops.
householderMatrix :: forall n r e. (KnownNat n, M.Manifest r e, Num e)
                  => Vector n r e -> e -> Matrix n n r e
householderMatrix v beta =
  makeMatrix @n @n @r $ \i j ->
    let ident = if i == j then 1 else 0
    in ident - beta * (v !. i) * (v !. j)
