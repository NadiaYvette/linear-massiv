{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE MagicHash #-}
{-# LANGUAGE UnboxedTuples #-}
{-# LANGUAGE BangPatterns #-}

-- |
-- Module      : Numeric.LinearAlgebra.Massiv.Solve.Cholesky
-- Copyright   : (c) Nadia Chambers 2026
-- License     : BSD-3-Clause
-- Maintainer  : nadia.chambers@iohk.io
-- Stability   : experimental
--
-- = Cholesky Factorization
--
-- Cholesky decomposition for symmetric positive definite (SPD) matrices,
-- following Golub & Van Loan, /Matrix Computations/, 4th edition (GVL4),
-- Section 4.2, pp. 163--169.
--
-- For any SPD matrix \(A \in \mathbb{R}^{n \times n}\), there exists a
-- unique lower triangular matrix \(G\) with positive diagonal entries such
-- that
--
-- \[
--   A = G G^T
-- \]
--
-- (GVL4 Theorem 4.2.1, p. 163). This is the /Cholesky factorization/.
-- Because it exploits symmetry, the Cholesky factorization requires only
-- half the work of a general LU factorization: \(O(n^3/3)\) flops vs.
-- \(O(2n^3/3)\) flops (GVL4 p. 165).
--
-- +-------------------+-------------------------------+---------------------------------+
-- | Function          | Algorithm                     | Reference                       |
-- +===================+===============================+=================================+
-- | 'cholesky'        | Outer-product Cholesky        | GVL4 Algorithm 4.2.1, p. 164    |
-- +-------------------+-------------------------------+---------------------------------+
-- | 'choleskyGaxpy'   | Gaxpy (column-oriented)       | GVL4 Algorithm 4.2.2, p. 165    |
-- +-------------------+-------------------------------+---------------------------------+
-- | 'choleskySolve'   | Solve via \(A = GG^T\)        | GVL4 Section 4.2, p. 166        |
-- +-------------------+-------------------------------+---------------------------------+
--
-- == Complexity
--
-- The factorization costs \(O(n^3/3)\) flops -- exactly half of LU
-- (GVL4 p. 165). The subsequent pair of triangular solves adds \(O(n^2)\)
-- flops.
--
-- == Type Safety
--
-- Matrix dimensions are tracked at the type level via 'KnownNat', so the
-- compiler statically ensures the coefficient matrix is square and the
-- right-hand side vector has a conforming length. Note that positive
-- definiteness is /not/ checked at the type level; if the input matrix is
-- not SPD, the algorithm may produce NaN values from taking the square root
-- of a negative number.
module Numeric.LinearAlgebra.Massiv.Solve.Cholesky
  ( -- * Cholesky factorization (\(A = GG^T\))
    cholesky
  , choleskyGaxpy
    -- * Solving with Cholesky (\(Ax = b\))
  , choleskySolve
  , choleskySolveP
  ) where

import qualified Data.Massiv.Array as M
import Data.Massiv.Array (Ix1, Ix2(..), Sz(..), unwrapByteArray, unwrapByteArrayOffset,
                          unwrapMutableByteArray, unwrapMutableByteArrayOffset)
import GHC.TypeNats (KnownNat)
import GHC.Exts
import GHC.ST (ST(..))
import Data.Primitive.ByteArray (ByteArray(..), MutableByteArray(..))

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.Solve.Triangular (forwardSub, backSub)
import Numeric.LinearAlgebra.Massiv.BLAS.Level3 (transpose)
import Numeric.LinearAlgebra.Massiv.Internal.Kernel
  (rawCholColumn, rawForwardSubCholPacked, rawBackSubCholTPacked)

-- | Outer-product Cholesky factorization (GVL4 Algorithm 4.2.1, p. 164).
--
-- Given a symmetric positive definite \(n \times n\) matrix \(A\), computes
-- the unique lower triangular matrix \(G\) with positive diagonal entries
-- such that
--
-- \[
--   A = G G^T
-- \]
--
-- Only the lower triangle of \(A\) is accessed; the upper triangle is
-- ignored.
--
-- The algorithm processes one column at a time using an /outer-product/
-- update of the trailing submatrix.
--
-- ==== Type-safety guarantees
--
-- 'KnownNat' \(n\) statically ensures \(A\) is square. The 'Floating'
-- constraint provides 'sqrt'. Positive definiteness is a run-time
-- precondition; violation may produce NaN from \(\sqrt{g_{jj}}\) when
-- \(g_{jj} < 0\).
--
-- ==== Complexity
--
-- \(O(n^3/3)\) flops (GVL4 p. 165).
--
-- ==== Reference
--
-- Golub & Van Loan, /Matrix Computations/, 4th ed., Algorithm 4.2.1
-- (Outer Product Cholesky), p. 164. Existence: Theorem 4.2.1, p. 163.
cholesky :: forall n r e. (KnownNat n, M.Manifest r e, Floating e, Ord e)
         => Matrix n n r e -> Matrix n n r e
cholesky (MkMatrix a) =
  let nn = dimVal @n
  in MkMatrix $ M.createArrayST_ (M.Sz2 nn nn) $ \mg -> do
    -- Initialize to zero
    mapM_ (\i -> mapM_ (\j -> M.write_ mg (i :. j) 0) [0..nn-1]) [0..nn-1]
    -- Copy lower triangle of A into working storage
    mapM_ (\j -> mapM_ (\i -> do
      let aij = M.index' a (i :. j)
      M.write_ mg (i :. j) aij
      ) [j..nn-1]) [0..nn-1]

    -- Outer product Cholesky
    mapM_ (\j -> do
      -- Subtract contributions from previous columns
      mapM_ (\k -> do
        gjk <- M.readM mg (j :. k)
        mapM_ (\i -> do
          gij <- M.readM mg (i :. j)
          gik <- M.readM mg (i :. k)
          M.write_ mg (i :. j) (gij - gik * gjk)
          ) [j..nn-1]
        ) [0..j-1]

      -- Scale column
      gjj <- M.readM mg (j :. j)
      let sjj = sqrt gjj
      M.write_ mg (j :. j) sjj
      mapM_ (\i -> do
        gij <- M.readM mg (i :. j)
        M.write_ mg (i :. j) (gij / sjj)
        ) [j+1..nn-1]
      ) [0..nn-1]

-- | Gaxpy (column-oriented) Cholesky factorization (GVL4 Algorithm 4.2.2,
-- p. 165).
--
-- Functionally equivalent to 'cholesky', but uses a /gaxpy/ (generalised
-- @y <- y - Gx@) inner loop that accumulates updates into each column
-- before normalising. This access pattern is advantageous for column-major
-- storage because it streams through contiguous memory.
--
-- ==== Mathematical definition
--
-- Computes the same \(G\) such that \(A = GG^T\) as 'cholesky'.
--
-- ==== Type-safety guarantees
--
-- Identical to 'cholesky'.
--
-- ==== Complexity
--
-- \(O(n^3/3)\) flops (GVL4 p. 165).
--
-- ==== Reference
--
-- Golub & Van Loan, /Matrix Computations/, 4th ed., Algorithm 4.2.2
-- (Gaxpy Cholesky), p. 165.
choleskyGaxpy :: forall n r e. (KnownNat n, M.Manifest r e, Floating e, Ord e)
              => Matrix n n r e -> Matrix n n r e
choleskyGaxpy (MkMatrix a) =
  let nn = dimVal @n
  in MkMatrix $ M.createArrayST_ (M.Sz2 nn nn) $ \mg -> do
    -- Initialize: copy lower triangle of A
    mapM_ (\i -> mapM_ (\j ->
      if i >= j
        then M.write_ mg (i :. j) (M.index' a (i :. j))
        else M.write_ mg (i :. j) 0
      ) [0..nn-1]) [0..nn-1]

    -- Column by column
    mapM_ (\j -> do
      -- Update column j using previous columns (gaxpy)
      mapM_ (\k -> do
        gjk <- M.readM mg (j :. k)
        mapM_ (\i -> do
          gij <- M.readM mg (i :. j)
          gik <- M.readM mg (i :. k)
          M.write_ mg (i :. j) (gij - gik * gjk)
          ) [j..nn-1]
        ) [0..j-1]

      -- Scale by 1/sqrt(g(j,j))
      gjj <- M.readM mg (j :. j)
      let sjj = sqrt gjj
      mapM_ (\i -> do
        gij <- M.readM mg (i :. j)
        M.write_ mg (i :. j) (gij / sjj)
        ) [j..nn-1]
      ) [0..nn-1]

-- | Solve \(Ax = b\) where \(A\) is symmetric positive definite, using
-- Cholesky factorization (GVL4 Section 4.2, p. 166).
--
-- The algorithm proceeds in three stages:
--
-- 1. Factor \(A = GG^T\) via 'cholesky' (Algorithm 4.2.1).
-- 2. Solve \(Gy = b\) by forward substitution ('forwardSub').
-- 3. Solve \(G^T x = y\) by back substitution ('backSub').
--
-- ==== Type-safety guarantees
--
-- 'KnownNat' \(n\) enforces that \(A\) is \(n \times n\) and \(b\) has
-- length \(n\) at compile time.
--
-- ==== Complexity
--
-- \(O(n^3/3)\) flops for the factorization plus \(O(n^2)\) flops for the
-- two triangular solves (GVL4 p. 166).
--
-- ==== Reference
--
-- Golub & Van Loan, /Matrix Computations/, 4th ed., Section 4.2,
-- pp. 163--169.
choleskySolve :: forall n r e. (KnownNat n, M.Manifest r e, Floating e, Ord e)
              => Matrix n n r e -> Vector n r e -> Vector n r e
choleskySolve a b =
  let g = cholesky a
      gt = transpose g
      y = forwardSub g b
  in backSub gt y

-- | Specialised Cholesky solve for @P Double@.
-- Does Cholesky factorisation + solve entirely using raw ByteArray# primops,
-- avoiding separate G/G^T matrix construction.
choleskySolveP :: forall n. KnownNat n
               => Matrix n n M.P Double -> Vector n M.P Double -> Vector n M.P Double
choleskySolveP (MkMatrix a) (MkVector b) =
  let nn = dimVal @n
  in createVector @n @M.P $ \mx -> do
    -- Allocate n×n working storage for G, copy lower triangle of A
    mg <- M.newMArray @M.P (M.Sz2 nn nn) (0 :: Double)
    let mbaG = unwrapMutableByteArray mg
        offG = unwrapMutableByteArrayOffset mg
    -- Copy lower triangle using raw primops
    copyLowerTriangle a mbaG offG nn

    -- Phase 1: Cholesky factorisation column by column
    mapM_ (rawCholColumn mbaG offG nn) [0..nn-1]

    -- Phase 2: Freeze G and prepare RHS
    frozenG <- M.freezeS mg
    let baG = unwrapByteArray frozenG
        offFG = unwrapByteArrayOffset frozenG

    -- Copy b into output vector
    let mbaX = unwrapMutableByteArray mx
        offX = unwrapMutableByteArrayOffset mx
    copyVectorRaw b mbaX offX nn

    -- Phase 3: Forward substitution (Gy = b)
    rawForwardSubCholPacked baG offFG nn mbaX offX

    -- Phase 4: Back substitution (G^T x = y), without forming G^T
    rawBackSubCholTPacked baG offFG nn mbaX offX
{-# NOINLINE choleskySolveP #-}

-- | Copy lower triangle of an immutable 2D P array into a mutable byte array.
copyLowerTriangle :: M.Array M.P Ix2 Double -> MutableByteArray s -> Int -> Int -> ST s ()
copyLowerTriangle src (MutableByteArray mba_dst) (I# off_dst) (I# n) = ST $ \s0 ->
  let ba_src = case unwrapByteArray src of ByteArray ba -> ba
      off_src = case unwrapByteArrayOffset src of I# o -> o
      goI i s
        | isTrue# (i >=# n) = s
        | otherwise = goI (i +# 1#) (goJ i 0# s)
      goJ i j s
        | isTrue# (j ># i) = s
        | otherwise =
            let v = indexDoubleArray# ba_src (off_src +# i *# n +# j)
            in case writeDoubleArray# mba_dst (off_dst +# i *# n +# j) v s of
                 s' -> goJ i (j +# 1#) s'
  in (# goI 0# s0, () #)
{-# INLINE copyLowerTriangle #-}

-- | Copy an immutable P vector into a mutable byte array.
copyVectorRaw :: M.Array M.P Ix1 Double -> MutableByteArray s -> Int -> Int -> ST s ()
copyVectorRaw src (MutableByteArray mba_dst) (I# off_dst) (I# n) = ST $ \s0 ->
  let ba_src = case unwrapByteArray src of ByteArray ba -> ba
      off_src = case unwrapByteArrayOffset src of I# o -> o
      go i s
        | isTrue# (i >=# n) = s
        | otherwise =
            let v = indexDoubleArray# ba_src (off_src +# i)
            in case writeDoubleArray# mba_dst (off_dst +# i) v s of
                 s' -> go (i +# 1#) s'
  in (# go 0# s0, () #)
{-# INLINE copyVectorRaw #-}
