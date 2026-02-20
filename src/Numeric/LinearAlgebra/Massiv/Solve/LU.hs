{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE MagicHash #-}
{-# LANGUAGE UnboxedTuples #-}
{-# LANGUAGE BangPatterns #-}

-- |
-- Module      : Numeric.LinearAlgebra.Massiv.Solve.LU
-- Copyright   : (c) Nadia Chambers 2026
-- License     : BSD-3-Clause
-- Maintainer  : nadia.chambers@iohk.io
-- Stability   : experimental
--
-- = LU Factorization
--
-- LU decomposition with and without partial pivoting, plus derived
-- operations (linear solve and determinant), following Golub & Van Loan,
-- /Matrix Computations/, 4th edition (GVL4), Sections 3.2 and 3.4,
-- pp. 114--131.
--
-- Given an \(n \times n\) matrix \(A\), the factorization produces
--
-- \[
--   PA = LU
-- \]
--
-- where \(P\) is a permutation matrix, \(L\) is unit lower triangular, and
-- \(U\) is upper triangular (GVL4 Theorem 3.4.1, p. 125).
--
-- Without pivoting (\(P = I\)) the factorization exists if and only if all
-- leading principal submatrices of \(A\) are nonsingular
-- (GVL4 Theorem 3.2.1, p. 116). Partial pivoting guarantees existence for
-- any nonsingular \(A\) and improves numerical stability by bounding the
-- growth factor (GVL4 Section 3.4.6).
--
-- +-------------------+----------------------------------+-------------------------------+
-- | Function          | Algorithm                        | Reference                     |
-- +===================+==================================+===============================+
-- | 'lu'              | LU with partial pivoting         | GVL4 Algorithm 3.4.1, p. 126  |
-- +-------------------+----------------------------------+-------------------------------+
-- | 'luNoPivot'       | Outer-product LU (no pivoting)   | GVL4 Algorithm 3.2.1, p. 115  |
-- +-------------------+----------------------------------+-------------------------------+
-- | 'luSolve'         | Solve via \(PA = LU\)            | GVL4 Section 3.2, p. 118      |
-- +-------------------+----------------------------------+-------------------------------+
-- | 'det'             | Determinant via \(PA = LU\)      | GVL4 Section 3.2, p. 120      |
-- +-------------------+----------------------------------+-------------------------------+
--
-- == Complexity
--
-- The factorization requires \(O(2n^3/3)\) flops (GVL4 p. 118). Each
-- subsequent triangular solve adds \(O(n^2)\) flops.
--
-- == Type Safety
--
-- Matrix dimensions are tracked at the type level via 'KnownNat', so the
-- compiler statically enforces that the coefficient matrix is square and
-- that right-hand side vectors have conforming length.
module Numeric.LinearAlgebra.Massiv.Solve.LU
  ( -- * LU factorization
    lu
  , luNoPivot
    -- * Solving with LU (\(Ax = b\))
  , luSolve
  , luSolveP
    -- * Determinant
  , det
  ) where

import qualified Data.Massiv.Array as M
import Data.Massiv.Array (Ix1, Ix2(..), Sz(..), unwrapByteArray, unwrapByteArrayOffset,
                          unwrapMutableByteArray, unwrapMutableByteArrayOffset)
import GHC.TypeNats (KnownNat)
import Data.Ord (comparing)
import Data.List (maximumBy)
import GHC.Exts
import GHC.ST (ST(..))
import Data.Primitive.ByteArray (ByteArray(..), MutableByteArray(..))

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.Solve.Triangular (forwardSubUnit, backSub)
import Numeric.LinearAlgebra.Massiv.Internal.Kernel
  (rawLUEliminateColumn, rawSwapRows, rawPivotSearch,
   rawForwardSubUnitPackedSIMD, rawBackSubPackedSIMD)

-- | LU factorization with partial pivoting (GVL4 Algorithm 3.4.1, p. 126).
--
-- Given an \(n \times n\) matrix \(A\), computes the factorization
--
-- \[
--   PA = LU
-- \]
--
-- where
--
-- * \(P\) is a permutation matrix (returned as a pivot-index vector of type
--   @Array P Ix1 Int@),
-- * \(L\) is unit lower triangular (stored /below/ the diagonal of the
--   returned packed matrix), and
-- * \(U\) is upper triangular (stored /on and above/ the diagonal).
--
-- Partial pivoting selects the entry of largest absolute value in the
-- current column as the pivot, guaranteeing existence for any nonsingular
-- \(A\) (GVL4 Theorem 3.4.1, p. 125).
--
-- ==== Type-safety guarantees
--
-- The 'KnownNat' constraint on \(n\) statically ensures the matrix is
-- square. The 'Ord' constraint is required for pivot selection.
--
-- ==== Complexity
--
-- \(O(2n^3/3)\) flops (GVL4 p. 118).
--
-- ==== Reference
--
-- Golub & Van Loan, /Matrix Computations/, 4th ed., Algorithm 3.4.1
-- (Outer Product LU with Partial Pivoting), p. 126.
lu :: forall n r e. (KnownNat n, M.Manifest r e, Fractional e, Ord e)
   => Matrix n n r e -> (Matrix n n r e, M.Array M.P Ix1 Int)
lu (MkMatrix a) =
  let nn = dimVal @n
      (pivArr, luArr) = M.withMArrayST a $ \ma -> do
        piv <- M.newMArray @M.P (M.Sz1 nn) 0
        mapM_ (\i -> M.write_ piv i i) [0..nn-1]

        mapM_ (\k -> do
          -- Find pivot row
          vals <- mapM (\i -> do
            v <- M.readM ma (i :. k)
            pure (i, abs v)
            ) [k..nn-1]
          let (pivRow, _) = maximumBy (comparing snd) vals

          -- Swap rows k and pivRow
          condM (pivRow /= k) $ do
            pk <- M.readM piv k
            pp <- M.readM piv pivRow
            M.write_ piv k pp
            M.write_ piv pivRow pk
            mapM_ (\j -> do
              akj <- M.readM ma (k :. j)
              apj <- M.readM ma (pivRow :. j)
              M.write_ ma (k :. j) apj
              M.write_ ma (pivRow :. j) akj
              ) [0..nn-1]

          -- Compute multipliers and update submatrix
          akk <- M.readM ma (k :. k)
          condM (akk /= 0) $
            mapM_ (\i -> do
              aik <- M.readM ma (i :. k)
              let mult = aik / akk
              M.write_ ma (i :. k) mult
              mapM_ (\j -> do
                aij <- M.readM ma (i :. j)
                akj <- M.readM ma (k :. j)
                M.write_ ma (i :. j) (aij - mult * akj)
                ) [k+1..nn-1]
              ) [k+1..nn-1]
          ) [0..nn-2]

        M.freezeS piv

  in (MkMatrix luArr, pivArr)

-- | LU factorization without pivoting (GVL4 Algorithm 3.2.1, p. 115).
--
-- Given an \(n \times n\) matrix \(A\) whose leading principal submatrices
-- are all nonsingular, computes the factorization \(A = LU\) where \(L\) is
-- unit lower triangular and \(U\) is upper triangular. Both factors are
-- packed into a single returned matrix: \(L\) occupies the strictly lower
-- triangle (the unit diagonal is implicit) and \(U\) occupies the upper
-- triangle including the diagonal.
--
-- __Precondition.__ All leading principal submatrices
-- \(A(1{:}k, 1{:}k)\), \(k = 1, \ldots, n\), must be nonsingular
-- (GVL4 Theorem 3.2.1, p. 116). Violating this precondition results in
-- division by zero.
--
-- ==== Type-safety guarantees
--
-- The 'KnownNat' constraint on \(n\) statically ensures the input is a
-- square matrix.
--
-- ==== Complexity
--
-- \(O(2n^3/3)\) flops (GVL4 p. 118).
--
-- ==== Reference
--
-- Golub & Van Loan, /Matrix Computations/, 4th ed., Algorithm 3.2.1
-- (Outer Product LU Factorization), p. 115.
luNoPivot :: forall n r e. (KnownNat n, M.Manifest r e, Fractional e)
          => Matrix n n r e -> Matrix n n r e
luNoPivot (MkMatrix a) =
  let nn = dimVal @n
  in MkMatrix $ snd $ M.withMArrayST a $ \ma ->
    mapM_ (\k -> do
      akk <- M.readM ma (k :. k)
      mapM_ (\i -> do
        aik <- M.readM ma (i :. k)
        M.write_ ma (i :. k) (aik / akk)
        ) [k+1..nn-1]
      mapM_ (\j ->
        mapM_ (\i -> do
          aij <- M.readM ma (i :. j)
          aik <- M.readM ma (i :. k)
          akj <- M.readM ma (k :. j)
          M.write_ ma (i :. j) (aij - aik * akj)
          ) [k+1..nn-1]
        ) [k+1..nn-1]
      ) [0..nn-2]

-- | Solve \(Ax = b\) using LU factorization with partial pivoting
-- (GVL4 Section 3.2, p. 118).
--
-- The algorithm proceeds in three stages:
--
-- 1. Factor \(PA = LU\) via 'lu' (Algorithm 3.4.1).
-- 2. Solve \(Ly = Pb\) by forward substitution ('forwardSubUnit').
-- 3. Solve \(Ux = y\) by back substitution ('backSub').
--
-- ==== Type-safety guarantees
--
-- 'KnownNat' \(n\) enforces that \(A\) is \(n \times n\) and \(b\) has
-- length \(n\) at compile time.
--
-- ==== Complexity
--
-- \(O(2n^3/3)\) flops for the factorization plus \(O(n^2)\) flops for each
-- of the two triangular solves (GVL4 p. 118).
--
-- ==== Reference
--
-- Golub & Van Loan, /Matrix Computations/, 4th ed., Section 3.2,
-- pp. 114--120.
luSolve :: forall n r e. (KnownNat n, M.Manifest r e, Fractional e, Ord e)
        => Matrix n n r e -> Vector n r e -> Vector n r e
luSolve a b =
  let (luMat, pivArr) = lu a
      -- Extract L (unit lower triangular)
      l = makeMatrix @n @n @r $ \i j ->
        if i == j then 1
        else if i > j then luMat ! (i, j)
        else 0
      -- Extract U (upper triangular)
      u = makeMatrix @n @n @r $ \i j ->
        if i <= j then luMat ! (i, j)
        else 0
      -- Apply permutation to b: pb = P·b
      pb = makeVector @n @r $ \i ->
        b !. M.index' pivArr i
      -- Solve Ly = Pb
      y = forwardSubUnit l pb
      -- Solve Ux = y
  in backSub u y

-- | Specialised LU solve for @P Double@.
-- Does LU factorisation + solve entirely using raw ByteArray# primops,
-- avoiding L/U matrix reconstruction and massiv's per-element overhead.
luSolveP :: forall n. KnownNat n
         => Matrix n n M.P Double -> Vector n M.P Double -> Vector n M.P Double
luSolveP (MkMatrix a) (MkVector b) =
  let nn = dimVal @n
  in createVector @n @M.P $ \mx -> do
    -- Thaw matrix for in-place LU factorisation
    ma <- M.thawS a
    let mbaA = unwrapMutableByteArray ma
        offA = unwrapMutableByteArrayOffset ma

    -- Phase 1: LU factorisation with partial pivoting using raw kernels
    pivots <- mapM (\k -> do
      pivRow <- rawPivotSearch mbaA offA nn k k
      condM (pivRow /= k) $
        rawSwapRows mbaA offA nn k pivRow 0
      rawLUEliminateColumn mbaA offA nn k
      pure (k, pivRow)
      ) [0..nn-2]

    -- Phase 2: Freeze LU and prepare RHS
    frozenLU <- M.freezeS ma
    let baLU = unwrapByteArray frozenLU
        offLU = unwrapByteArrayOffset frozenLU

    -- Copy b into the output vector mx
    let mbaX = unwrapMutableByteArray mx
        offX = unwrapMutableByteArrayOffset mx
    copyVector b mbaX offX nn

    -- Apply pivot permutation to x
    applyPivotsForward mbaX offX pivots

    -- Phase 3: Forward substitution (unit lower triangular, SIMD dot-product)
    rawForwardSubUnitPackedSIMD baLU offLU nn mbaX offX

    -- Phase 4: Back substitution (upper triangular, SIMD dot-product)
    rawBackSubPackedSIMD baLU offLU nn mbaX offX
{-# NOINLINE luSolveP #-}

-- | Copy an immutable P vector into a mutable byte array.
copyVector :: M.Array M.P Ix1 Double -> MutableByteArray s -> Int -> Int -> ST s ()
copyVector src (MutableByteArray mba_dst) (I# off_dst) (I# n) = ST $ \s0 ->
  let ba_src = case unwrapByteArray src of ByteArray ba -> ba
      off_src = case unwrapByteArrayOffset src of I# o -> o
      go i s
        | isTrue# (i >=# n) = s
        | otherwise =
            let v = indexDoubleArray# ba_src (off_src +# i)
            in case writeDoubleArray# mba_dst (off_dst +# i) v s of
                 s' -> go (i +# 1#) s'
  in (# go 0# s0, () #)
{-# INLINE copyVector #-}

-- | Apply pivot permutation to a vector (forward direction).
applyPivotsForward :: MutableByteArray s -> Int -> [(Int, Int)] -> ST s ()
applyPivotsForward (MutableByteArray mba) (I# off) pivots = ST $ \s0 ->
  let go [] s = s
      go ((I# k, I# pivRow) : rest) s
        | isTrue# (k ==# pivRow) = go rest s
        | otherwise =
            case readDoubleArray# mba (off +# k) s of
              (# s1, vk #) ->
                case readDoubleArray# mba (off +# pivRow) s1 of
                  (# s2, vp #) ->
                    case writeDoubleArray# mba (off +# k) vp s2 of
                      s3 -> case writeDoubleArray# mba (off +# pivRow) vk s3 of
                              s4 -> go rest s4
  in (# go pivots s0, () #)
{-# INLINE applyPivotsForward #-}

-- | Compute the determinant of an \(n \times n\) matrix via LU factorization
-- (GVL4 Section 3.2, p. 120).
--
-- ==== Mathematical definition
--
-- From \(PA = LU\) it follows that
--
-- \[
--   \det(A) = (-1)^s \prod_{i=1}^{n} u_{ii}
-- \]
--
-- where \(s\) is the number of row transpositions performed during partial
-- pivoting.
--
-- ==== Type-safety guarantees
--
-- 'KnownNat' \(n\) ensures \(A\) is square at compile time.
--
-- ==== Complexity
--
-- \(O(2n^3/3)\) flops, dominated by the LU factorization.
--
-- ==== Reference
--
-- Golub & Van Loan, /Matrix Computations/, 4th ed., Section 3.2, p. 120.
det :: forall n r e. (KnownNat n, M.Manifest r e, Fractional e, Ord e)
    => Matrix n n r e -> e
det a =
  let nn = dimVal @n
      (luMat, pivArr) = lu a
      -- Product of U diagonal
      diagProd = foldl' (\acc i -> acc * (luMat ! (i, i))) 1 [0..nn-1]
      -- Count transpositions: number of i where piv[i] /= i
      pivList = map (M.index' pivArr) [0..nn-1]
      nswaps = countSwaps pivList
      sign = if even nswaps then 1 else -1
  in sign * diagProd

-- | Count the number of swaps in a permutation.
countSwaps :: [Int] -> Int
countSwaps perm = go (zip [0..] perm) 0
  where
    go [] n = n
    go ((i,p):rest) n
      | i == p = go rest n
      | otherwise =
          -- Swap p into position i by finding where i is
          let rest' = map (\(idx, v) -> if v == i then (idx, p) else (idx, v)) rest
          in go rest' (n + 1)

-- | Conditional monadic action.
condM :: Applicative m => Bool -> m () -> m ()
condM True act = act
condM False _ = pure ()
