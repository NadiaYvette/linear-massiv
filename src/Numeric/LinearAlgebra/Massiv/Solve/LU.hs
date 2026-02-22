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
import Control.Monad (when, forM)
import GHC.Exts
import GHC.Prim
import GHC.ST (ST(..))
import Data.Primitive.ByteArray (ByteArray(..), MutableByteArray(..), newByteArray, unsafeFreezeByteArray)

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.Solve.Triangular (forwardSubUnit, backSub)
import Numeric.LinearAlgebra.Massiv.Internal.Kernel
  (rawLUEliminateColumn, rawLUEliminateColumnTo, rawSwapRows, rawPivotSearch,
   rawForwardSubUnitPackedSIMD, rawBackSubPackedSIMD,
   rawGemmKernel, rawZeroDoubles)

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
-- For n >= 64, uses panel (blocked) LU factorisation with GEMM trailing update.
luSolveP :: forall n. KnownNat n
         => Matrix n n M.P Double -> Vector n M.P Double -> Vector n M.P Double
luSolveP (MkMatrix a) (MkVector b) =
  let nn = dimVal @n
  in createVector @n @M.P $ \mx -> do
    -- Thaw matrix for in-place LU factorisation
    ma <- M.thawS a
    let mbaA = unwrapMutableByteArray ma
        offA = unwrapMutableByteArrayOffset ma

    -- Phase 1: LU factorisation with partial pivoting
    pivots <- if nn >= 64
      then panelLUFactor mbaA offA nn 32
      else mapM (\k -> do
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

-- | Panel (blocked) LU factorisation with GEMM trailing update.
-- Processes columns in panels of width @nb@.  Within each panel, elimination
-- is restricted to the panel columns; the trailing submatrix is updated via a
-- single GEMM call, converting O(n) column-by-column Level-2 updates into one
-- cache-friendly Level-3 GEMM.
panelLUFactor :: MutableByteArray s -> Int -> Int -> Int -> ST s [(Int, Int)]
panelLUFactor mbaA offA nn nb = go 0 []
  where
    go !k0 !pivAcc
      | k0 >= nn - 1 = pure (reverse pivAcc)
      | otherwise = do
          let !panelEnd = min (k0 + nb) nn
              !actualNb = panelEnd - k0

          -- Factor panel columns k0..panelEnd-1 with restricted trailing update
          panelPivs <- forM [k0..panelEnd-1] $ \k -> do
            pivRow <- rawPivotSearch mbaA offA nn k k
            condM (pivRow /= k) $
              rawSwapRows mbaA offA nn k pivRow 0
            rawLUEliminateColumnTo mbaA offA nn k panelEnd
            pure (k, pivRow)

          -- Apply panel's L to trailing columns: triangular solve for U12
          when (panelEnd < nn) $ do
            rawTriSolvePanelTrail mbaA offA nn k0 panelEnd

            -- GEMM update: A22 -= L21 × U12
            let !mTrail = nn - panelEnd
                !nTrail = nn - panelEnd

            -- Copy L21 to dense buffer (mTrail × actualNb)
            bufL <- newByteArray (mTrail * actualNb * 8)
            rawCopySubmatrixToDense mbaA offA nn panelEnd k0 mTrail actualNb bufL 0

            -- Copy U12 to dense buffer (actualNb × nTrail)
            bufU <- newByteArray (actualNb * nTrail * 8)
            rawCopySubmatrixToDense mbaA offA nn k0 panelEnd actualNb nTrail bufU 0

            -- Freeze for immutable GEMM inputs
            baL <- unsafeFreezeByteArray bufL
            baU <- unsafeFreezeByteArray bufU

            -- GEMM: C = L21 × U12
            bufC <- newByteArray (mTrail * nTrail * 8)
            rawZeroDoubles bufC 0 (mTrail * nTrail)
            rawGemmKernel baL 0 baU 0 bufC 0 mTrail actualNb nTrail

            -- Subtract C from A22
            baC <- unsafeFreezeByteArray bufC
            rawSubtractFromStrided baC 0 nTrail mbaA offA nn panelEnd panelEnd mTrail nTrail

          go panelEnd (reverse panelPivs ++ pivAcc)

-- | Apply unit lower triangular solve from the panel to trailing columns.
-- After factoring panel columns [k0..panelEnd-1] with restricted updates,
-- the trailing columns [panelEnd..n-1] need: for each k in the panel, apply
-- the multipliers to rows k+1..panelEnd-1 of the trailing columns.
-- i.e.  A[i,j] -= A[i,k] * A[k,j]  for  k0 <= k < panelEnd, k < i < panelEnd, panelEnd <= j < n
rawTriSolvePanelTrail :: MutableByteArray s -> Int -> Int -> Int -> Int -> ST s ()
rawTriSolvePanelTrail (MutableByteArray mba) (I# off) (I# n) (I# k0) (I# panelEnd) = ST $ \s0 ->
  let -- For each column k in the panel
      goK k s
        | isTrue# (k >=# panelEnd) = s
        | otherwise =
            let kRowOff = off +# k *# n
            in goI k (k +# 1#) kRowOff s
        where
          -- For each row i in [k+1..panelEnd-1]
          goI k_ i kRowOff s_
            | isTrue# (i >=# panelEnd) = goK (k_ +# 1#) s_
            | otherwise =
                let iRowOff = off +# i *# n
                in case readDoubleArray# mba (iRowOff +# k_) s_ of
                     (# s1, lik #) ->
                       let negLik = negateDouble# lik
                           negLikV = broadcastDoubleX4# negLik
                           jSpan = n -# panelEnd
                           j4End = panelEnd +# (jSpan -# (jSpan `remInt#` 4#))
                           -- SIMD j-loop
                           goJSimd j s__
                             | isTrue# (j >=# j4End) = s__
                             | otherwise =
                                 case readDoubleArrayAsDoubleX4# mba (iRowOff +# j) s__ of
                                   (# s2, aij #) ->
                                     case readDoubleArrayAsDoubleX4# mba (kRowOff +# j) s2 of
                                       (# s3, akj #) ->
                                         let aij' = fmaddDoubleX4# negLikV akj aij
                                         in case writeDoubleArrayAsDoubleX4# mba (iRowOff +# j) aij' s3 of
                                              s4 -> goJSimd (j +# 4#) s4
                           -- Scalar cleanup
                           goJScalar j s__
                             | isTrue# (j >=# n) = s__
                             | otherwise =
                                 case readDoubleArray# mba (iRowOff +# j) s__ of
                                   (# s2, aij #) ->
                                     case readDoubleArray# mba (kRowOff +# j) s2 of
                                       (# s3, akj #) ->
                                         case writeDoubleArray# mba (iRowOff +# j) (aij +## negLik *## akj) s3 of
                                           s4 -> goJScalar (j +# 1#) s4
                       in goI k_ (i +# 1#) kRowOff (goJScalar j4End (goJSimd panelEnd s1))
  in (# goK k0 s0, () #)
{-# INLINE rawTriSolvePanelTrail #-}

-- | Copy a submatrix A[rowStart..rowStart+m-1, colStart..colStart+k-1] (stride n)
-- into a dense buffer (stride k).
rawCopySubmatrixToDense :: MutableByteArray s -> Int -> Int  -- src, offset, n
                        -> Int -> Int -> Int -> Int          -- rowStart, colStart, m, k
                        -> MutableByteArray s -> Int          -- dst, dstOffset
                        -> ST s ()
rawCopySubmatrixToDense (MutableByteArray mba_src) (I# off_src) (I# n)
                        (I# rowStart) (I# colStart) (I# m) (I# k)
                        (MutableByteArray mba_dst) (I# off_dst) = ST $ \s0 ->
  let goI i s
        | isTrue# (i >=# m) = s
        | otherwise =
            let srcRow = off_src +# (rowStart +# i) *# n +# colStart
                dstRow = off_dst +# i *# k
                span_ = k -# (k `remInt#` 4#)
                goSimd j s_
                  | isTrue# (j >=# span_) = s_
                  | otherwise =
                      case readDoubleArrayAsDoubleX4# mba_src (srcRow +# j) s_ of
                        (# s1, v #) ->
                          case writeDoubleArrayAsDoubleX4# mba_dst (dstRow +# j) v s1 of
                            s2 -> goSimd (j +# 4#) s2
                goScalar j s_
                  | isTrue# (j >=# k) = s_
                  | otherwise =
                      case readDoubleArray# mba_src (srcRow +# j) s_ of
                        (# s1, v #) ->
                          case writeDoubleArray# mba_dst (dstRow +# j) v s1 of
                            s2 -> goScalar (j +# 1#) s2
            in goI (i +# 1#) (goScalar span_ (goSimd 0# s))
  in (# goI 0# s0, () #)
{-# INLINE rawCopySubmatrixToDense #-}

-- | Subtract a dense matrix C (m × k, stride srcStride) from a strided submatrix
-- A[rowStart..rowStart+m-1, colStart..colStart+k-1] (stride n).
rawSubtractFromStrided :: ByteArray -> Int -> Int              -- src, srcOffset, srcStride
                       -> MutableByteArray s -> Int -> Int     -- dst, dstOffset, n
                       -> Int -> Int -> Int -> Int             -- rowStart, colStart, m, k
                       -> ST s ()
rawSubtractFromStrided (ByteArray ba_src) (I# off_src) (I# srcStride)
                       (MutableByteArray mba_dst) (I# off_dst) (I# n)
                       (I# rowStart) (I# colStart) (I# m) (I# k) = ST $ \s0 ->
  let goI i s
        | isTrue# (i >=# m) = s
        | otherwise =
            let srcRow = off_src +# i *# srcStride
                dstRow = off_dst +# (rowStart +# i) *# n +# colStart
                span_ = k -# (k `remInt#` 4#)
                goSimd j s_
                  | isTrue# (j >=# span_) = s_
                  | otherwise =
                      case readDoubleArrayAsDoubleX4# mba_dst (dstRow +# j) s_ of
                        (# s1, aij #) ->
                          let cij = indexDoubleArrayAsDoubleX4# ba_src (srcRow +# j)
                              aij' = plusDoubleX4# aij (negateDoubleX4# cij)
                          in case writeDoubleArrayAsDoubleX4# mba_dst (dstRow +# j) aij' s1 of
                               s2 -> goSimd (j +# 4#) s2
                goScalar j s_
                  | isTrue# (j >=# k) = s_
                  | otherwise =
                      case readDoubleArray# mba_dst (dstRow +# j) s_ of
                        (# s1, aij #) ->
                          let cij = indexDoubleArray# ba_src (srcRow +# j)
                          in case writeDoubleArray# mba_dst (dstRow +# j) (aij -## cij) s1 of
                               s2 -> goScalar (j +# 1#) s2
            in goI (i +# 1#) (goScalar span_ (goSimd 0# s))
  in (# goI 0# s0, () #)
{-# INLINE rawSubtractFromStrided #-}

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
