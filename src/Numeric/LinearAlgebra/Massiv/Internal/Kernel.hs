{-# LANGUAGE MagicHash #-}
{-# LANGUAGE UnboxedTuples #-}
{-# LANGUAGE BangPatterns #-}

-- | Raw ByteArray# + AVX2 SIMD kernels for hot inner loops.
--
-- These functions bypass massiv's per-element abstraction layer
-- (M.readM / M.write_ / mapM_) and operate directly on the underlying
-- ByteArray# / MutableByteArray# storage, using GHC's DoubleX4# primops
-- for 256-bit AVX2 SIMD where possible.
module Numeric.LinearAlgebra.Massiv.Internal.Kernel
  ( -- * BLAS-1: dot product
    rawDot
    -- * BLAS-2: matrix-vector multiply
  , rawGemv
    -- * BLAS-3: matrix multiply (GEMM)
  , rawGemmKernel
  , rawGemmBISlice
    -- * QR helpers (immutable)
  , rawSumSqRange
  , rawSumProdRange
  , rawHouseholderApplyCol
  , rawQAccumCol
    -- * Eigen / tridiag helpers
  , rawApplyGivensRows
  , rawSymRank2Update
    -- * LU kernels
  , rawLUEliminateColumn
  , rawSwapRows
  , rawPivotSearch
  , rawForwardSubUnitPacked
  , rawBackSubPacked
    -- * Cholesky kernels
  , rawCholColumn
  , rawCholColumnSIMD
  , rawForwardSubCholPacked
  , rawBackSubCholTPacked
    -- * QR mutable kernels
  , rawMutSumSqColumn
  , rawMutSumProdColumns
  , rawMutHouseholderApply
  , rawMutQAccum
    -- * Tridiagonalisation mutable kernels
  , rawMutSymMatvecSub
  , rawMutSymRank2Update
  , rawMutTridiagQAccum
    -- * Eigen mutable kernels
  , rawMutApplyGivensColumns
    -- * SVD / bidiagonalisation kernels
  , rawMutHouseholderApplyRow
  , rawMutSumSqRow
  ) where

import GHC.Exts
import GHC.Prim
import GHC.ST (ST(..))
import GHC.Types (Double(..), Int(..))
import Data.Primitive.ByteArray (ByteArray(..), MutableByteArray(..))

-- --------------------------------------------------------------------------
-- BLAS-1: dot product  (DoubleX4# FMA, scalar cleanup)
-- --------------------------------------------------------------------------

-- | Dot product of two Double vectors stored in ByteArrays.
-- @rawDot ba1 off1 ba2 off2 n@ computes Σ ba1[off1+i] * ba2[off2+i] for i in [0..n-1].
rawDot :: ByteArray -> Int -> ByteArray -> Int -> Int -> Double
rawDot (ByteArray ba1) (I# off1) (ByteArray ba2) (I# off2) (I# n) =
  D# (rawDot# ba1 off1 ba2 off2 n)
{-# INLINE rawDot #-}

rawDot# :: ByteArray# -> Int# -> ByteArray# -> Int# -> Int# -> Double#
rawDot# ba1 off1 ba2 off2 n =
  let n4 = n -# (n `remInt#` 4#)
      -- SIMD phase: accumulate 4 doubles at a time
      goSimd i acc
        | isTrue# (i >=# n4) = acc
        | otherwise =
            let va = indexDoubleArrayAsDoubleX4# ba1 (off1 +# i)
                vb = indexDoubleArrayAsDoubleX4# ba2 (off2 +# i)
            in goSimd (i +# 4#) (fmaddDoubleX4# va vb acc)
      acc4 = goSimd 0# (broadcastDoubleX4# 0.0##)
      !(# a, b, c, d #) = unpackDoubleX4# acc4
      simdSum = a +## b +## c +## d
      -- Scalar cleanup for remainder
      goScalar i acc
        | isTrue# (i >=# n) = acc
        | otherwise =
            let x = indexDoubleArray# ba1 (off1 +# i)
                y = indexDoubleArray# ba2 (off2 +# i)
            in goScalar (i +# 1#) (acc +## x *## y)
  in goScalar n4 simdSum
{-# NOINLINE rawDot# #-}

-- --------------------------------------------------------------------------
-- BLAS-2: matrix-vector multiply
-- --------------------------------------------------------------------------

-- | @rawGemv ba_a off_a n_cols ba_x off_x mba_y off_y n_rows@ computes
-- y[i] = Σ_j A[i,j] * x[j] for i in [0..n_rows-1], j in [0..n_cols-1].
-- A is row-major with stride = n_cols.
rawGemv :: ByteArray -> Int -> Int
        -> ByteArray -> Int
        -> MutableByteArray s -> Int -> Int
        -> ST s ()
rawGemv (ByteArray ba_a) (I# off_a) (I# ncols)
        (ByteArray ba_x) (I# off_x)
        (MutableByteArray mba_y) (I# off_y) (I# nrows) = ST $ \s0 ->
  let go i s
        | isTrue# (i >=# nrows) = s
        | otherwise =
            let rowOff = off_a +# i *# ncols
                dot = rawDot# ba_a rowOff ba_x off_x ncols
            in case writeDoubleArray# mba_y (off_y +# i) dot s of
                 s' -> go (i +# 1#) s'
  in (# go 0# s0, () #)
{-# INLINE rawGemv #-}

-- --------------------------------------------------------------------------
-- BLAS-3: tiled ikj GEMM kernel
-- --------------------------------------------------------------------------

-- | @rawGemmKernel ba_a off_a ba_b off_b mba_c off_c m k n@ computes
-- C += A * B where A is m×k, B is k×n, C is m×n (all row-major).
-- C must be pre-initialised (e.g. to zero, or to β*C for gemm).
rawGemmKernel :: ByteArray -> Int -> ByteArray -> Int
              -> MutableByteArray s -> Int
              -> Int -> Int -> Int -> ST s ()
rawGemmKernel ba offA bb offB mc offC m k n =
  rawGemmBISlice ba offA bb offB mc offC 0 m m k n
{-# INLINE rawGemmKernel #-}

-- | @rawGemmBISlice ba_a off_a ba_b off_b mba_c off_c biStart biEnd m k n@
-- computes C[biStart..biEnd-1, :] += A[biStart..biEnd-1, :] * B
-- where A is m×k, B is k×n, C is m×n (all row-major).
-- Only the rows [biStart, biEnd) of C are written; all of B is read.
-- This enables parallel GEMM by partitioning the row range across threads.
rawGemmBISlice :: ByteArray -> Int -> ByteArray -> Int
               -> MutableByteArray s -> Int
               -> Int -> Int -> Int -> Int -> Int -> ST s ()
rawGemmBISlice (ByteArray ba_a) (I# off_a) (ByteArray ba_b) (I# off_b)
               (MutableByteArray mba_c) (I# off_c)
               (I# biStart) (I# biEnd) (I# _m) (I# k) (I# n) = ST $ \s0 ->
  let bs = 64#

      goBI bi s
        | isTrue# (bi >=# biEnd) = s
        | otherwise =
            let iEnd = minI (bi +# bs) biEnd
            in goBI (bi +# bs) (goBK bi iEnd 0# s)

      goBK bi iEnd bk s
        | isTrue# (bk >=# k) = s
        | otherwise =
            let kEnd = minI (bk +# bs) k
            in goBK bi iEnd (bk +# bs) (goBJ bi iEnd bk kEnd 0# s)

      goBJ bi iEnd bk kEnd bj s
        | isTrue# (bj >=# n) = s
        | otherwise =
            let jEnd = minI (bj +# bs) n
            in goBJ bi iEnd bk kEnd (bj +# bs) (innerIKJ bi iEnd bk kEnd bj jEnd s)

      innerIKJ bi iEnd bk kEnd bj jEnd s0_ = goI bi s0_
        where
          goI i s
            | isTrue# (i >=# iEnd) = s
            | otherwise = goI (i +# 1#) (goK i bk s)

          goK i kk s
            | isTrue# (kk >=# kEnd) = s
            | otherwise =
                let aik = indexDoubleArray# ba_a (off_a +# i *# k +# kk)
                    aikV = broadcastDoubleX4# aik
                    bRowOff = off_b +# kk *# n
                    cRowOff = off_c +# i *# n
                    jSpan = jEnd -# bj
                    j4End = bj +# (jSpan -# (jSpan `remInt#` 4#))
                in goK i (kk +# 1#) (goJSimd bj j4End aikV bRowOff cRowOff
                                      (goJScalar j4End jEnd aik bRowOff cRowOff s))

          goJSimd j j4End aikV bRowOff cRowOff s
            | isTrue# (j >=# j4End) = s
            | otherwise =
                let bv = indexDoubleArrayAsDoubleX4# ba_b (bRowOff +# j)
                in case readDoubleArrayAsDoubleX4# mba_c (cRowOff +# j) s of
                     (# s', cv #) ->
                       let cv' = fmaddDoubleX4# aikV bv cv
                       in case writeDoubleArrayAsDoubleX4# mba_c (cRowOff +# j) cv' s' of
                            s'' -> goJSimd (j +# 4#) j4End aikV bRowOff cRowOff s''

          goJScalar j jEnd_ aik bRowOff cRowOff s
            | isTrue# (j >=# jEnd_) = s
            | otherwise =
                let bkj = indexDoubleArray# ba_b (bRowOff +# j)
                in case readDoubleArray# mba_c (cRowOff +# j) s of
                     (# s', cij #) ->
                       case writeDoubleArray# mba_c (cRowOff +# j) (cij +## aik *## bkj) s' of
                         s'' -> goJScalar (j +# 1#) jEnd_ aik bRowOff cRowOff s''

  in (# goBI biStart s0, () #)
{-# INLINE rawGemmBISlice #-}

-- --------------------------------------------------------------------------
-- QR helpers
-- --------------------------------------------------------------------------

-- | Sum of squares: Σ arr[off+i]^2 for i in [from..to-1]
rawSumSqRange :: ByteArray -> Int -> Int -> Int -> Double
rawSumSqRange (ByteArray ba) (I# off) (I# from_) (I# to) =
  D# (goSumSq from_ 0.0##)
  where
    goSumSq i acc
      | isTrue# (i >=# to) = acc
      | otherwise =
          let x = indexDoubleArray# ba (off +# i)
          in goSumSq (i +# 1#) (acc +## x *## x)
{-# INLINE rawSumSqRange #-}

-- | Dot product of a column slice: Σ arr1[off1+i*stride1] * arr2[off2+i*stride2]
-- for i in [from..to-1]. Used for column-wise access patterns in QR.
rawSumProdRange :: ByteArray -> Int -> Int
                -> ByteArray -> Int -> Int
                -> Int -> Int -> Double
rawSumProdRange (ByteArray ba1) (I# off1) (I# stride1)
                (ByteArray ba2) (I# off2) (I# stride2)
                (I# from_) (I# to) =
  D# (goSumProd from_ 0.0##)
  where
    goSumProd i acc
      | isTrue# (i >=# to) = acc
      | otherwise =
          let x = indexDoubleArray# ba1 (off1 +# i *# stride1)
              y = indexDoubleArray# ba2 (off2 +# i *# stride2)
          in goSumProd (i +# 1#) (acc +## x *## y)
{-# INLINE rawSumProdRange #-}

-- | Apply Householder reflector to one column of a mutable matrix.
-- @rawHouseholderApplyCol mba_r off_r ncols ba_v off_v beta from to col@
-- computes w = β * Σ_{i=from}^{to-1} v[i] * R[i, col], then
-- R[i, col] -= v[i] * w for i in [from..to-1].
rawHouseholderApplyCol :: MutableByteArray s -> Int -> Int
                       -> ByteArray -> Int -> Double
                       -> Int -> Int -> Int -> ST s ()
rawHouseholderApplyCol (MutableByteArray mba_r) (I# off_r) (I# ncols)
                       (ByteArray ba_v) (I# off_v) (D# beta)
                       (I# from_) (I# to) (I# col) = ST $ \s0 ->
  -- Phase 1: compute w = β * Σ v[i] * R[i, col]
  let goSum i acc s
        | isTrue# (i >=# to) = (# s, beta *## acc #)
        | otherwise =
            let vi = indexDoubleArray# ba_v (off_v +# i)
            in case readDoubleArray# mba_r (off_r +# i *# ncols +# col) s of
                 (# s', rij #) -> goSum (i +# 1#) (acc +## vi *## rij) s'
  in case goSum from_ 0.0## s0 of
       (# s1, w #) ->
         -- Phase 2: R[i, col] -= v[i] * w
         let goUpdate i s
               | isTrue# (i >=# to) = s
               | otherwise =
                   let vi = indexDoubleArray# ba_v (off_v +# i)
                   in case readDoubleArray# mba_r (off_r +# i *# ncols +# col) s of
                        (# s', rij #) ->
                          case writeDoubleArray# mba_r (off_r +# i *# ncols +# col) (rij -## vi *## w) s' of
                            s'' -> goUpdate (i +# 1#) s''
         in (# goUpdate from_ s1, () #)
{-# INLINE rawHouseholderApplyCol #-}

-- | Update one row of Q during backward accumulation.
-- @rawQAccumCol mba_q off_q ncols ba_v off_v beta from to row@
-- computes qi = β * Σ_{k=from}^{to-1} Q[row, k] * v[k], then
-- Q[row, k] -= qi * v[k] for all k in [from..to-1].
rawQAccumCol :: MutableByteArray s -> Int -> Int
             -> ByteArray -> Int -> Double
             -> Int -> Int -> Int -> ST s ()
rawQAccumCol (MutableByteArray mba_q) (I# off_q) (I# ncols)
             (ByteArray ba_v) (I# off_v) (D# beta)
             (I# from_) (I# to) (I# row) = ST $ \s0 ->
  -- Phase 1: compute qi = β * Σ Q[row, k] * v[k]
  let goSum k acc s
        | isTrue# (k >=# to) = (# s, beta *## acc #)
        | otherwise =
            let vk = indexDoubleArray# ba_v (off_v +# k)
            in case readDoubleArray# mba_q (off_q +# row *# ncols +# k) s of
                 (# s', qrk #) -> goSum (k +# 1#) (acc +## vk *## qrk) s'
  in case goSum from_ 0.0## s0 of
       (# s1, qi #) ->
         -- Phase 2: Q[row, k] -= qi * v[k]
         let goUpdate k s
               | isTrue# (k >=# to) = s
               | otherwise =
                   let vk = indexDoubleArray# ba_v (off_v +# k)
                   in case readDoubleArray# mba_q (off_q +# row *# ncols +# k) s of
                        (# s', qrk #) ->
                          case writeDoubleArray# mba_q (off_q +# row *# ncols +# k) (qrk -## qi *## vk) s' of
                            s'' -> goUpdate (k +# 1#) s''
         in (# goUpdate from_ s1, () #)
{-# INLINE rawQAccumCol #-}

-- --------------------------------------------------------------------------
-- Eigen / tridiag helpers
-- --------------------------------------------------------------------------

-- | Apply Givens rotation to two rows of a mutable matrix.
-- @rawApplyGivensRows mba off ncols cosθ sinθ row_p row_q from to@
-- For each column j in [from..to-1]:
--   tmp        =  c * M[row_p, j] + s * M[row_q, j]
--   M[row_q, j] = -s * M[row_p, j] + c * M[row_q, j]
--   M[row_p, j] = tmp
rawApplyGivensRows :: MutableByteArray s -> Int -> Int
                   -> Double -> Double -> Int -> Int
                   -> Int -> Int -> ST s ()
rawApplyGivensRows (MutableByteArray mba) (I# off) (I# ncols)
                   (D# c_) (D# s_) (I# row_p) (I# row_q)
                   (I# from_) (I# to) = ST $ \s0 ->
  let pOff = off +# row_p *# ncols
      qOff = off +# row_q *# ncols
      jSpan = to -# from_
      j4End = from_ +# (jSpan -# (jSpan `remInt#` 4#))
      cV = broadcastDoubleX4# c_
      sV = broadcastDoubleX4# s_
      nsV = negateDoubleX4# sV

      goSimd j s
        | isTrue# (j >=# j4End) = s
        | otherwise =
            case readDoubleArrayAsDoubleX4# mba (pOff +# j) s of
              (# s1, pv #) ->
                case readDoubleArrayAsDoubleX4# mba (qOff +# j) s1 of
                  (# s2, qv #) ->
                    let tmp = fmaddDoubleX4# cV pv (timesDoubleX4# sV qv)
                        q'  = fmaddDoubleX4# nsV pv (timesDoubleX4# cV qv)
                    in case writeDoubleArrayAsDoubleX4# mba (pOff +# j) tmp s2 of
                         s3 -> case writeDoubleArrayAsDoubleX4# mba (qOff +# j) q' s3 of
                                 s4 -> goSimd (j +# 4#) s4

      goScalar j s
        | isTrue# (j >=# to) = s
        | otherwise =
            case readDoubleArray# mba (pOff +# j) s of
              (# s1, pj #) ->
                case readDoubleArray# mba (qOff +# j) s1 of
                  (# s2, qj #) ->
                    let tmp = c_ *## pj +## s_ *## qj
                        qj' = negateDouble# s_ *## pj +## c_ *## qj
                    in case writeDoubleArray# mba (pOff +# j) tmp s2 of
                         s3 -> case writeDoubleArray# mba (qOff +# j) qj' s3 of
                                 s4 -> goScalar (j +# 1#) s4

  in (# goScalar j4End (goSimd from_ s0), () #)
{-# INLINE rawApplyGivensRows #-}

-- | Symmetric rank-2 update on a mutable matrix.
-- @rawSymRank2Update mba off n ba_v off_v ba_w off_w from to@
-- For i in [from..to-1], j in [from..i]:
--   T[i,j] -= v[i]*w[j] + w[i]*v[j]
--   T[j,i] = T[i,j]   (maintain symmetry)
rawSymRank2Update :: MutableByteArray s -> Int -> Int
                  -> ByteArray -> Int -> ByteArray -> Int
                  -> Int -> Int -> ST s ()
rawSymRank2Update (MutableByteArray mba) (I# off) (I# n)
                  (ByteArray ba_v) (I# off_v) (ByteArray ba_w) (I# off_w)
                  (I# from_) (I# to) = ST $ \s0 ->
  let goI i s
        | isTrue# (i >=# to) = s
        | otherwise =
            let vi = indexDoubleArray# ba_v (off_v +# i)
                wi = indexDoubleArray# ba_w (off_w +# i)
            in goI (i +# 1#) (goJ i vi wi from_ s)

      goJ i vi wi j s
        | isTrue# (j ># i) = s
        | otherwise =
            let vj = indexDoubleArray# ba_v (off_v +# j)
                wj = indexDoubleArray# ba_w (off_w +# j)
                delta = vi *## wj +## wi *## vj
                ij = off +# i *# n +# j
                ji = off +# j *# n +# i
            in case readDoubleArray# mba ij s of
                 (# s1, tij #) ->
                   let tij' = tij -## delta
                   in case writeDoubleArray# mba ij tij' s1 of
                        s2 | isTrue# (i ==# j) -> goJ i vi wi (j +# 1#) s2
                           | otherwise ->
                               case writeDoubleArray# mba ji tij' s2 of
                                 s3 -> goJ i vi wi (j +# 1#) s3

  in (# goI from_ s0, () #)
{-# INLINE rawSymRank2Update #-}

-- --------------------------------------------------------------------------
-- Tridiagonalisation mutable kernels
-- --------------------------------------------------------------------------

-- | Symmetric submatrix-vector product for tridiagonalisation.
-- Computes p[i-from] = Σ_{j=from}^{to-1} T[i,j] * v[j-from]
-- for i in [from..to-1].
-- T is read from MutableByteArray (being modified in-place),
-- v is read from MutableByteArray (temporary vector),
-- p is written to MutableByteArray (temporary vector).
rawMutSymMatvecSub :: MutableByteArray s -> Int -> Int
                   -> MutableByteArray s -> Int
                   -> MutableByteArray s -> Int
                   -> Int -> Int -> ST s ()
rawMutSymMatvecSub (MutableByteArray mba_t) (I# off_t) (I# ncols)
                   (MutableByteArray mba_v) (I# off_v)
                   (MutableByteArray mba_p) (I# off_p)
                   (I# from_) (I# to) = ST $ \s0 ->
  let goI i s
        | isTrue# (i >=# to) = s
        | otherwise = case goJ i from_ 0.0## s of
            (# s', acc #) ->
              case writeDoubleArray# mba_p (off_p +# (i -# from_)) acc s' of
                s'' -> goI (i +# 1#) s''

      goJ i j acc s
        | isTrue# (j >=# to) = (# s, acc #)
        | otherwise =
            case readDoubleArray# mba_t (off_t +# i *# ncols +# j) s of
              (# s1, tij #) ->
                case readDoubleArray# mba_v (off_v +# (j -# from_)) s1 of
                  (# s2, vj #) -> goJ i (j +# 1#) (acc +## tij *## vj) s2

  in (# goI from_ s0, () #)
{-# INLINE rawMutSymMatvecSub #-}

-- | Symmetric rank-2 update reading v, w from MutableByteArrays.
-- For i in [from..to-1], j in [from..i]:
--   T[i,j] -= v[i-from]*w[j-from] + w[i-from]*v[j-from]
--   T[j,i] = T[i,j]   (maintain symmetry)
-- v and w are indexed relative to from (i.e., v[0] corresponds to row 'from').
rawMutSymRank2Update :: MutableByteArray s -> Int -> Int
                     -> MutableByteArray s -> Int
                     -> MutableByteArray s -> Int
                     -> Int -> Int -> ST s ()
rawMutSymRank2Update (MutableByteArray mba_t) (I# off_t) (I# n)
                     (MutableByteArray mba_v) (I# off_v)
                     (MutableByteArray mba_w) (I# off_w)
                     (I# from_) (I# to) = ST $ \s0 ->
  let goI i s
        | isTrue# (i >=# to) = s
        | otherwise =
            case readDoubleArray# mba_v (off_v +# (i -# from_)) s of
              (# s1, vi #) ->
                case readDoubleArray# mba_w (off_w +# (i -# from_)) s1 of
                  (# s2, wi #) -> goI (i +# 1#) (goJ i vi wi from_ s2)

      goJ i vi wi j s
        | isTrue# (j ># i) = s
        | otherwise =
            case readDoubleArray# mba_v (off_v +# (j -# from_)) s of
              (# s1, vj #) ->
                case readDoubleArray# mba_w (off_w +# (j -# from_)) s1 of
                  (# s2, wj #) ->
                    let delta = vi *## wj +## wi *## vj
                        ij = off_t +# i *# n +# j
                        ji = off_t +# j *# n +# i
                    in case readDoubleArray# mba_t ij s2 of
                         (# s3, tij #) ->
                           let tij' = tij -## delta
                           in case writeDoubleArray# mba_t ij tij' s3 of
                                s4 | isTrue# (i ==# j) -> goJ i vi wi (j +# 1#) s4
                                   | otherwise ->
                                       case writeDoubleArray# mba_t ji tij' s4 of
                                         s5 -> goJ i vi wi (j +# 1#) s5

  in (# goI from_ s0, () #)
{-# INLINE rawMutSymRank2Update #-}

-- | Q accumulation for tridiagonalisation.
-- Householder vectors are stored in column hvCol of frozen T,
-- with implicit v[qCol] = 1.0.
-- Phase 1: wi = beta * (Q[row,qCol] + Σ_{l=qCol+1}^{endRow-1} Q[row,l] * T[l,hvCol])
-- Phase 2: Q[row,qCol] -= wi; Q[row,l] -= wi * T[l,hvCol]
rawMutTridiagQAccum :: MutableByteArray s -> Int -> Int
                    -> ByteArray -> Int -> Int
                    -> Double -> Int -> Int -> Int -> Int -> ST s ()
rawMutTridiagQAccum (MutableByteArray mba_q) (I# off_q) (I# qcols)
                    (ByteArray ba_t) (I# off_t) (I# tcols)
                    (D# beta) (I# qCol) (I# hvCol) (I# endRow) (I# row) = ST $ \s0 ->
  -- Phase 1: wi = beta * (Q[row,qCol] + Σ Q[row,l] * T[l,hvCol])
  case readDoubleArray# mba_q (off_q +# row *# qcols +# qCol) s0 of
    (# s1, qrk #) ->
      let goSum l acc s
            | isTrue# (l >=# endRow) = (# s, beta *## (qrk +## acc) #)
            | otherwise =
                let vl = indexDoubleArray# ba_t (off_t +# l *# tcols +# hvCol)
                in case readDoubleArray# mba_q (off_q +# row *# qcols +# l) s of
                     (# s', qrl #) -> goSum (l +# 1#) (acc +## qrl *## vl) s'
      in case goSum (qCol +# 1#) 0.0## s1 of
           (# s2, wi #) ->
             -- Phase 2: Q[row,qCol] -= wi
             case writeDoubleArray# mba_q (off_q +# row *# qcols +# qCol) (qrk -## wi) s2 of
               s3 ->
                 let goUpdate l s
                       | isTrue# (l >=# endRow) = s
                       | otherwise =
                           let vl = indexDoubleArray# ba_t (off_t +# l *# tcols +# hvCol)
                           in case readDoubleArray# mba_q (off_q +# row *# qcols +# l) s of
                                (# s', qrl #) ->
                                  case writeDoubleArray# mba_q (off_q +# row *# qcols +# l) (qrl -## wi *## vl) s' of
                                    s'' -> goUpdate (l +# 1#) s''
                 in (# goUpdate (qCol +# 1#) s3, () #)
{-# INLINE rawMutTridiagQAccum #-}

-- --------------------------------------------------------------------------
-- LU kernels
-- --------------------------------------------------------------------------

-- | In-place LU elimination for column k of an n×n row-major matrix.
-- Computes multipliers and updates the trailing submatrix.
-- Inner j-loop uses DoubleX4# SIMD (contiguous row access).
rawLUEliminateColumn :: MutableByteArray s -> Int -> Int -> Int -> ST s ()
rawLUEliminateColumn (MutableByteArray mba) (I# off) (I# n) (I# k) = ST $ \s0 ->
  -- Read A[k,k] (the pivot)
  case readDoubleArray# mba (off +# k *# n +# k) s0 of
    (# s1, akk #) ->
      let goI i s
            | isTrue# (i >=# n) = s
            | otherwise =
                -- Read A[i,k], compute multiplier
                case readDoubleArray# mba (off +# i *# n +# k) s of
                  (# s', aik #) ->
                    let mult = aik /## akk
                        iRowOff = off +# i *# n
                        kRowOff = off +# k *# n
                        jSpan = n -# k -# 1#
                        jStart = k +# 1#
                        j4End = jStart +# (jSpan -# (jSpan `remInt#` 4#))
                        negMultV = broadcastDoubleX4# (negateDouble# mult)
                    -- Store multiplier at A[i,k]
                    in case writeDoubleArray# mba (off +# i *# n +# k) mult s' of
                         s'' ->
                           -- SIMD j-loop: A[i,j] -= mult * A[k,j]  =  A[i,j] + (-mult)*A[k,j]
                           let goJSimd j s_
                                 | isTrue# (j >=# j4End) = s_
                                 | otherwise =
                                     case readDoubleArrayAsDoubleX4# mba (iRowOff +# j) s_ of
                                       (# s1_, aij #) ->
                                         case readDoubleArrayAsDoubleX4# mba (kRowOff +# j) s1_ of
                                              (# s2_, akjV_ #) ->
                                                let aij' = fmaddDoubleX4# negMultV akjV_ aij
                                                in case writeDoubleArrayAsDoubleX4# mba (iRowOff +# j) aij' s2_ of
                                                     s3_ -> goJSimd (j +# 4#) s3_
                               -- Scalar cleanup
                               goJScalar j s_
                                 | isTrue# (j >=# n) = s_
                                 | otherwise =
                                     case readDoubleArray# mba (iRowOff +# j) s_ of
                                       (# s1_, aij #) ->
                                         case readDoubleArray# mba (kRowOff +# j) s1_ of
                                           (# s2_, akj #) ->
                                             case writeDoubleArray# mba (iRowOff +# j) (aij -## mult *## akj) s2_ of
                                               s3_ -> goJScalar (j +# 1#) s3_
                           in goI (i +# 1#) (goJScalar j4End (goJSimd jStart s''))
      in (# goI (k +# 1#) s1, () #)
{-# INLINE rawLUEliminateColumn #-}

-- | Swap elements in columns [fromCol..n-1] between two rows of an n-wide matrix.
-- Uses DoubleX4# SIMD for the contiguous row data.
rawSwapRows :: MutableByteArray s -> Int -> Int -> Int -> Int -> Int -> ST s ()
rawSwapRows (MutableByteArray mba) (I# off) (I# n) (I# row1) (I# row2) (I# fromCol) = ST $ \s0 ->
  let r1Off = off +# row1 *# n
      r2Off = off +# row2 *# n
      jSpan = n -# fromCol
      j4End = fromCol +# (jSpan -# (jSpan `remInt#` 4#))

      goSimd j s
        | isTrue# (j >=# j4End) = s
        | otherwise =
            case readDoubleArrayAsDoubleX4# mba (r1Off +# j) s of
              (# s1, v1 #) ->
                case readDoubleArrayAsDoubleX4# mba (r2Off +# j) s1 of
                  (# s2, v2 #) ->
                    case writeDoubleArrayAsDoubleX4# mba (r1Off +# j) v2 s2 of
                      s3 -> case writeDoubleArrayAsDoubleX4# mba (r2Off +# j) v1 s3 of
                              s4 -> goSimd (j +# 4#) s4

      goScalar j s
        | isTrue# (j >=# n) = s
        | otherwise =
            case readDoubleArray# mba (r1Off +# j) s of
              (# s1, v1 #) ->
                case readDoubleArray# mba (r2Off +# j) s1 of
                  (# s2, v2 #) ->
                    case writeDoubleArray# mba (r1Off +# j) v2 s2 of
                      s3 -> case writeDoubleArray# mba (r2Off +# j) v1 s3 of
                              s4 -> goScalar (j +# 1#) s4

  in (# goScalar j4End (goSimd fromCol s0), () #)
{-# INLINE rawSwapRows #-}

-- | Find row with maximum |A[i,k]| for i in [fromRow..n-1].
-- Returns the row index of the pivot.
rawPivotSearch :: MutableByteArray s -> Int -> Int -> Int -> Int -> ST s Int
rawPivotSearch (MutableByteArray mba) (I# off) (I# n) (I# k) (I# fromRow) = ST $ \s0 ->
  let go i bestIdx bestVal s
        | isTrue# (i >=# n) = (# s, I# bestIdx #)
        | otherwise =
            case readDoubleArray# mba (off +# i *# n +# k) s of
              (# s', v #) ->
                let av = if isTrue# (v >=## 0.0##) then v else negateDouble# v
                in if isTrue# (av >## bestVal)
                   then go (i +# 1#) i av s'
                   else go (i +# 1#) bestIdx bestVal s'
  in case readDoubleArray# mba (off +# fromRow *# n +# k) s0 of
       (# s1, v0 #) ->
         let av0 = if isTrue# (v0 >=## 0.0##) then v0 else negateDouble# v0
         in go (fromRow +# 1#) fromRow av0 s1
{-# INLINE rawPivotSearch #-}

-- | In-place forward substitution using packed LU matrix (unit lower triangular).
-- Solves Ly = b where L is stored in the strictly lower part of ba_lu.
-- The solution overwrites mba_x.
rawForwardSubUnitPacked :: ByteArray -> Int -> Int -> MutableByteArray s -> Int -> ST s ()
rawForwardSubUnitPacked (ByteArray ba_lu) (I# off_lu) (I# n)
                        (MutableByteArray mba_x) (I# off_x) = ST $ \s0 ->
  let goJ j s
        | isTrue# (j >=# n) = s
        | otherwise =
            case readDoubleArray# mba_x (off_x +# j) s of
              (# s', xj #) ->
                let goI i s_
                      | isTrue# (i >=# n) = s_
                      | otherwise =
                          let lij = indexDoubleArray# ba_lu (off_lu +# i *# n +# j)
                          in case readDoubleArray# mba_x (off_x +# i) s_ of
                               (# s1, xi #) ->
                                 case writeDoubleArray# mba_x (off_x +# i) (xi -## lij *## xj) s1 of
                                   s2 -> goI (i +# 1#) s2
                in goJ (j +# 1#) (goI (j +# 1#) s')
  in (# goJ 0# s0, () #)
{-# INLINE rawForwardSubUnitPacked #-}

-- | In-place back substitution using packed LU matrix (upper triangular).
-- Solves Ux = y where U is stored in the upper part of ba_lu.
-- The solution overwrites mba_x.
rawBackSubPacked :: ByteArray -> Int -> Int -> MutableByteArray s -> Int -> ST s ()
rawBackSubPacked (ByteArray ba_lu) (I# off_lu) (I# n)
                 (MutableByteArray mba_x) (I# off_x) = ST $ \s0 ->
  let goJ j s
        | isTrue# (j <# 0#) = s
        | otherwise =
            -- x[j] /= U[j,j]
            let ujj = indexDoubleArray# ba_lu (off_lu +# j *# n +# j)
            in case readDoubleArray# mba_x (off_x +# j) s of
                 (# s', xj_ #) ->
                   let xj = xj_ /## ujj
                   in case writeDoubleArray# mba_x (off_x +# j) xj s' of
                        s'' ->
                          -- for i = 0..j-1: x[i] -= U[i,j] * x[j]
                          let goI i s_
                                | isTrue# (i >=# j) = s_
                                | otherwise =
                                    let uij = indexDoubleArray# ba_lu (off_lu +# i *# n +# j)
                                    in case readDoubleArray# mba_x (off_x +# i) s_ of
                                         (# s1, xi #) ->
                                           case writeDoubleArray# mba_x (off_x +# i) (xi -## uij *## xj) s1 of
                                             s2 -> goI (i +# 1#) s2
                          in goJ (j -# 1#) (goI 0# s'')
  in (# goJ (n -# 1#) s0, () #)
{-# INLINE rawBackSubPacked #-}

-- --------------------------------------------------------------------------
-- Cholesky kernels
-- --------------------------------------------------------------------------

-- | Process one column j of Cholesky factorisation in-place.
-- For k in [0..j-1]: subtract G[i,k]*G[j,k] from G[i,j] for i in [j..n-1].
-- Then scale: G[j,j] = sqrt(G[j,j]); G[i,j] /= G[j,j] for i > j.
rawCholColumn :: MutableByteArray s -> Int -> Int -> Int -> ST s ()
rawCholColumn (MutableByteArray mba) (I# off) (I# n) (I# j) = ST $ \s0 ->
  -- Phase 1: subtract contributions from previous columns
  let goK k s
        | isTrue# (k >=# j) = s
        | otherwise =
            -- Read G[j,k]
            case readDoubleArray# mba (off +# j *# n +# k) s of
              (# s', gjk #) ->
                -- For i in [j..n-1]: G[i,j] -= G[i,k] * gjk
                let goI i s_
                      | isTrue# (i >=# n) = s_
                      | otherwise =
                          case readDoubleArray# mba (off +# i *# n +# j) s_ of
                            (# s1, gij #) ->
                              case readDoubleArray# mba (off +# i *# n +# k) s1 of
                                (# s2, gik #) ->
                                  case writeDoubleArray# mba (off +# i *# n +# j) (gij -## gik *## gjk) s2 of
                                    s3 -> goI (i +# 1#) s3
                in goK (k +# 1#) (goI j s')
  in case goK 0# s0 of
       s1 ->
         -- Phase 2: scale column
         case readDoubleArray# mba (off +# j *# n +# j) s1 of
           (# s2, gjj #) ->
             let sjj = sqrtDouble# gjj
             in case writeDoubleArray# mba (off +# j *# n +# j) sjj s2 of
                  s3 ->
                    let goScale i s_
                          | isTrue# (i >=# n) = s_
                          | otherwise =
                              case readDoubleArray# mba (off +# i *# n +# j) s_ of
                                (# s4, gij #) ->
                                  case writeDoubleArray# mba (off +# i *# n +# j) (gij /## sjj) s4 of
                                    s5 -> goScale (i +# 1#) s5
                    in (# goScale (j +# 1#) s3, () #)
{-# INLINE rawCholColumn #-}

-- | Forward substitution with Cholesky factor G (lower triangular, non-unit diagonal).
-- Solves Gy = b, overwrites mba_x with y.
rawForwardSubCholPacked :: ByteArray -> Int -> Int -> MutableByteArray s -> Int -> ST s ()
rawForwardSubCholPacked (ByteArray ba_g) (I# off_g) (I# n)
                        (MutableByteArray mba_x) (I# off_x) = ST $ \s0 ->
  let goJ j s
        | isTrue# (j >=# n) = s
        | otherwise =
            let gjj = indexDoubleArray# ba_g (off_g +# j *# n +# j)
            in case readDoubleArray# mba_x (off_x +# j) s of
                 (# s', xj_ #) ->
                   let xj = xj_ /## gjj
                   in case writeDoubleArray# mba_x (off_x +# j) xj s' of
                        s'' ->
                          let goI i s_
                                | isTrue# (i >=# n) = s_
                                | otherwise =
                                    let gij = indexDoubleArray# ba_g (off_g +# i *# n +# j)
                                    in case readDoubleArray# mba_x (off_x +# i) s_ of
                                         (# s1, xi #) ->
                                           case writeDoubleArray# mba_x (off_x +# i) (xi -## gij *## xj) s1 of
                                             s2 -> goI (i +# 1#) s2
                          in goJ (j +# 1#) (goI (j +# 1#) s'')
  in (# goJ 0# s0, () #)
{-# INLINE rawForwardSubCholPacked #-}

-- | Back substitution with G^T (upper triangular) WITHOUT forming G^T.
-- Solves G^T x = y, overwrites mba_x with x.
-- Uses G^T[i,j] = G[j,i] to read from the lower triangle.
rawBackSubCholTPacked :: ByteArray -> Int -> Int -> MutableByteArray s -> Int -> ST s ()
rawBackSubCholTPacked (ByteArray ba_g) (I# off_g) (I# n)
                      (MutableByteArray mba_x) (I# off_x) = ST $ \s0 ->
  let goJ j s
        | isTrue# (j <# 0#) = s
        | otherwise =
            -- G^T[j,j] = G[j,j]
            let gjj = indexDoubleArray# ba_g (off_g +# j *# n +# j)
            in case readDoubleArray# mba_x (off_x +# j) s of
                 (# s', xj_ #) ->
                   let xj = xj_ /## gjj
                   in case writeDoubleArray# mba_x (off_x +# j) xj s' of
                        s'' ->
                          -- for i = 0..j-1: x[i] -= G^T[i,j] * x[j] = G[j,i] * x[j]
                          let goI i s_
                                | isTrue# (i >=# j) = s_
                                | otherwise =
                                    -- G^T[i,j] = G[j,i]
                                    let gji = indexDoubleArray# ba_g (off_g +# j *# n +# i)
                                    in case readDoubleArray# mba_x (off_x +# i) s_ of
                                         (# s1, xi #) ->
                                           case writeDoubleArray# mba_x (off_x +# i) (xi -## gji *## xj) s1 of
                                             s2 -> goI (i +# 1#) s2
                          in goJ (j -# 1#) (goI 0# s'')
  in (# goJ (n -# 1#) s0, () #)
{-# INLINE rawBackSubCholTPacked #-}

-- --------------------------------------------------------------------------
-- QR mutable kernels
-- --------------------------------------------------------------------------

-- | Sum of squares of a column slice in a mutable row-major matrix.
-- Σ A[i,col]² for i in [startRow..endRow-1].
rawMutSumSqColumn :: MutableByteArray s -> Int -> Int -> Int -> Int -> Int -> ST s Double
rawMutSumSqColumn (MutableByteArray mba) (I# off) (I# ncols) (I# startRow) (I# endRow) (I# col) = ST $ \s0 ->
  let go i acc s
        | isTrue# (i >=# endRow) = (# s, D# acc #)
        | otherwise =
            case readDoubleArray# mba (off +# i *# ncols +# col) s of
              (# s', v #) -> go (i +# 1#) (acc +## v *## v) s'
  in go startRow 0.0## s0
{-# INLINE rawMutSumSqColumn #-}

-- | Dot product of two column slices in a mutable row-major matrix.
-- Σ A[i,col1] * A[i,col2] for i in [startRow..endRow-1].
rawMutSumProdColumns :: MutableByteArray s -> Int -> Int -> Int -> Int -> Int -> Int -> ST s Double
rawMutSumProdColumns (MutableByteArray mba) (I# off) (I# ncols) (I# startRow) (I# endRow) (I# col1) (I# col2) = ST $ \s0 ->
  let go i acc s
        | isTrue# (i >=# endRow) = (# s, D# acc #)
        | otherwise =
            case readDoubleArray# mba (off +# i *# ncols +# col1) s of
              (# s1, v1 #) ->
                case readDoubleArray# mba (off +# i *# ncols +# col2) s1 of
                  (# s2, v2 #) -> go (i +# 1#) (acc +## v1 *## v2) s2
  in go startRow 0.0## s0
{-# INLINE rawMutSumProdColumns #-}

-- | Apply Householder reflector stored in column k (rows k+1..endRow-1,
-- with v[k]=1 implicit) to targetCol of a mutable row-major matrix.
-- Phase 1: w = beta * (R[k,targetCol] + Σ_{i=k+1}^{endRow-1} v[i]*R[i,targetCol])
-- Phase 2: R[k,targetCol] -= w; R[i,targetCol] -= v[i]*w
rawMutHouseholderApply :: MutableByteArray s -> Int -> Int -> Double
                       -> Int -> Int -> Int -> ST s ()
rawMutHouseholderApply (MutableByteArray mba) (I# off) (I# ncols) (D# beta)
                       (I# k) (I# endRow) (I# targetCol) = ST $ \s0 ->
  -- Phase 1: compute dot product
  case readDoubleArray# mba (off +# k *# ncols +# targetCol) s0 of
    (# s1, rkj #) ->
      let goSum i acc s
            | isTrue# (i >=# endRow) = (# s, beta *## (rkj +## acc) #)
            | otherwise =
                case readDoubleArray# mba (off +# i *# ncols +# k) s of
                  (# s', vi #) ->
                    case readDoubleArray# mba (off +# i *# ncols +# targetCol) s' of
                      (# s'', rij #) -> goSum (i +# 1#) (acc +## vi *## rij) s''
      in case goSum (k +# 1#) 0.0## s1 of
           (# s2, w #) ->
             -- Phase 2: update R[k,targetCol]
             case writeDoubleArray# mba (off +# k *# ncols +# targetCol) (rkj -## w) s2 of
               s3 ->
                 let goUpdate i s
                       | isTrue# (i >=# endRow) = s
                       | otherwise =
                           case readDoubleArray# mba (off +# i *# ncols +# k) s of
                             (# s', vi #) ->
                               case readDoubleArray# mba (off +# i *# ncols +# targetCol) s' of
                                 (# s'', rij #) ->
                                   case writeDoubleArray# mba (off +# i *# ncols +# targetCol) (rij -## vi *## w) s'' of
                                     s''' -> goUpdate (i +# 1#) s'''
                 in (# goUpdate (k +# 1#) s3, () #)
{-# INLINE rawMutHouseholderApply #-}

-- | Apply stored Householder reflector to one row of Q during accumulation.
-- v is stored in the subdiagonal of frozen R (column k, rows k+1..endRow-1).
-- Phase 1: wi = beta * (Q[row,k] + Σ_{l=k+1}^{endRow-1} Q[row,l] * v[l])
-- Phase 2: Q[row,k] -= wi; Q[row,l] -= wi * v[l]
rawMutQAccum :: MutableByteArray s -> Int -> Int
             -> ByteArray -> Int -> Int
             -> Double -> Int -> Int -> Int -> ST s ()
rawMutQAccum (MutableByteArray mba_q) (I# off_q) (I# qcols)
             (ByteArray ba_r) (I# off_r) (I# rcols)
             (D# beta) (I# k) (I# endRow) (I# row) = ST $ \s0 ->
  -- Phase 1: compute wi = beta * (Q[row,k] + Σ Q[row,l] * v[l])
  case readDoubleArray# mba_q (off_q +# row *# qcols +# k) s0 of
    (# s1, qrk #) ->
      let goSum l acc s
            | isTrue# (l >=# endRow) = (# s, beta *## (qrk +## acc) #)
            | otherwise =
                let vl = indexDoubleArray# ba_r (off_r +# l *# rcols +# k)
                in case readDoubleArray# mba_q (off_q +# row *# qcols +# l) s of
                     (# s', qrl #) -> goSum (l +# 1#) (acc +## qrl *## vl) s'
      in case goSum (k +# 1#) 0.0## s1 of
           (# s2, wi #) ->
             -- Phase 2: Q[row,k] -= wi
             case writeDoubleArray# mba_q (off_q +# row *# qcols +# k) (qrk -## wi) s2 of
               s3 ->
                 let goUpdate l s
                       | isTrue# (l >=# endRow) = s
                       | otherwise =
                           let vl = indexDoubleArray# ba_r (off_r +# l *# rcols +# k)
                           in case readDoubleArray# mba_q (off_q +# row *# qcols +# l) s of
                                (# s', qrl #) ->
                                  case writeDoubleArray# mba_q (off_q +# row *# qcols +# l) (qrl -## wi *## vl) s' of
                                    s'' -> goUpdate (l +# 1#) s''
                 in (# goUpdate (k +# 1#) s3, () #)
{-# INLINE rawMutQAccum #-}

-- --------------------------------------------------------------------------
-- Eigen mutable kernels
-- --------------------------------------------------------------------------

-- | Apply Givens rotation to two columns of a mutable matrix.
-- For each row in [0..nrows-1]:
--   tmp = c * M[row,col_p] + s * M[row,col_q]
--   M[row,col_q] = -s * M[row,col_p] + c * M[row,col_q]
--   M[row,col_p] = tmp
rawMutApplyGivensColumns :: MutableByteArray s -> Int -> Int
                         -> Double -> Double -> Int -> Int -> Int -> ST s ()
rawMutApplyGivensColumns (MutableByteArray mba) (I# off) (I# ncols)
                         (D# c_) (D# s_) (I# col_p) (I# col_q) (I# nrows) = ST $ \s0 ->
  let go row s
        | isTrue# (row >=# nrows) = s
        | otherwise =
            let pIdx = off +# row *# ncols +# col_p
                qIdx = off +# row *# ncols +# col_q
            in case readDoubleArray# mba pIdx s of
                 (# s1, mp #) ->
                   case readDoubleArray# mba qIdx s1 of
                     (# s2, mq #) ->
                       let tmp = c_ *## mp +## s_ *## mq
                           qnew = negateDouble# s_ *## mp +## c_ *## mq
                       in case writeDoubleArray# mba pIdx tmp s2 of
                            s3 -> case writeDoubleArray# mba qIdx qnew s3 of
                                    s4 -> go (row +# 1#) s4
  in (# go 0# s0, () #)
{-# INLINE rawMutApplyGivensColumns #-}

-- --------------------------------------------------------------------------
-- SVD / bidiagonalisation kernels
-- --------------------------------------------------------------------------

-- | Apply a right Householder reflector to one row of a mutable matrix.
-- The Householder vector v is stored in row hvRow of the matrix,
-- columns [hvStart..hvEnd-1], with implicit v[hvStart] = 1.0.
-- Updates row targetRow: R[targetRow, hvStart..hvEnd-1] -= w * v
-- where w = beta * (R[targetRow,hvStart] + Σ_{l=hvStart+1}^{hvEnd-1} R[targetRow,l] * R[hvRow,l])
rawMutHouseholderApplyRow :: MutableByteArray s -> Int -> Int
                          -> Double -> Int -> Int -> Int -> Int -> ST s ()
rawMutHouseholderApplyRow (MutableByteArray mba) (I# off) (I# ncols) (D# beta)
                          (I# hvRow) (I# hvStart) (I# hvEnd) (I# targetRow) = ST $ \s0 ->
  let trOff = off +# targetRow *# ncols
      hvOff = off +# hvRow *# ncols
  -- Phase 1: w = beta * (R[targetRow,hvStart] + Σ R[targetRow,l] * R[hvRow,l])
  in case readDoubleArray# mba (trOff +# hvStart) s0 of
       (# s1, r0 #) ->
         let goSum l acc s
               | isTrue# (l >=# hvEnd) = (# s, beta *## (r0 +## acc) #)
               | otherwise =
                   case readDoubleArray# mba (trOff +# l) s of
                     (# s', rl #) ->
                       case readDoubleArray# mba (hvOff +# l) s' of
                         (# s'', vl #) -> goSum (l +# 1#) (acc +## rl *## vl) s''
         in case goSum (hvStart +# 1#) 0.0## s1 of
              (# s2, w #) ->
                -- Phase 2: R[targetRow,hvStart] -= w (implicit v[hvStart]=1)
                case writeDoubleArray# mba (trOff +# hvStart) (r0 -## w) s2 of
                  s3 ->
                    let goUpdate l s
                          | isTrue# (l >=# hvEnd) = s
                          | otherwise =
                              case readDoubleArray# mba (hvOff +# l) s of
                                (# s', vl #) ->
                                  case readDoubleArray# mba (trOff +# l) s' of
                                    (# s'', rl #) ->
                                      case writeDoubleArray# mba (trOff +# l) (rl -## w *## vl) s'' of
                                        s''' -> goUpdate (l +# 1#) s'''
                    in (# goUpdate (hvStart +# 1#) s3, () #)
{-# INLINE rawMutHouseholderApplyRow #-}

-- | Sum of squares of a row slice in a mutable row-major matrix.
-- Σ A[row,j]² for j in [startCol..endCol-1].
rawMutSumSqRow :: MutableByteArray s -> Int -> Int -> Int -> Int -> Int -> ST s Double
rawMutSumSqRow (MutableByteArray mba) (I# off) (I# ncols) (I# row) (I# startCol) (I# endCol) = ST $ \s0 ->
  let rowOff = off +# row *# ncols
      go j acc s
        | isTrue# (j >=# endCol) = (# s, D# acc #)
        | otherwise =
            case readDoubleArray# mba (rowOff +# j) s of
              (# s', v #) -> go (j +# 1#) (acc +## v *## v) s'
  in go startCol 0.0## s0
{-# INLINE rawMutSumSqRow #-}

-- --------------------------------------------------------------------------
-- Cholesky SIMD kernel
-- --------------------------------------------------------------------------

-- | SIMD-vectorised Cholesky column kernel.
-- Restructures the inner loop as a dot product of contiguous row segments:
--   G[i,j] -= Σ_{k=0}^{j-1} G[i,k] * G[j,k]
-- which is a dot product of row[i][0..j-1] and row[j][0..j-1].
-- Since rows are contiguous in row-major storage, this enables DoubleX4# SIMD.
rawCholColumnSIMD :: MutableByteArray s -> Int -> Int -> Int -> ST s ()
rawCholColumnSIMD (MutableByteArray mba) (I# off) (I# n) (I# j) = ST $ \s0 ->
  let jRowOff = off +# j *# n
      -- For each row i in [j..n-1], subtract dot(row[i][0..j-1], row[j][0..j-1])
      j4 = j -# (j `remInt#` 4#)  -- SIMD boundary for dot of length j

      goI i s
        | isTrue# (i >=# n) = s
        | otherwise =
            let iRowOff = off +# i *# n
            in case mutRowDot iRowOff jRowOff 0# j4 j s of
                 (# s', dot #) ->
                   case readDoubleArray# mba (iRowOff +# j) s' of
                     (# s'', gij #) ->
                       case writeDoubleArray# mba (iRowOff +# j) (gij -## dot) s'' of
                         s''' -> goI (i +# 1#) s'''

      -- Dot product of two mutable row segments using SIMD
      mutRowDot r1 r2 k k4End kEnd s
        -- SIMD phase
        | isTrue# (k <# k4End) =
            case goSimd r1 r2 k k4End (broadcastDoubleX4# 0.0##) s of
              (# s', acc4 #) ->
                let !(# a, b, c, d #) = unpackDoubleX4# acc4
                    simdSum = a +## b +## c +## d
                in mutRowDotScalar r1 r2 k4End kEnd simdSum s'
        | otherwise = mutRowDotScalar r1 r2 k kEnd 0.0## s

      goSimd r1 r2 k k4End acc s
        | isTrue# (k >=# k4End) = (# s, acc #)
        | otherwise =
            case readDoubleArrayAsDoubleX4# mba (r1 +# k) s of
              (# s1, v1 #) ->
                case readDoubleArrayAsDoubleX4# mba (r2 +# k) s1 of
                  (# s2, v2 #) -> goSimd r1 r2 (k +# 4#) k4End (fmaddDoubleX4# v1 v2 acc) s2

      mutRowDotScalar r1 r2 k kEnd acc s
        | isTrue# (k >=# kEnd) = (# s, acc #)
        | otherwise =
            case readDoubleArray# mba (r1 +# k) s of
              (# s1, v1 #) ->
                case readDoubleArray# mba (r2 +# k) s1 of
                  (# s2, v2 #) -> mutRowDotScalar r1 r2 (k +# 1#) kEnd (acc +## v1 *## v2) s2

  in case goI j s0 of
       s1 ->
         -- Scale column: G[j,j] = sqrt(G[j,j]); G[i,j] /= G[j,j] for i > j
         case readDoubleArray# mba (jRowOff +# j) s1 of
           (# s2, gjj #) ->
             let sjj = sqrtDouble# gjj
             in case writeDoubleArray# mba (jRowOff +# j) sjj s2 of
                  s3 ->
                    let goScale i s
                          | isTrue# (i >=# n) = s
                          | otherwise =
                              case readDoubleArray# mba (off +# i *# n +# j) s of
                                (# s4, gij #) ->
                                  case writeDoubleArray# mba (off +# i *# n +# j) (gij /## sjj) s4 of
                                    s5 -> goScale (i +# 1#) s5
                    in (# goScale (j +# 1#) s3, () #)
{-# INLINE rawCholColumnSIMD #-}

-- --------------------------------------------------------------------------
-- Utilities
-- --------------------------------------------------------------------------

minI :: Int# -> Int# -> Int#
minI a b = if isTrue# (a <=# b) then a else b
{-# INLINE minI #-}
