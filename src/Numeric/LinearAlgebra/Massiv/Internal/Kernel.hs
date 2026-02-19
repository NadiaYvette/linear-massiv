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
    -- * QR helpers
  , rawSumSqRange
  , rawSumProdRange
  , rawHouseholderApplyCol
  , rawQAccumCol
    -- * Eigen / tridiag helpers
  , rawApplyGivensRows
  , rawSymRank2Update
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
rawGemmKernel (ByteArray ba_a) (I# off_a) (ByteArray ba_b) (I# off_b)
              (MutableByteArray mba_c) (I# off_c)
              (I# m) (I# k) (I# n) = ST $ \s0 ->
  let bs = 64#

      goBI bi s
        | isTrue# (bi >=# m) = s
        | otherwise =
            let iEnd = minI (bi +# bs) m
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

  in (# goBI 0# s0, () #)
{-# INLINE rawGemmKernel #-}

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
-- Utilities
-- --------------------------------------------------------------------------

minI :: Int# -> Int# -> Int#
minI a b = if isTrue# (a <=# b) then a else b
{-# INLINE minI #-}
