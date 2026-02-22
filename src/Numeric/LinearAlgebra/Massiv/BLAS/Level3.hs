{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE MagicHash #-}
{-# LANGUAGE UnboxedTuples #-}

-- |
-- Module      : Numeric.LinearAlgebra.Massiv.BLAS.Level3
-- Copyright   : (c) Nadia Chambers 2026
-- License     : BSD-3-Clause
-- Maintainer  : nadia.chambers@iohk.io
-- Stability   : experimental
--
-- = BLAS Level 3: Matrix–Matrix Operations
--
-- This module provides type-safe, dimension-indexed wrappers around the
-- standard BLAS Level 3 kernels.  These are /matrix–matrix/ operations
-- whose arithmetic cost is \(O(m \, n \, k)\) for an
-- \(m \times k\) by \(k \times n\) multiply, one level above the
-- \(O(m \, n)\) matrix–vector operations of Level 2.
--
-- The algorithms implemented here correspond to the six loop orderings
-- of the triple-loop matrix multiplication described in:
--
-- * Golub, G. H. & Van Loan, C. F. (2013). /Matrix Computations/,
--   4th edition (GVL4). Johns Hopkins University Press.
--   __Chapter 1, Sections 1.1.5–1.1.8__, pp. 12–18.
--
-- Specifically:
--
-- * __Algorithm 1.1.5__ (ijk variant, p. 12) — The "row-oriented
--   inner-product" form.  For each entry \(c_{ij}\) the inner product
--   of row \(i\) of \(A\) with column \(j\) of \(B\) is computed.
--   This is the variant implemented by 'gemm', 'matMul', and
--   'matMulComp'.
--
-- * __Algorithms 1.1.6–1.1.8__ (pp. 13–15) — The jki (Gaxpy),
--   kji (outer-product), and other orderings.  These alternatives
--   differ in data-access pattern but compute the same result.  The
--   present implementation uses the ijk ordering; cache-oblivious or
--   blocked variants can be added in the future.
--
-- The module also provides elementary matrix arithmetic (addition,
-- subtraction, scaling, transpose) and a triangular matrix–matrix
-- multiply ('trmmLeft') that exploits the triangular structure to
-- halve the work.
--
-- == Complexity
--
-- * 'gemm', 'matMul', 'matMulComp': \(O(m \, k \, n)\).
-- * 'transpose', 'mAdd', 'mSub', 'mScale': \(O(m \, n)\).
-- * 'trmmLeft': \(O(n^{3}/2)\) (triangular, in-place structure).
module Numeric.LinearAlgebra.Massiv.BLAS.Level3
  ( -- * Matrix–matrix multiply (Algorithm 1.1.5, GVL4 p. 12)
    gemm
  , matMul
  , matMulP
  , matMulPPar
  , matMulComp
    -- * Elementary matrix operations (GVL4 Section 1.1, pp. 4–18)
  , transpose
  , transposeP
  , matMulAtAP
  , mAdd
  , mSub
  , mScale
    -- * Triangular matrix multiply (GVL4 Section 1.1.8, p. 15)
  , trmmLeft
  ) where

import qualified Data.Massiv.Array as M
import Data.Massiv.Array (Ix2(..), Sz(..), Comp(..), unwrapByteArray, unwrapByteArrayOffset,
                          unwrapMutableByteArray, unwrapMutableByteArrayOffset)
import GHC.TypeNats (KnownNat)
import GHC.Exts (Int(..), isTrue#, (>=#), (*#), (+#))
import GHC.Prim
import GHC.ST (ST(..))
import Data.Array.Byte (MutableByteArray(..))
import Data.Primitive.ByteArray (ByteArray(..), newByteArray, unsafeFreezeByteArray)
import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Control.Concurrent (forkOn, newEmptyMVar, putMVar, takeMVar, getNumCapabilities)
import System.IO.Unsafe (unsafePerformIO)
import GHC.IO (stToIO)
import Numeric.LinearAlgebra.Massiv.Internal.Kernel (rawGemmKernel, rawGemmBISlice, rawSyrkLowerKernel)

-- | Block size for cache-tiled GEMM (generic fallback path).
gemmBlockSize :: Int
gemmBlockSize = 32
{-# INLINE gemmBlockSize #-}

-- | General matrix multiply (BLAS @GEMM@).
--
-- __GVL4 Reference:__ Algorithm 1.1.5 (ijk matrix multiply), p. 12.
--
-- Given \(A \in \mathbb{R}^{m \times k}\),
-- \(B \in \mathbb{R}^{k \times n}\),
-- \(C \in \mathbb{R}^{m \times n}\), and scalars
-- \(\alpha, \beta\), computes
--
-- \[
--   C \;\leftarrow\; \alpha \, A \, B \;+\; \beta \, C
-- \]
gemm :: forall m k n r e. (KnownNat m, KnownNat k, KnownNat n, M.Manifest r e, Num e)
     => e -> Matrix m k r e -> Matrix k n r e -> e -> Matrix m n r e -> Matrix m n r e
gemm alpha a b beta c =
  let mm = dimVal @m
      kk = dimVal @k
      nn = dimVal @n
      bs = gemmBlockSize
  in createMatrix @m @n @r $ \mc -> do
    -- Initialize C with β·C₀
    mapM_ (\i -> mapM_ (\j -> do
      M.write_ mc (i :. j) (beta * (c ! (i, j)))
      ) [0..nn-1]) [0..mm-1]
    -- Tiled ikj loop: for each block triple, accumulate α·A·B
    let go_bi bi = do
          let iEnd = min (bi + bs) mm
          let go_bk bk = do
                let kEnd = min (bk + bs) kk
                let go_bj bj = do
                      let jEnd = min (bj + bs) nn
                      -- Inner micro-kernel: ikj within the block
                      mapM_ (\i -> mapM_ (\l -> do
                        let aik = alpha * (a ! (i, l))
                        mapM_ (\j -> do
                          cij <- M.readM mc (i :. j)
                          M.write_ mc (i :. j) (cij + aik * (b ! (l, j)))
                          ) [bj..jEnd-1]
                        ) [bk..kEnd-1]) [bi..iEnd-1]
                mapM_ go_bj [0, bs .. nn-1]
          mapM_ go_bk [0, bs .. kk-1]
    mapM_ go_bi [0, bs .. mm-1]

-- | Simple matrix multiply (specialisation of 'gemm').
--
-- __GVL4 Reference:__ Algorithm 1.1.5 (ijk matrix multiply), p. 12,
-- with \(\alpha = 1\) and \(\beta = 0\).
--
-- Given \(A \in \mathbb{R}^{m \times k}\) and
-- \(B \in \mathbb{R}^{k \times n}\), computes
--
-- \[
--   C \;=\; A \, B
-- \]
matMul :: forall m k n r e. (KnownNat m, KnownNat k, KnownNat n, M.Manifest r e, Num e)
       => Matrix m k r e -> Matrix k n r e -> Matrix m n r e
matMul a b = matMulGeneric a b
{-# NOINLINE [1] matMul #-}

-- | Generic fallback for non-P or non-Double representations.
matMulGeneric :: forall m k n r e. (KnownNat m, KnownNat k, KnownNat n, M.Manifest r e, Num e)
              => Matrix m k r e -> Matrix k n r e -> Matrix m n r e
matMulGeneric a b =
  let mm = dimVal @m
      kk = dimVal @k
      nn = dimVal @n
      bs = gemmBlockSize
  in createMatrix @m @n @r $ \mc -> do
    -- Initialize C to zero
    mapM_ (\i -> mapM_ (\j ->
      M.write_ mc (i :. j) 0
      ) [0..nn-1]) [0..mm-1]
    -- Tiled ikj loop
    let go_bi bi = do
          let iEnd = min (bi + bs) mm
          let go_bk bk = do
                let kEnd = min (bk + bs) kk
                let go_bj bj = do
                      let jEnd = min (bj + bs) nn
                      mapM_ (\i -> mapM_ (\l -> do
                        let aik = a ! (i, l)
                        mapM_ (\j -> do
                          cij <- M.readM mc (i :. j)
                          M.write_ mc (i :. j) (cij + aik * (b ! (l, j)))
                          ) [bj..jEnd-1]
                        ) [bk..kEnd-1]) [bi..iEnd-1]
                mapM_ go_bj [0, bs .. nn-1]
          mapM_ go_bk [0, bs .. kk-1]
    mapM_ go_bi [0, bs .. mm-1]

-- | Specialised raw-array GEMM for @P Double@.
-- Bypasses massiv's per-element abstraction and uses raw ByteArray# primops
-- with AVX2 DoubleX4# SIMD for the inner kernel.
matMulP :: forall m k n. (KnownNat m, KnownNat k, KnownNat n)
        => Matrix m k M.P Double -> Matrix k n M.P Double -> Matrix m n M.P Double
matMulP (MkMatrix arrA) (MkMatrix arrB) =
  let mm = dimVal @m
      kk = dimVal @k
      nn = dimVal @n
      baA = unwrapByteArray arrA
      offA = unwrapByteArrayOffset arrA
      baB = unwrapByteArray arrB
      offB = unwrapByteArrayOffset arrB
  in createMatrix @m @n @M.P $ \mc -> do
    -- Zero-initialise C
    let mbaC = unwrapMutableByteArray mc
        offC = unwrapMutableByteArrayOffset mc
        !(I# mm#) = mm
        !(I# nn#) = nn
    ST $ \s0 ->
      let go i s
            | isTrue# (i >=# (mm# *# nn#)) = s
            | otherwise = case writeDoubleArray# (unMBA# mbaC) (unI offC +# i) 0.0## s of
                            s' -> go (i +# 1#) s'
      in (# go 0# s0, () #)
    -- Run the raw SIMD kernel
    rawGemmKernel baA offA baB offB mbaC offC mm kk nn
{-# NOINLINE matMulP #-}

-- | Parallel specialised GEMM for @P Double@.
-- Partitions the row range across OS threads, each calling 'rawGemmBISlice'
-- on non-overlapping row ranges.  Uses @forkIO@ + @MVar@ barrier.
-- Falls back to single-threaded 'matMulP' when @getNumCapabilities == 1@.
matMulPPar :: forall m k n. (KnownNat m, KnownNat k, KnownNat n)
           => Matrix m k M.P Double -> Matrix k n M.P Double -> Matrix m n M.P Double
matMulPPar a b = unsafePerformIO $ do
  let !mm = dimVal @m
      !kk = dimVal @k
      !nn = dimVal @n
      !baA = unwrapByteArray (unMatrix a)
      !offA = unwrapByteArrayOffset (unMatrix a)
      !baB = unwrapByteArray (unMatrix b)
      !offB = unwrapByteArrayOffset (unMatrix b)
  caps <- getNumCapabilities
  -- Adaptive thread count: ensure each thread gets enough rows to
  -- amortize fork/join overhead (minimum 16 rows per thread).
  let !minRowsPerThread = 16
      !maxThreads = max 1 (mm `div` minRowsPerThread)
      !numThreads = min caps (min mm maxThreads)
  if numThreads <= 1
    then pure (matMulP a b)
    else do
      -- Allocate mutable C, zero-initialise
      mc <- stToIO $ M.newMArray (Sz (mm :. nn)) (0 :: Double)
      let !mbaC = unwrapMutableByteArray mc
          !offC = unwrapMutableByteArrayOffset mc
          !chunkSize = (mm + numThreads - 1) `div` numThreads
      -- Fork threads pinned to capabilities (avoids OS migration)
      mvars <- mapM (\t -> do
        let !biStart = t * chunkSize
            !biEnd = min (biStart + chunkSize) mm
        mv <- newEmptyMVar
        _ <- forkOn t $ do
          stToIO $ rawGemmBISlice baA offA baB offB mbaC offC biStart biEnd mm kk nn
          putMVar mv ()
        pure mv
        ) [0..numThreads-1]
      -- Wait for all threads
      mapM_ takeMVar mvars
      -- Freeze and wrap
      arr <- stToIO $ M.freezeS mc
      pure (MkMatrix arr)
{-# NOINLINE matMulPPar #-}

{-# RULES
"matMul/P/Double" forall (a :: Matrix m k M.P Double)
                         (b :: Matrix k n M.P Double).
    matMul a b = matMulP a b
  #-}

-- | Matrix multiply with explicit computation strategy.
matMulComp :: forall m k n r e. (KnownNat m, KnownNat k, KnownNat n, M.Manifest r e, Num e)
           => Comp -> Matrix m k r e -> Matrix k n r e -> Matrix m n r e
matMulComp comp a b =
  case comp of
    Seq -> matMul a b
    _   -> -- For parallel: use delayed array with ikj-reordered inner product
           let kk = dimVal @k
           in makeMatrixComp @m @n @r comp $ \i j ->
             foldl' (\acc l -> acc + (a ! (i, l)) * (b ! (l, j))) 0 [0..kk-1]

-- | Matrix transpose.
transpose :: forall m n r e. (KnownNat m, KnownNat n, M.Manifest r e)
          => Matrix m n r e -> Matrix n m r e
transpose (MkMatrix arr) = MkMatrix $ M.compute $ M.transposeInner arr

-- | P-specialised raw-primop matrix transpose.
-- Avoids per-element overhead of massiv's delayed transpose.
transposeP :: forall m n. (KnownNat m, KnownNat n)
           => Matrix m n M.P Double -> Matrix n m M.P Double
transposeP (MkMatrix a) =
  let !mm = dimVal @m
      !nn = dimVal @n
      !(I# mm#) = mm
      !(I# nn#) = nn
      !(ByteArray ba#) = unwrapByteArray a
      !(I# off#) = unwrapByteArrayOffset a
  in createMatrix @n @m @M.P $ \mu ->
    let !(MutableByteArray mba#) = unwrapMutableByteArray mu
        !(I# offR#) = unwrapMutableByteArrayOffset mu
    in ST $ \s0 ->
      -- Iterate source rows (sequential read, strided write)
      let goRow j s
            | isTrue# (j >=# mm#) = s
            | otherwise =
                let goCol i s1
                      | isTrue# (i >=# nn#) = s1
                      | otherwise =
                          let !src = off# +# j *# nn# +# i
                              !dst = offR# +# i *# mm# +# j
                              !v = indexDoubleArray# ba# src
                          in case writeDoubleArray# mba# dst v s1 of
                               s2 -> goCol (i +# 1#) s2
                in goRow (j +# 1#) (goCol 0# s)
      in (# goRow 0# s0, () #)
{-# NOINLINE transposeP #-}

-- | P-specialised A^T * A without materialising A^T.
-- Computes C = A^T * A using a fused DSYRK kernel that processes only the
-- lower triangle (halving flops) and mirrors to the upper triangle.
-- Avoids materialising A^T entirely — one allocation, no transpose pass.
matMulAtAP :: forall m n. (KnownNat m, KnownNat n)
           => Matrix m n M.P Double -> Matrix n n M.P Double
matMulAtAP (MkMatrix arrA) =
  let mm = dimVal @m
      nn = dimVal @n
      baA = unwrapByteArray arrA
      offA = unwrapByteArrayOffset arrA
  in createMatrix @n @n @M.P $ \mc -> do
    -- Zero-initialise C (n×n)
    let mbaC = unwrapMutableByteArray mc
        offC = unwrapMutableByteArrayOffset mc
        !(I# nn#) = nn
    ST $ \s0 ->
      let go i s
            | isTrue# (i >=# (nn# *# nn#)) = s
            | otherwise = case writeDoubleArray# (unMBA# mbaC) (unI offC +# i) 0.0## s of
                            s' -> go (i +# 1#) s'
      in (# go 0# s0, () #)
    -- Run the fused SYRK kernel: C = A^T * A (lower triangle + mirror)
    rawSyrkLowerKernel baA offA mbaC offC mm nn
{-# NOINLINE matMulAtAP #-}

-- | Element-wise matrix addition.
mAdd :: (KnownNat m, KnownNat n, M.Manifest r e, Num e)
     => Matrix m n r e -> Matrix m n r e -> Matrix m n r e
mAdd (MkMatrix a) (MkMatrix b) = MkMatrix $ M.compute $ M.zipWith (+) a b

-- | Element-wise matrix subtraction.
mSub :: (KnownNat m, KnownNat n, M.Manifest r e, Num e)
     => Matrix m n r e -> Matrix m n r e -> Matrix m n r e
mSub (MkMatrix a) (MkMatrix b) = MkMatrix $ M.compute $ M.zipWith (-) a b

-- | Scalar–matrix multiply.
mScale :: (KnownNat m, KnownNat n, M.Manifest r e, Num e)
       => e -> Matrix m n r e -> Matrix m n r e
mScale alpha (MkMatrix a) = MkMatrix $ M.compute $ M.map (* alpha) a

-- | Left-multiply by a lower-triangular matrix (BLAS @TRMM@, left side).
trmmLeft :: forall n r e. (KnownNat n, M.Manifest r e, Num e)
         => Matrix n n r e -> Matrix n n r e -> Matrix n n r e
trmmLeft l b =
  let nn = dimVal @n
  in makeMatrix @n @n @r $ \i j ->
    foldl' (\acc k -> acc + (l ! (i, k)) * (b ! (k, j))) 0 [0..min i (nn-1)]

-- Helpers to unwrap newtypes for raw primop access
unMBA# :: MutableByteArray s -> MutableByteArray# s
unMBA# (MutableByteArray mba) = mba
{-# INLINE unMBA# #-}

unI :: Int -> Int#
unI (I# i) = i
{-# INLINE unI #-}
