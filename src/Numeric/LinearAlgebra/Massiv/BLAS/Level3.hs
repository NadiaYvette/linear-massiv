{-# LANGUAGE AllowAmbiguousTypes #-}

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
  , matMulComp
    -- * Elementary matrix operations (GVL4 Section 1.1, pp. 4–18)
  , transpose
  , mAdd
  , mSub
  , mScale
    -- * Triangular matrix multiply (GVL4 Section 1.1.8, p. 15)
  , trmmLeft
  ) where

import qualified Data.Massiv.Array as M
import Data.Massiv.Array (Ix2(..), Sz(..), Comp(..))
import GHC.TypeNats (KnownNat)
import Control.Monad.ST (ST)

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal

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
--
-- Each entry of the result is
-- \(c_{ij} \leftarrow \alpha \sum_{l=1}^{k} a_{il} \, b_{lj} + \beta \, c_{ij}\),
-- using the ijk (row-oriented inner-product) loop ordering of
-- Algorithm 1.1.5.
--
-- Uses massiv's delayed array computation for natural parallelism
-- when called with a parallel 'Comp' strategy.
--
-- ==== Type-safety guarantees
--
-- * \(A\) is @m x k@, \(B\) is @k x n@, \(C\) and the result are @m x n@.
-- * The shared inner dimension @k@ is enforced at compile time, making
--   a dimension mismatch a type error.
--
-- ==== Complexity
--
-- \(O(m \, k \, n)\) — two floating-point operations per
-- \((i, l, j)\) triple in the inner product, plus \(O(m \, n)\)
-- work for the \(\alpha / \beta\) scaling.
-- | Block size for cache-tiled GEMM.  Chosen so that three bs×bs blocks
-- of Double (3 × bs² × 8 bytes) fit comfortably in L1 cache (typically 32–48 KiB).
gemmBlockSize :: Int
gemmBlockSize = 32
{-# INLINE gemmBlockSize #-}

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
--
-- This is a convenience wrapper equivalent to @'gemm' 1 a b 0 zero@
-- but avoids requiring an initial \(C\) matrix.
--
-- ==== Type-safety guarantees
--
-- * \(A\) is @m x k@, \(B\) is @k x n@, the result is @m x n@.
-- * The inner dimension @k@ is checked at compile time.
--
-- ==== Complexity
--
-- \(O(m \, k \, n)\).
matMul :: forall m k n r e. (KnownNat m, KnownNat k, KnownNat n, M.Manifest r e, Num e)
       => Matrix m k r e -> Matrix k n r e -> Matrix m n r e
matMul a b =
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

-- | Matrix multiply with explicit computation strategy.
--
-- __GVL4 Reference:__ Algorithm 1.1.5 (ijk matrix multiply), p. 12.
--
-- Identical to 'matMul' but accepts a 'Comp' argument that controls
-- the massiv computation strategy:
--
-- * @Seq@ — sequential evaluation (single-threaded).
-- * @Par@ — automatic parallel evaluation; massiv will partition the
--   work across available cores.
--
-- Given \(A \in \mathbb{R}^{m \times k}\) and
-- \(B \in \mathbb{R}^{k \times n}\), computes
--
-- \[
--   C \;=\; A \, B
-- \]
--
-- ==== Type-safety guarantees
--
-- * \(A\) is @m x k@, \(B\) is @k x n@, the result is @m x n@.
-- * The inner dimension @k@ is checked at compile time.
--
-- ==== Complexity
--
-- \(O(m \, k \, n)\).  With @Par@, wall-clock time scales roughly
-- as \(O(m \, k \, n / p)\) for \(p\) cores, subject to memory
-- bandwidth.
matMulComp :: forall m k n r e. (KnownNat m, KnownNat k, KnownNat n, M.Manifest r e, Num e)
           => Comp -> Matrix m k r e -> Matrix k n r e -> Matrix m n r e
matMulComp comp a b =
  case comp of
    Seq -> matMul a b
    _   -> -- For parallel: use delayed array with ikj-reordered inner product
           -- (each element still computed independently for parallelism)
           let kk = dimVal @k
           in makeMatrixComp @m @n @r comp $ \i j ->
             foldl' (\acc l -> acc + (a ! (i, l)) * (b ! (l, j))) 0 [0..kk-1]

-- | Matrix transpose.
--
-- __GVL4 Reference:__ Section 1.1, p. 5 (transpose notation and
-- properties).
--
-- Given \(A \in \mathbb{R}^{m \times n}\), produces
-- \(A^{T} \in \mathbb{R}^{n \times m}\) where
--
-- \[
--   (A^{T})_{ij} \;=\; a_{ji}
-- \]
--
-- Implemented via massiv's 'M.transposeInner', which performs the
-- index permutation lazily and materialises the result in a single
-- pass.
--
-- ==== Type-safety guarantees
--
-- The input is @m x n@ and the result is @n x m@; the dimension swap
-- is reflected at the type level.
--
-- ==== Complexity
--
-- \(O(m \, n)\) — one element copy per entry.
transpose :: forall m n r e. (KnownNat m, KnownNat n, M.Manifest r e)
          => Matrix m n r e -> Matrix n m r e
transpose (MkMatrix arr) = MkMatrix $ M.compute $ M.transposeInner arr

-- | Element-wise matrix addition.
--
-- __GVL4 Reference:__ Section 1.1, p. 6 (matrix arithmetic).
--
-- Given \(A, B \in \mathbb{R}^{m \times n}\), computes
--
-- \[
--   C \;=\; A + B, \qquad c_{ij} = a_{ij} + b_{ij}
-- \]
--
-- ==== Type-safety guarantees
--
-- Both operands and the result share the same compile-time dimensions
-- @m x n@.
--
-- ==== Complexity
--
-- \(O(m \, n)\).
mAdd :: (KnownNat m, KnownNat n, M.Manifest r e, Num e)
     => Matrix m n r e -> Matrix m n r e -> Matrix m n r e
mAdd (MkMatrix a) (MkMatrix b) = MkMatrix $ M.compute $ M.zipWith (+) a b

-- | Element-wise matrix subtraction.
--
-- __GVL4 Reference:__ Section 1.1, p. 6 (matrix arithmetic).
--
-- Given \(A, B \in \mathbb{R}^{m \times n}\), computes
--
-- \[
--   C \;=\; A - B, \qquad c_{ij} = a_{ij} - b_{ij}
-- \]
--
-- ==== Type-safety guarantees
--
-- Both operands and the result share the same compile-time dimensions
-- @m x n@.
--
-- ==== Complexity
--
-- \(O(m \, n)\).
mSub :: (KnownNat m, KnownNat n, M.Manifest r e, Num e)
     => Matrix m n r e -> Matrix m n r e -> Matrix m n r e
mSub (MkMatrix a) (MkMatrix b) = MkMatrix $ M.compute $ M.zipWith (-) a b

-- | Scalar–matrix multiply.
--
-- __GVL4 Reference:__ Section 1.1, p. 6 (scalar–matrix operations).
--
-- Given a scalar \(\alpha\) and \(A \in \mathbb{R}^{m \times n}\),
-- computes
--
-- \[
--   B \;=\; \alpha \, A, \qquad b_{ij} = \alpha \, a_{ij}
-- \]
--
-- ==== Type-safety guarantees
--
-- The result retains the same compile-time dimensions @m x n@ as the
-- input.
--
-- ==== Complexity
--
-- \(O(m \, n)\) — one multiplication per entry.
mScale :: (KnownNat m, KnownNat n, M.Manifest r e, Num e)
       => e -> Matrix m n r e -> Matrix m n r e
mScale alpha (MkMatrix a) = MkMatrix $ M.compute $ M.map (* alpha) a

-- | Left-multiply by a lower-triangular matrix (BLAS @TRMM@, left side).
--
-- __GVL4 Reference:__ Section 1.1.8, p. 15 (triangular matrix
-- operations).  See also Section 3.1.1 (Forward Substitution,
-- p. 106) for context on triangular structure exploitation.
--
-- Given a lower-triangular matrix \(L \in \mathbb{R}^{n \times n}\)
-- and a general matrix \(B \in \mathbb{R}^{n \times n}\), computes
--
-- \[
--   B \;\leftarrow\; L \, B
-- \]
--
-- Only the lower-triangular entries of \(L\) (i.e.,
-- \(l_{ik}\) with \(k \leq i\)) are accessed, so any values stored
-- in the strict upper triangle are ignored.  For each entry of the
-- result:
--
-- \[
--   c_{ij} \;=\; \sum_{k=0}^{i} l_{ik} \, b_{kj}
-- \]
--
-- ==== Type-safety guarantees
--
-- Both \(L\) and \(B\) are @n x n@ (square), and the result is
-- @n x n@.  The square constraint is enforced at compile time.
--
-- ==== Complexity
--
-- \(O(n^{3} / 2)\) — roughly half the work of a full @GEMM@ because
-- the upper-triangular entries of \(L\) are zero and skipped.
trmmLeft :: forall n r e. (KnownNat n, M.Manifest r e, Num e)
         => Matrix n n r e -> Matrix n n r e -> Matrix n n r e
trmmLeft l b =
  let nn = dimVal @n
  in makeMatrix @n @n @r $ \i j ->
    foldl' (\acc k -> acc + (l ! (i, k)) * (b ! (k, j))) 0 [0..min i (nn-1)]
