{-# LANGUAGE AllowAmbiguousTypes #-}

-- |
-- Module      : Numeric.LinearAlgebra.Massiv.Internal
-- Copyright   : (c) Nadia Chambers 2026
-- License     : BSD-3-Clause
-- Maintainer  : nadia.chambers@iohk.io
-- Stability   : experimental
--
-- Internal module providing unsafe constructors, dimension reification helpers,
-- and array creation utilities. This module is re-exported by
-- "Numeric.LinearAlgebra.Massiv" for convenience, but the unsafe constructors
-- ('unsafeMatrix', 'unsafeVector') bypass the dimension checks provided by the
-- smart constructors in "Numeric.LinearAlgebra.Massiv.Types".
--
-- = Array creation patterns
--
-- Two primary patterns are provided for constructing dimensioned arrays:
--
-- 1. __Pure indexed construction__ via 'makeMatrix' and 'makeVector': supply a
--    pure function @(Int -> Int -> e)@ or @(Int -> e)@ that computes each element
--    from its indices. These use massiv's delayed (@D@) intermediate representation
--    and then @'Data.Massiv.Array.compute'@ to materialise the result.
--
-- 2. __Mutable ST construction__ via 'createMatrix' and 'createVector': supply an
--    @ST@ action operating on a mutable array. This is essential for algorithms
--    that require in-place updates (e.g., LU factorization, Cholesky).
--
-- Both patterns produce arrays with the sequential ('Data.Massiv.Array.Seq')
-- computation strategy by default. Use the @Comp@ variants ('makeMatrixComp',
-- 'makeVectorComp') for parallel construction.
module Numeric.LinearAlgebra.Massiv.Internal
  ( -- * Unsafe constructors
    unsafeMatrix
  , unsafeVector
    -- * Dimension reification
  , dimVal
  , dimVal2
  , reifyDim
  , reifyDim2
    -- * Array creation helpers (pure, sequential)
  , makeMatrix
  , makeVector
    -- * Array creation helpers (pure, with Comp)
  , makeMatrixComp
  , makeVectorComp
    -- * Array creation helpers (mutable ST)
  , createMatrix
  , createVector
  , createMatrixComp
  , createVectorComp
    -- * Mutable modification helpers
  , withMutableMatrix
  , withMutableVector
  , withMutableMatrix_
  , withMutableVector_
    -- * Indexing helpers
  , (!)
  , (!.)
    -- * Identity and zero
  , identityMatrix
  , zeroMatrix
  , zeroVector
  ) where

import Data.Massiv.Array (Array, Ix2(..), Sz(..), Ix1, Comp(..), D)
import qualified Data.Massiv.Array as M
import GHC.TypeNats (Nat, KnownNat, natVal, SomeNat(..), someNatVal)
import Data.Proxy (Proxy(..))
import Control.Monad.ST (ST)

import Numeric.LinearAlgebra.Massiv.Types

-- | Unsafe matrix constructor — wraps a massiv array with /no/ dimension check.
--
-- __Precondition__: the array must have exactly \(m\) rows and \(n\) columns.
-- Violating this precondition leads to index-out-of-bounds errors at runtime.
--
-- Prefer the safe 'matrix' constructor from "Numeric.LinearAlgebra.Massiv.Types"
-- unless you can guarantee correctness (e.g., the array was just constructed
-- with the correct dimensions).
unsafeMatrix :: Array r Ix2 e -> Matrix m n r e
unsafeMatrix = MkMatrix

-- | Unsafe vector constructor — wraps a massiv array with /no/ size check.
--
-- __Precondition__: the array must have exactly \(n\) elements.
unsafeVector :: Array r Ix1 e -> Vector n r e
unsafeVector = MkVector

-- | Get the runtime value of a type-level dimension.
--
-- @
-- dimVal \@3  ==  3
-- dimVal \@100  ==  100
-- @
dimVal :: forall n. KnownNat n => Int
dimVal = fromIntegral (natVal (Proxy @n))

-- | Get both dimensions of a matrix type as a tuple.
dimVal2 :: forall m n. (KnownNat m, KnownNat n) => (Int, Int)
dimVal2 = (dimVal @m, dimVal @n)

-- | Index into a matrix (0-based, unchecked).
--
-- @mat '!' (i, j)@ returns the element at row \(i\), column \(j\).
--
-- __Warning__: No bounds checking is performed. Out-of-bounds access
-- results in undefined behaviour for unboxed\/primitive representations.
(!) :: M.Manifest r e => Matrix m n r e -> (Int, Int) -> e
(!) (MkMatrix arr) (i, j) = M.index' arr (i :. j)

-- | Index into a vector (0-based, unchecked).
--
-- @vec '!.' i@ returns the element at position \(i\).
(!.) :: M.Manifest r e => Vector n r e -> Int -> e
(!.) (MkVector arr) i = M.index' arr i

-- | Reify a runtime 'Int' as a type-level 'GHC.TypeNats.Nat'.
--
-- The continuation receives a 'Proxy' carrying the reified type.
-- This is useful for working with matrices of runtime-determined size.
reifyDim :: Int -> (forall n. KnownNat n => Proxy n -> a) -> a
reifyDim n f = case someNatVal (fromIntegral n) of
  SomeNat p -> f p

-- | Reify two runtime 'Int's as type-level 'GHC.TypeNats.Nat's.
reifyDim2 :: Int -> Int -> (forall m n. (KnownNat m, KnownNat n) => Proxy m -> Proxy n -> a) -> a
reifyDim2 m n f = reifyDim m $ \pm -> reifyDim n $ \pn -> f pm pn

-- | Create a matrix using a pure indexing function (sequential computation).
--
-- @
-- makeMatrix \@3 \@3 \@P $ \\i j -> fromIntegral (i * 3 + j)
-- @
--
-- Internally uses massiv's @'Data.Massiv.Array.Delayed'@ representation as
-- an intermediate before computing to the target representation @r@.
makeMatrix :: forall m n r e. (KnownNat m, KnownNat n, M.Manifest r e)
           => (Int -> Int -> e) -> Matrix m n r e
makeMatrix f =
  let r = dimVal @m
      c = dimVal @n
  in MkMatrix $ M.compute @r $ M.makeArray @D Seq (M.Sz2 r c) (\(i :. j) -> f i j)

-- | Create a vector using a pure indexing function (sequential computation).
makeVector :: forall n r e. (KnownNat n, M.Manifest r e)
           => (Int -> e) -> Vector n r e
makeVector f =
  let sz = dimVal @n
  in MkVector $ M.compute @r $ M.makeArray @D Seq (M.Sz1 sz) f

-- | Create a matrix using a pure indexing function with specified
-- 'Data.Massiv.Array.Comp' strategy.
--
-- Use @Par@ for parallel construction of large matrices:
--
-- @
-- makeMatrixComp \@1000 \@1000 \@P Par $ \\i j -> ...
-- @
makeMatrixComp :: forall m n r e. (KnownNat m, KnownNat n, M.Manifest r e)
               => Comp -> (Int -> Int -> e) -> Matrix m n r e
makeMatrixComp comp f =
  let r = dimVal @m
      c = dimVal @n
  in MkMatrix $ M.compute @r $ M.makeArray @D comp (M.Sz2 r c) (\(i :. j) -> f i j)

-- | Create a vector using a pure indexing function with specified 'Comp'.
makeVectorComp :: forall n r e. (KnownNat n, M.Manifest r e)
               => Comp -> (Int -> e) -> Vector n r e
makeVectorComp comp f =
  let sz = dimVal @n
  in MkVector $ M.compute @r $ M.makeArray @D comp (M.Sz1 sz) f

-- | Create a matrix using a mutable 'ST' computation.
--
-- The action receives a pre-allocated mutable array of the correct size.
-- All writes must be within bounds. The mutable array is frozen after
-- the action completes.
--
-- This is the primary mechanism for implementing algorithms with in-place
-- updates (e.g., LU factorization, Cholesky decomposition).
createMatrix :: forall m n r e. (KnownNat m, KnownNat n, M.Manifest r e)
             => (forall s. M.MArray s r Ix2 e -> ST s ()) -> Matrix m n r e
createMatrix action =
  let r = dimVal @m
      c = dimVal @n
      arr = M.createArrayST_ (M.Sz2 r c) action
  in MkMatrix arr

-- | Create a vector using a mutable 'ST' computation.
createVector :: forall n r e. (KnownNat n, M.Manifest r e)
             => (forall s. M.MArray s r Ix1 e -> ST s ()) -> Vector n r e
createVector action =
  let sz = dimVal @n
      arr = M.createArrayST_ (M.Sz1 sz) action
  in MkVector arr

-- | Create a matrix using a mutable computation with specified 'Comp'.
createMatrixComp :: forall m n r e. (KnownNat m, KnownNat n, M.Manifest r e)
                 => Comp -> (forall s. M.MArray s r Ix2 e -> ST s ()) -> Matrix m n r e
createMatrixComp _comp action =
  -- Note: createArrayST_ always runs sequentially; Comp is for delayed computations
  createMatrix @m @n action

-- | Create a vector using a mutable computation with specified 'Comp'.
createVectorComp :: forall n r e. (KnownNat n, M.Manifest r e)
                 => Comp -> (forall s. M.MArray s r Ix1 e -> ST s ()) -> Vector n r e
createVectorComp _comp action = createVector @n action

-- | Run a mutable operation on a /copy/ of the matrix, returning both the
-- action's result and the modified matrix. The original matrix is not modified.
withMutableMatrix :: (M.Manifest r e)
                  => Matrix m n r e
                  -> (forall s. M.MArray s r Ix2 e -> ST s a)
                  -> (a, Matrix m n r e)
withMutableMatrix (MkMatrix arr) action =
  let (result, arr') = M.withMArrayST arr action
  in (result, MkMatrix arr')

-- | Run a mutable operation on a /copy/ of the vector.
withMutableVector :: (M.Manifest r e)
                  => Vector n r e
                  -> (forall s. M.MArray s r Ix1 e -> ST s a)
                  -> (a, Vector n r e)
withMutableVector (MkVector arr) action =
  let (result, arr') = M.withMArrayST arr action
  in (result, MkVector arr')

-- | Like 'withMutableMatrix' but discards the action's result.
withMutableMatrix_ :: (M.Manifest r e)
                   => Matrix m n r e
                   -> (forall s. M.MArray s r Ix2 e -> ST s ())
                   -> Matrix m n r e
withMutableMatrix_ mat action = snd $ withMutableMatrix mat action

-- | Like 'withMutableVector' but discards the action's result.
withMutableVector_ :: (M.Manifest r e)
                   => Vector n r e
                   -> (forall s. M.MArray s r Ix1 e -> ST s ())
                   -> Vector n r e
withMutableVector_ vec action = snd $ withMutableVector vec action

-- | The \(n \times n\) identity matrix \(I_n\).
--
-- \[
--   I_{ij} = \begin{cases} 1 & \text{if } i = j \\ 0 & \text{otherwise} \end{cases}
-- \]
identityMatrix :: forall n r e. (KnownNat n, M.Manifest r e, Num e)
               => Matrix n n r e
identityMatrix = makeMatrix @n @n @r $ \i j -> if i == j then 1 else 0

-- | The \(m \times n\) zero matrix.
zeroMatrix :: forall m n r e. (KnownNat m, KnownNat n, M.Manifest r e, Num e)
           => Matrix m n r e
zeroMatrix = makeMatrix @m @n @r $ \_ _ -> 0

-- | The zero vector of dimension \(n\).
zeroVector :: forall n r e. (KnownNat n, M.Manifest r e, Num e)
           => Vector n r e
zeroVector = makeVector @n @r $ const 0
