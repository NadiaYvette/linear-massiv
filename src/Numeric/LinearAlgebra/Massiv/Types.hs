{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE AllowAmbiguousTypes #-}

-- |
-- Module      : Numeric.LinearAlgebra.Massiv.Types
-- Copyright   : (c) Nadia Chambers 2026
-- License     : BSD-3-Clause
-- Maintainer  : nadia.chambers@iohk.io
-- Stability   : experimental
--
-- Core type definitions for type-safe dimensioned matrices and vectors
-- backed by <https://hackage.haskell.org/package/massiv massiv> arrays.
--
-- The central types are 'Matrix' and 'Vector', which wrap massiv arrays
-- with phantom type-level natural number parameters encoding their dimensions.
-- This ensures that dimensionally incorrect operations (e.g., adding matrices
-- of different sizes, or multiplying matrices with incompatible inner dimensions)
-- are caught at compile time by GHC's type checker.
--
-- = Type-level dimension encoding
--
-- Dimensions are encoded as GHC @DataKinds@ promoted @'GHC.TypeNats.Nat'@
-- values. The constraint @'GHC.TypeNats.KnownNat' n@ provides access to the
-- runtime value via @'GHC.TypeNats.natVal'@.
--
-- @
-- -- A 3x4 matrix of Doubles using Primitive representation
-- type MyMatrix = Matrix 3 4 P Double
--
-- -- A 5-element vector of Doubles using Unboxed representation
-- type MyVector = Vector 5 U Double
-- @
--
-- = Existential wrappers
--
-- For situations where dimensions are not known until runtime (e.g., reading
-- a matrix from a file), use 'SomeMatrix' and 'SomeVector'. These existentially
-- quantify the dimension parameters while retaining 'KnownNat' evidence.
--
-- See "Numeric.LinearAlgebra.Massiv.Internal" for construction helpers.
module Numeric.LinearAlgebra.Massiv.Types
  ( -- * Dimensioned matrix type
    Matrix(..)
    -- * Dimensioned vector type
  , Vector(..)
    -- * Smart constructors
  , matrix
  , vector
    -- * Existential wrappers
  , SomeMatrix(..)
  , SomeVector(..)
  , someMatrix
  , someVector
    -- * Dimension queries
  , rows
  , cols
  , size
    -- * Type-level helpers
  , type KnownDims
  ) where

import Data.Massiv.Array (Array, Ix2(..), Sz(..), Ix1, Comp(..))
import qualified Data.Massiv.Array as M
import GHC.TypeNats (Nat, KnownNat, natVal, SomeNat(..), someNatVal)
import Data.Proxy (Proxy(..))
import Control.DeepSeq (NFData(..))

-- | Constraint synonym for two known dimensions.
--
-- @KnownDims m n@ is equivalent to @(KnownNat m, KnownNat n)@.
type KnownDims m n = (KnownNat m, KnownNat n)

-- | A matrix with compile-time known dimensions \(m\) (rows) \(\times\) \(n\) (cols).
--
-- Wraps a massiv @'Data.Massiv.Array.Array' r 'Data.Massiv.Array.Ix2' e@.
-- The phantom type parameters \(m\) and \(n\) enforce dimensional conformance
-- at compile time. For example, matrix multiplication via 'matMul' requires
-- the inner dimensions to unify:
--
-- @
-- 'matMul' :: Matrix m __k__ r e -> Matrix __k__ n r e -> Matrix m n r e
-- @
--
-- The representation parameter @r@ selects the massiv array backend:
--
-- * @'Data.Massiv.Array.P'@ — Primitive (best for 'Double', 'Int'; pinned memory)
-- * @'Data.Massiv.Array.U'@ — Unboxed (via @Data.Vector.Unboxed@)
-- * @'Data.Massiv.Array.S'@ — Storable (via @Foreign.ForeignPtr@; useful for FFI)
-- * @'Data.Massiv.Array.B'@ — Boxed (polymorphic but slower; GC overhead)
newtype Matrix (m :: Nat) (n :: Nat) r e = MkMatrix { unMatrix :: Array r Ix2 e }

deriving instance Show (Array r Ix2 e) => Show (Matrix m n r e)
deriving instance Eq (Array r Ix2 e) => Eq (Matrix m n r e)

instance NFData (Array r Ix2 e) => NFData (Matrix m n r e) where
  rnf (MkMatrix arr) = rnf arr

-- | A vector with compile-time known dimension \(n\).
--
-- Wraps a massiv @'Data.Massiv.Array.Array' r 'Data.Massiv.Array.Ix1' e@.
-- The phantom parameter \(n\) ensures that vector operations (e.g., 'dot',
-- 'axpy') are only applied to vectors of matching dimension.
newtype Vector (n :: Nat) r e = MkVector { unVector :: Array r Ix1 e }

deriving instance Show (Array r Ix1 e) => Show (Vector n r e)
deriving instance Eq (Array r Ix1 e) => Eq (Vector n r e)

instance NFData (Array r Ix1 e) => NFData (Vector n r e) where
  rnf (MkVector arr) = rnf arr

-- | Smart constructor for matrices. Checks at runtime that the array
-- dimensions match the type-level dimensions \(m\) and \(n\).
--
-- Returns 'Nothing' if the dimensions do not match.
--
-- @
-- let arr = M.makeArray Seq (Sz2 3 4) (\\(i :. j) -> fromIntegral (i + j))
-- matrix \@3 \@4 arr  -- Just (MkMatrix arr)
-- matrix \@2 \@4 arr  -- Nothing
-- @
matrix :: forall m n r e. (KnownDims m n, M.Size r)
       => Array r Ix2 e -> Maybe (Matrix m n r e)
matrix arr
  | M.Sz2 r c <- M.size arr
  , r == fromIntegral (natVal (Proxy @m))
  , c == fromIntegral (natVal (Proxy @n))
  = Just (MkMatrix arr)
  | otherwise = Nothing

-- | Smart constructor for vectors. Checks at runtime that the array
-- size matches the type-level dimension \(n\).
--
-- Returns 'Nothing' if the size does not match.
vector :: forall n r e. (KnownNat n, M.Size r)
       => Array r Ix1 e -> Maybe (Vector n r e)
vector arr
  | M.Sz1 n <- M.size arr
  , n == fromIntegral (natVal (Proxy @n))
  = Just (MkVector arr)
  | otherwise = Nothing

-- | Get the number of rows at the value level. \(O(1)\).
rows :: forall m n r e. KnownNat m => Matrix m n r e -> Int
rows _ = fromIntegral (natVal (Proxy @m))

-- | Get the number of columns at the value level. \(O(1)\).
cols :: forall m n r e. KnownNat n => Matrix m n r e -> Int
cols _ = fromIntegral (natVal (Proxy @n))

-- | Get the size of a vector at the value level. \(O(1)\).
size :: forall n r e. KnownNat n => Vector n r e -> Int
size _ = fromIntegral (natVal (Proxy @n))

-- | Existentially quantified matrix with runtime-determined dimensions.
--
-- Use 'someMatrix' to wrap a massiv array whose dimensions are not known
-- at compile time. Pattern matching on 'SomeMatrix' brings 'KnownNat'
-- evidence into scope:
--
-- @
-- case someMatrix arr of
--   SomeMatrix (mat :: Matrix m n r e) -> ...
--   -- m and n are now in scope as KnownNat
-- @
data SomeMatrix r e where
  SomeMatrix :: (KnownNat m, KnownNat n) => Matrix m n r e -> SomeMatrix r e

-- | Existentially quantified vector with runtime-determined dimensions.
data SomeVector r e where
  SomeVector :: KnownNat n => Vector n r e -> SomeVector r e

-- | Wrap a massiv 2D array into an existentially typed matrix.
someMatrix :: M.Size r => Array r Ix2 e -> SomeMatrix r e
someMatrix arr =
  let M.Sz2 r c = M.size arr
  in case someNatVal (fromIntegral r) of
    SomeNat (_ :: Proxy m) -> case someNatVal (fromIntegral c) of
      SomeNat (_ :: Proxy n) -> SomeMatrix @m @n (MkMatrix arr)

-- | Wrap a massiv 1D array into an existentially typed vector.
someVector :: M.Size r => Array r Ix1 e -> SomeVector r e
someVector arr =
  let M.Sz1 n = M.size arr
  in case someNatVal (fromIntegral n) of
    SomeNat (_ :: Proxy n) -> SomeVector @n (MkVector arr)
