{-# LANGUAGE AllowAmbiguousTypes #-}

-- |
-- Module      : Numeric.LinearAlgebra.Massiv.Linear
-- Copyright   : (c) Nadia Chambers 2026
-- License     : BSD-3-Clause
-- Maintainer  : nadia.chambers@iohk.io
-- Stability   : experimental
--
-- = Integration with the @linear@ library
--
-- This module provides conversion functions between the
-- <https://hackage.haskell.org/package/linear linear> library's types
-- (@'Linear.V.V'@, @'Linear.V2.V2'@, @'Linear.V3.V3'@, @'Linear.V4.V4'@)
-- and our dimensioned @'Vector'@ \/ @'Matrix'@ types.
--
-- == Why not typeclass instances?
--
-- The @linear@ library's typeclasses ('Linear.Additive.Additive',
-- 'Linear.Metric.Metric', 'Linear.Trace.Trace') expect types of kind
-- @* -> *@ (i.e., functors over the element type). Our @Vector n r e@ and
-- @Matrix m n r e@ carry additional type parameters (@n@, @r@) before @e@,
-- making direct functor-based instances impractical without additional
-- newtype wrappers.
--
-- Instead, equivalent operations are provided as standalone functions:
--
-- * "Numeric.LinearAlgebra.Massiv.BLAS.Level1" — 'dot', 'axpy', 'scal', 'nrm2'
-- * "Numeric.LinearAlgebra.Massiv.BLAS.Level3" — 'mAdd', 'mSub', 'mScale', 'matMul', 'transpose'
-- * "Numeric.LinearAlgebra.Massiv.Norms" — 'normFrob', 'norm1', 'normInf'
--
-- == Conversion semantics
--
-- 'fromLinearV' converts to the @'Data.Massiv.Array.B'@ (boxed) representation
-- because @linear@'s @V@ stores elements in a boxed @Data.Vector.Vector@.
-- For small fixed-size types (@V2@, @V3@, @V4@), 'fromV2' etc. produce
-- vectors in any representation @r@ satisfying @Manifest r e@.
module Numeric.LinearAlgebra.Massiv.Linear
  ( -- * Conversions with @linear@'s @V n a@
    fromLinearV
  , toLinearV
    -- * Small fixed-size vector conversions
  , fromV2
  , fromV3
  , fromV4
    -- * List-based matrix I\/O
  , toListMatrix
  , fromListMatrix
  ) where

import qualified Data.Massiv.Array as M
import Data.Massiv.Array (Ix1, Ix2(..), Sz(..))
import qualified Data.Vector as BV
import GHC.TypeNats (KnownNat, natVal)
import Data.Proxy (Proxy(..))

import qualified Linear.V as L
import Linear.V2 (V2(..))
import Linear.V3 (V3(..))
import Linear.V4 (V4(..))

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal

-- | Convert a @linear@ @'Linear.V.V' n a@ to our @'Vector' n 'Data.Massiv.Array.B' a@.
--
-- The result uses the boxed (@B@) representation since @linear@'s @V@ is
-- backed by a boxed @Data.Vector.Vector@.
fromLinearV :: forall n a. KnownNat n => L.V n a -> Vector n M.B a
fromLinearV lv =
  let bv = L.toVector lv
      nn = fromIntegral (natVal (Proxy @n))
      arr = M.compute @M.B $ M.makeArray @M.D M.Seq (M.Sz1 nn) (bv BV.!)
  in MkVector arr

-- | Convert our @'Vector' n 'Data.Massiv.Array.B' a@ to a @linear@ @'Linear.V.V' n a@.
--
-- This is the inverse of 'fromLinearV'. The dimension is checked by @linear@'s
-- 'Linear.V.fromVector' (which returns 'Maybe'); the 'error' case is unreachable
-- given correct type-level dimensions.
toLinearV :: forall n a. KnownNat n => Vector n M.B a -> L.V n a
toLinearV (MkVector arr) =
  let nn = fromIntegral (natVal (Proxy @n))
      bv = BV.generate nn (\i -> M.index' arr i)
  in case L.fromVector bv of
    Just v  -> v
    Nothing -> error "toLinearV: impossible dimension mismatch"

-- | Convert a @linear@ @'Linear.V2.V2'@ to a 2-element 'Vector'.
fromV2 :: M.Manifest r e => V2 e -> Vector 2 r e
fromV2 (V2 x y) = makeVector @2 $ \i -> case i of { 0 -> x; _ -> y }

-- | Convert a @linear@ @'Linear.V3.V3'@ to a 3-element 'Vector'.
fromV3 :: M.Manifest r e => V3 e -> Vector 3 r e
fromV3 (V3 x y z) = makeVector @3 $ \i -> case i of { 0 -> x; 1 -> y; _ -> z }

-- | Convert a @linear@ @'Linear.V4.V4'@ to a 4-element 'Vector'.
fromV4 :: M.Manifest r e => V4 e -> Vector 4 r e
fromV4 (V4 x y z w) = makeVector @4 $ \i -> case i of { 0 -> x; 1 -> y; 2 -> z; _ -> w }

-- | Convert a matrix to a list of lists (row-major order).
--
-- @
-- toListMatrix mat  ==  [[mat '!' (i,j) | j <- [0..n-1]] | i <- [0..m-1]]
-- @
toListMatrix :: forall m n r e. (KnownNat m, KnownNat n, M.Manifest r e)
             => Matrix m n r e -> [[e]]
toListMatrix mat =
  let r = dimVal @m
      c = dimVal @n
  in [[mat ! (i, j) | j <- [0..c-1]] | i <- [0..r-1]]

-- | Create a matrix from a list of lists (row-major order).
--
-- Returns 'Nothing' if the list dimensions do not match the type-level
-- dimensions \(m\) and \(n\).
fromListMatrix :: forall m n r e. (KnownNat m, KnownNat n, M.Manifest r e)
               => [[e]] -> Maybe (Matrix m n r e)
fromListMatrix rows_
  | length rows_ /= dimVal @m = Nothing
  | any (\row -> length row /= dimVal @n) rows_ = Nothing
  | otherwise = Just $ makeMatrix @m @n @r $ \i j -> (rows_ !! i) !! j
