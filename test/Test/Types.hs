{-# LANGUAGE AllowAmbiguousTypes #-}

-- | Test helpers: Arbitrary instances, approximate equality, matrix generators.
module Test.Types
  ( -- * Approximate equality
    (~=)
  , matApproxEq
  , vecApproxEq
    -- * Matrix generators
  , genMatrix
  , genVector
  , genSPDMatrix
  , genUpperTriangular
  , genLowerTriangular
  , genSymmetric
    -- * Tolerance
  , defaultTol
  ) where

import qualified Data.Massiv.Array as M
import Data.Massiv.Array (Ix2(..), Sz(..), Comp(..))
import GHC.TypeNats (KnownNat, natVal)
import Data.Proxy (Proxy(..))
import Test.QuickCheck

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal

-- | Default tolerance for floating-point comparisons.
defaultTol :: Double
defaultTol = 1e-8

-- | Approximate equality for scalars.
(~=) :: Double -> Double -> Bool
x ~= y = abs (x - y) < defaultTol * (1 + abs x + abs y)

-- | Approximate equality for matrices.
matApproxEq :: forall m n. (KnownNat m, KnownNat n)
            => Matrix m n M.P Double -> Matrix m n M.P Double -> Bool
matApproxEq a b =
  let r = dimVal @m
      c = dimVal @n
  in all (\(i, j) -> (a ! (i, j)) ~= (b ! (i, j)))
     [(i, j) | i <- [0..r-1], j <- [0..c-1]]

-- | Approximate equality for vectors.
vecApproxEq :: forall n. KnownNat n
            => Vector n M.P Double -> Vector n M.P Double -> Bool
vecApproxEq a b =
  let nn = dimVal @n
  in all (\i -> (a !. i) ~= (b !. i)) [0..nn-1]

-- | Generate a random m×n matrix with entries in [-10, 10].
genMatrix :: forall m n. (KnownNat m, KnownNat n) => Gen (Matrix m n M.P Double)
genMatrix = do
  let r = fromIntegral (natVal (Proxy @m))
      c = fromIntegral (natVal (Proxy @n))
  entries <- vectorOf (r * c) (choose (-10, 10))
  pure $ makeMatrix @m @n @M.P $ \i j -> entries !! (i * c + j)

-- | Generate a random n-element vector with entries in [-10, 10].
genVector :: forall n. KnownNat n => Gen (Vector n M.P Double)
genVector = do
  let nn = fromIntegral (natVal (Proxy @n))
  entries <- vectorOf nn (choose (-10, 10))
  pure $ makeVector @n @M.P $ \i -> entries !! i

-- | Generate a symmetric positive definite matrix: A = BBᵀ + εI.
genSPDMatrix :: forall n. KnownNat n => Gen (Matrix n n M.P Double)
genSPDMatrix = do
  let nn = fromIntegral (natVal (Proxy @n))
  entries <- vectorOf (nn * nn) (choose (-5, 5))
  let b = makeMatrix @n @n @M.P $ \i j -> entries !! (i * nn + j)
      -- BBᵀ + εI
      epsilon = 0.1 :: Double
  pure $ makeMatrix @n @n @M.P $ \i j ->
    let bbT = foldl' (\acc k -> acc + (b ! (i, k)) * (b ! (j, k))) 0 [0..nn-1]
    in bbT + if i == j then epsilon else 0

-- | Generate an upper triangular matrix with nonzero diagonal.
genUpperTriangular :: forall n. KnownNat n => Gen (Matrix n n M.P Double)
genUpperTriangular = do
  let nn = fromIntegral (natVal (Proxy @n))
  entries <- vectorOf (nn * nn) (choose (-10, 10))
  diags <- vectorOf nn (choose (1, 10))  -- Ensure nonzero diagonal
  pure $ makeMatrix @n @n @M.P $ \i j ->
    if i == j then diags !! i
    else if i < j then entries !! (i * nn + j)
    else 0

-- | Generate a lower triangular matrix with nonzero diagonal.
genLowerTriangular :: forall n. KnownNat n => Gen (Matrix n n M.P Double)
genLowerTriangular = do
  let nn = fromIntegral (natVal (Proxy @n))
  entries <- vectorOf (nn * nn) (choose (-10, 10))
  diags <- vectorOf nn (choose (1, 10))
  pure $ makeMatrix @n @n @M.P $ \i j ->
    if i == j then diags !! i
    else if i > j then entries !! (i * nn + j)
    else 0

-- | Generate a symmetric matrix.
genSymmetric :: forall n. KnownNat n => Gen (Matrix n n M.P Double)
genSymmetric = do
  let nn = fromIntegral (natVal (Proxy @n))
  entries <- vectorOf (nn * nn) (choose (-10, 10))
  pure $ makeMatrix @n @n @M.P $ \i j ->
    if i <= j then entries !! (i * nn + j)
    else entries !! (j * nn + i)
