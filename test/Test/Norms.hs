{-# LANGUAGE AllowAmbiguousTypes #-}

module Test.Norms (normTests) where

import Test.Tasty
import Test.Tasty.QuickCheck
import Test.Tasty.HUnit
import qualified Data.Massiv.Array as M

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.Norms
import Numeric.LinearAlgebra.Massiv.BLAS.Level3 (transpose, matMul)
import Test.Types

normTests :: TestTree
normTests = testGroup "Norms"
  [ testGroup "Vector norms"
    [ testProperty "vnorm2 non-negative" prop_vnorm2NonNeg
    , testProperty "vnorm1 non-negative" prop_vnorm1NonNeg
    , testProperty "vnormInf <= vnorm1" prop_vnormInfLeqVnorm1
    , testCase "known vector norm" test_knownVnorm
    ]
  , testGroup "Matrix norms"
    [ testProperty "normFrob non-negative" prop_normFrobNonNeg
    , testProperty "norm1 = normInf of transpose" prop_norm1TransposeNormInf
    , testProperty "triangle inequality (Frobenius)" prop_triangleInequality
    , testCase "known Frobenius norm" test_knownFrobNorm
    ]
  ]

-- Vector norms

prop_vnorm2NonNeg :: Property
prop_vnorm2NonNeg = forAll (genVector @4) $ \x ->
  vnorm2 x >= (0 :: Double)

prop_vnorm1NonNeg :: Property
prop_vnorm1NonNeg = forAll (genVector @4) $ \x ->
  vnorm1 x >= (0 :: Double)

prop_vnormInfLeqVnorm1 :: Property
prop_vnormInfLeqVnorm1 = forAll (genVector @4) $ \x ->
  vnormInf x <= vnorm1 x + 1e-12

test_knownVnorm :: Assertion
test_knownVnorm = do
  let v = makeVector @3 @M.P $ \i -> case i of { 0 -> 3; 1 -> 4; _ -> 0 } :: Double
  assertBool "vnorm2 of [3,4,0] = 5" $ abs (vnorm2 v - 5) < 1e-12
  assertBool "vnorm1 of [3,4,0] = 7" $ abs (vnorm1 v - 7) < 1e-12
  assertBool "vnormInf of [3,4,0] = 4" $ abs (vnormInf v - 4) < 1e-12

-- Matrix norms

prop_normFrobNonNeg :: Property
prop_normFrobNonNeg = forAll (genMatrix @3 @3) $ \a ->
  normFrob a >= (0 :: Double)

prop_norm1TransposeNormInf :: Property
prop_norm1TransposeNormInf = forAll (genMatrix @3 @4) $ \a ->
  let at = transpose a
  in abs (norm1 a - normInf at) < 1e-8

prop_triangleInequality :: Property
prop_triangleInequality = forAll ((,) <$> genMatrix @3 @3 <*> genMatrix @3 @3) $ \(a, b) ->
  let ab = makeMatrix @3 @3 @M.P $ \i j -> (a ! (i,j)) + (b ! (i,j))
  in normFrob ab <= normFrob a + normFrob b + 1e-10

test_knownFrobNorm :: Assertion
test_knownFrobNorm = do
  -- ‖[[1,2],[3,4]]‖_F = sqrt(1+4+9+16) = sqrt(30)
  let a = makeMatrix @2 @2 @M.P $ \i j -> case (i,j) of
            (0,0) -> 1; (0,1) -> 2; (1,0) -> 3; (1,1) -> 4; _ -> 0 :: Double
  assertBool "Frobenius norm" $ abs (normFrob a - sqrt 30) < 1e-12
