{-# LANGUAGE AllowAmbiguousTypes #-}

module Test.BLAS (blasTests) where

import Test.Tasty
import Test.Tasty.QuickCheck
import Test.Tasty.HUnit
import qualified Data.Massiv.Array as M

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.BLAS.Level1
import Numeric.LinearAlgebra.Massiv.BLAS.Level2
import Numeric.LinearAlgebra.Massiv.BLAS.Level3
import Test.Types

blasTests :: TestTree
blasTests = testGroup "BLAS"
  [ testGroup "Level 1"
    [ testProperty "dot product commutative" prop_dotCommutative
    , testProperty "dot with zero vector" prop_dotZero
    , testProperty "axpy identity" prop_axpyIdentity
    , testProperty "scal by 1" prop_scalIdentity
    , testProperty "nrm2 non-negative" prop_nrm2NonNeg
    ]
  , testGroup "Level 2"
    [ testProperty "matvec with identity" prop_matvecIdentity
    , testProperty "gemv alpha=1 beta=0" prop_gemvSimple
    ]
  , testGroup "Level 3"
    [ testProperty "matMul with identity (left)" prop_matMulIdentityLeft
    , testProperty "matMul with identity (right)" prop_matMulIdentityRight
    , testProperty "transpose involution" prop_transposeInvolution
    , testProperty "mAdd commutative" prop_mAddCommutative
    , testCase "3x3 matmul known" test_matMulKnown
    ]
  ]

-- Level 1 properties

prop_dotCommutative :: Property
prop_dotCommutative = forAll ((,) <$> genVector @4 <*> genVector @4) $ \(x, y) ->
  dot x y ~= dot y x

prop_dotZero :: Property
prop_dotZero = forAll (genVector @4) $ \x ->
  let z = zeroVector @4 @M.P :: Vector 4 M.P Double
  in dot x z ~= 0

prop_axpyIdentity :: Property
prop_axpyIdentity = forAll (genVector @4) $ \x ->
  let z = zeroVector @4 @M.P :: Vector 4 M.P Double
  in vecApproxEq @4 (axpy 1 z x) x

prop_scalIdentity :: Property
prop_scalIdentity = forAll (genVector @4) $ \x ->
  vecApproxEq @4 (scal 1 x) x

prop_nrm2NonNeg :: Property
prop_nrm2NonNeg = forAll (genVector @4) $ \x ->
  nrm2 x >= (0 :: Double)

-- Level 2 properties

prop_matvecIdentity :: Property
prop_matvecIdentity = forAll (genVector @3) $ \x ->
  let eye = identityMatrix @3 @M.P :: Matrix 3 3 M.P Double
  in vecApproxEq @3 (matvec eye x) x

prop_gemvSimple :: Property
prop_gemvSimple = forAll ((,) <$> genMatrix @3 @3 <*> genVector @3) $ \(a, x) ->
  let z = zeroVector @3 @M.P :: Vector 3 M.P Double
      result = gemv 1.0 a x 0.0 z
      expected = matvec a x
  in vecApproxEq @3 result expected

-- Level 3 properties

prop_matMulIdentityLeft :: Property
prop_matMulIdentityLeft = forAll (genMatrix @3 @3) $ \a ->
  let eye = identityMatrix @3 @M.P :: Matrix 3 3 M.P Double
  in matApproxEq @3 @3 (matMul eye a) a

prop_matMulIdentityRight :: Property
prop_matMulIdentityRight = forAll (genMatrix @3 @3) $ \a ->
  let eye = identityMatrix @3 @M.P :: Matrix 3 3 M.P Double
  in matApproxEq @3 @3 (matMul a eye) a

prop_transposeInvolution :: Property
prop_transposeInvolution = forAll (genMatrix @3 @4) $ \a ->
  matApproxEq @3 @4 (transpose (transpose a)) a

prop_mAddCommutative :: Property
prop_mAddCommutative = forAll ((,) <$> genMatrix @3 @3 <*> genMatrix @3 @3) $ \(a, b) ->
  matApproxEq @3 @3 (mAdd a b) (mAdd b a)

-- Known-value test for 3×3 matrix multiplication
test_matMulKnown :: Assertion
test_matMulKnown = do
  -- A = [[1,2],[3,4]], B = [[5,6],[7,8]]
  -- AB = [[19,22],[43,50]]
  let a = makeMatrix @2 @2 @M.P $ \i j -> case (i,j) of
            (0,0) -> 1; (0,1) -> 2; (1,0) -> 3; (1,1) -> 4; _ -> 0 :: Double
      b = makeMatrix @2 @2 @M.P $ \i j -> case (i,j) of
            (0,0) -> 5; (0,1) -> 6; (1,0) -> 7; (1,1) -> 8; _ -> 0 :: Double
      c = matMul a b
  assertBool "c(0,0) = 19" $ (c ! (0,0)) ~= 19
  assertBool "c(0,1) = 22" $ (c ! (0,1)) ~= 22
  assertBool "c(1,0) = 43" $ (c ! (1,0)) ~= 43
  assertBool "c(1,1) = 50" $ (c ! (1,1)) ~= 50
