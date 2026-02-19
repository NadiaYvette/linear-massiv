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
import Numeric.LinearAlgebra.Massiv.Norms (normFrob)
import Test.Types (genMatrix, genVector, hilbertMatrix, (~=), matApproxEq, vecApproxEq)
import Test.Residuals (machineEps)

blasTests :: TestTree
blasTests = testGroup "BLAS"
  [ testGroup "Level 1"
    [ testProperty "dot product commutative" prop_dotCommutative
    , testProperty "dot with zero vector" prop_dotZero
    , testCase "dot Hilbert columns" test_dotHilbertColumns
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
    , testProperty "matMul associative 5x5" prop_matMulAssociative5
    , testProperty "matMul with identity 10x10" prop_matMulIdentity10
    , testCase "3x3 matmul known" test_matMulKnown
    , testCase "gemm larger 3x3" test_gemmLarger
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

-- | (AB)C ≈ A(BC) for random 5×5 matrices.
prop_matMulAssociative5 :: Property
prop_matMulAssociative5 =
  forAll ((,,) <$> genMatrix @5 @5 <*> genMatrix @5 @5 <*> genMatrix @5 @5) $ \(a, b, c) ->
    let lhs = matMul (matMul a b) c
        rhs = matMul a (matMul b c)
    in normFrob (mSub lhs rhs) / (normFrob lhs + 1e-15) < 1e-6

-- | I·A = A for 10×10 matrices.
prop_matMulIdentity10 :: Property
prop_matMulIdentity10 = forAll (genMatrix @10 @10) $ \a ->
  let eye = identityMatrix @10 @M.P :: Matrix 10 10 M.P Double
  in matApproxEq @10 @10 (matMul eye a) a

-- | Dot product of columns 0 and 1 of hilbertMatrix @5.
-- col0 = [1, 1/2, 1/3, 1/4, 1/5]
-- col1 = [1/2, 1/3, 1/4, 1/5, 1/6]
-- dot  = sum_{k=0}^{4} 1/((k+1)*(k+2)) = 1/2 + 1/6 + 1/12 + 1/20 + 1/30 = 50/60 = 5/6
test_dotHilbertColumns :: Assertion
test_dotHilbertColumns = do
  let h = hilbertMatrix @5 :: Matrix 5 5 M.P Double
      col0 = makeVector @5 @M.P $ \k -> h ! (k, 0)
      col1 = makeVector @5 @M.P $ \k -> h ! (k, 1)
      result = dot col0 col1
      expected = 5 / 6 :: Double
  assertBool ("dot of Hilbert cols 0,1 = 5/6, got " ++ show result)
    $ abs (result - expected) < 1e-12

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

-- | GEMM with α=2.0, β=0.5 on known 3×3 matrices.
-- A = [[1,2,3],[4,5,6],[7,8,9]]
-- B = [[9,8,7],[6,5,4],[3,2,1]]
-- C = [[1,0,0],[0,1,0],[0,0,1]]
-- Result = α*A*B + β*C
--
-- A*B = [[30,24,18],[84,69,54],[138,114,90]]
-- α*A*B = [[60,48,36],[168,138,108],[276,228,180]]
-- β*C   = [[0.5,0,0],[0,0.5,0],[0,0,0.5]]
-- Final = [[60.5,48,36],[168,138.5,108],[276,228,180.5]]
test_gemmLarger :: Assertion
test_gemmLarger = do
  let a = makeMatrix @3 @3 @M.P $ \i j ->
            fromIntegral (i * 3 + j + 1) :: Double
      b = makeMatrix @3 @3 @M.P $ \i j ->
            fromIntegral (9 - (i * 3 + j)) :: Double
      c = identityMatrix @3 @M.P :: Matrix 3 3 M.P Double
      result = gemm 2.0 a b 0.5 c
  assertBool "gemm(0,0) = 60.5" $ (result ! (0,0)) ~= 60.5
  assertBool "gemm(0,1) = 48"   $ (result ! (0,1)) ~= 48
  assertBool "gemm(0,2) = 36"   $ (result ! (0,2)) ~= 36
  assertBool "gemm(1,0) = 168"  $ (result ! (1,0)) ~= 168
  assertBool "gemm(1,1) = 138.5"$ (result ! (1,1)) ~= 138.5
  assertBool "gemm(1,2) = 108"  $ (result ! (1,2)) ~= 108
  assertBool "gemm(2,0) = 276"  $ (result ! (2,0)) ~= 276
  assertBool "gemm(2,1) = 228"  $ (result ! (2,1)) ~= 228
  assertBool "gemm(2,2) = 180.5"$ (result ! (2,2)) ~= 180.5
