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
import Numeric.LinearAlgebra.Massiv.Eigen.SVD (singularValues)
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
    , testProperty "norm consistency with sigma_max" prop_normConsistency
    , testProperty "submultiplicativity" prop_submultiplicativity
    , testCase "known Frobenius norm" test_knownFrobNorm
    , testCase "Hilbert matrix norm relationships" test_hilbertNorms
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

-- | For 5×5 matrices: sigma_max ≤ normFrob A ≤ sqrt(5) * sigma_max.
--
-- The Frobenius norm satisfies ‖A‖_F = sqrt(sum sigma_i^2), so
-- sigma_max ≤ ‖A‖_F ≤ sqrt(n) * sigma_max.
prop_normConsistency :: Property
prop_normConsistency = forAll (genMatrix @5 @5) $ \a ->
  let sv = singularValues a
      sigmaMax = sv !. 0
      frobA = normFrob a
  in sigmaMax <= frobA + 1e-10
     && frobA <= sqrt 5 * sigmaMax + 1e-10

-- | Submultiplicativity: ‖A·B‖_F ≤ ‖A‖_F · ‖B‖_F for 3×3 matrices.
prop_submultiplicativity :: Property
prop_submultiplicativity =
  forAll ((,) <$> genMatrix @3 @3 <*> genMatrix @3 @3) $ \(a, b) ->
    normFrob (matMul a b) <= normFrob a * normFrob b + 1e-10

-- | Hilbert matrix norm relationships for hilbertMatrix @4:
-- norm1 ≤ normFrob * sqrt(n) and normFrob ≤ sqrt(n) * normInf.
test_hilbertNorms :: Assertion
test_hilbertNorms = do
  let h = hilbertMatrix @4 :: Matrix 4 4 M.P Double
      frobH = normFrob h
      n1H = norm1 h
      niH = normInf h
      sqrtN = sqrt 4 :: Double
  assertBool ("norm1 <= normFrob * sqrt(4): norm1=" ++ show n1H
              ++ " normFrob*2=" ++ show (frobH * sqrtN))
    $ n1H <= frobH * sqrtN + 1e-12
  assertBool ("normFrob <= sqrt(4) * normInf: normFrob=" ++ show frobH
              ++ " 2*normInf=" ++ show (sqrtN * niH))
    $ frobH <= sqrtN * niH + 1e-12

test_knownFrobNorm :: Assertion
test_knownFrobNorm = do
  -- ‖[[1,2],[3,4]]‖_F = sqrt(1+4+9+16) = sqrt(30)
  let a = makeMatrix @2 @2 @M.P $ \i j -> case (i,j) of
            (0,0) -> 1; (0,1) -> 2; (1,0) -> 3; (1,1) -> 4; _ -> 0 :: Double
  assertBool "Frobenius norm" $ abs (normFrob a - sqrt 30) < 1e-12
