{-# LANGUAGE AllowAmbiguousTypes #-}

module Test.Orthogonal (orthogonalTests) where

import Test.Tasty
import Test.Tasty.QuickCheck
import Test.Tasty.HUnit
import qualified Data.Massiv.Array as M

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.BLAS.Level3 (matMul, transpose, mSub)
import Numeric.LinearAlgebra.Massiv.Orthogonal.Householder
import Numeric.LinearAlgebra.Massiv.Orthogonal.Givens
import Numeric.LinearAlgebra.Massiv.Orthogonal.QR (qr, qrGivens)
import Numeric.LinearAlgebra.Massiv.Orthogonal.LeastSquares
import Numeric.LinearAlgebra.Massiv.BLAS.Level2 (matvec)
import Numeric.LinearAlgebra.Massiv.Norms (normFrob, vnorm2)
import Test.Types
import Test.Residuals

orthogonalTests :: TestTree
orthogonalTests = testGroup "Orthogonal"
  [ testGroup "Householder"
    [ testProperty "Householder zeros subdiagonal" prop_householderZeros
    , testProperty "Householder matrix orthogonal" prop_householderOrthogonal
    ]
  , testGroup "Givens"
    [ testProperty "Givens zeros target" prop_givensZeros
    , testProperty "Givens preserves norm" prop_givensNorm
    ]
  , testGroup "QR"
    [ testProperty "A = QR reconstruction" prop_qrReconstruction
    , testProperty "Q orthogonal" prop_qOrthogonal
    , testProperty "R upper triangular" prop_rUpperTriangular
    , testCase "3x3 QR known" test_qrKnown
    ]
  , testGroup "Least Squares"
    [ testCase "overdetermined system" test_leastSquares
    ]
  , testGroup "QR scaled residuals"
    [ testProperty "QR scaled residual 5x5" prop_qrScaledResidual5
    , testProperty "Q orthogonality residual 5x5" prop_qOrthogonalityResidual5
    , testProperty "Q orthogonality residual 10x10" prop_qOrthogonalityResidual10
    ]
  , testGroup "QR reconstruction (larger)"
    [ testProperty "QR reconstruction 10x10" prop_qrReconstruction10x10
    , testProperty "QR reconstruction 10x5 rectangular" prop_qrReconstruction10x5
    ]
  , testGroup "Givens vs Householder"
    [ testProperty "Givens matches Householder 5x5" prop_givensMatchesHouseholder
    ]
  , testGroup "Least squares (property)"
    [ testProperty "normal equations 10x5" prop_leastSquares10x5
    ]
  ]

-- Householder tests

prop_householderZeros :: Property
prop_householderZeros = forAll (genVector @4) $ \x ->
  let (v, beta) = householderVector x
      h = householderMatrix @4 @M.P v beta
      hx = matvec h x
      -- All entries below first should be ~0
  in all (\i -> abs (hx !. i) < 1e-6) [1..3]

prop_householderOrthogonal :: Property
prop_householderOrthogonal = forAll (genVector @3) $ \x ->
  let (v, beta) = householderVector x
      h = householderMatrix @3 @M.P v beta
      ht = transpose h
      hht = matMul h ht
      eye = identityMatrix @3 @M.P :: Matrix 3 3 M.P Double
  in matApproxEq @3 @3 hht eye

-- Givens tests

prop_givensZeros :: Property
prop_givensZeros = forAll ((,) <$> choose (-10, 10) <*> choose (-10, 10)) $ \(a, b) ->
  let (c, s) = givensRotation a b
      -- Convention: [c, -s; s, c] * [a; b] = [r; 0]
      zero = s * a + c * b  -- should be ~0
  in abs zero < (1e-10 :: Double)

prop_givensNorm :: Property
prop_givensNorm = forAll ((,) <$> choose (-10, 10) <*> choose (-10, 10)) $ \(a, b) ->
  let (c, s) = givensRotation a (b :: Double)
      r = c * a - s * b  -- r = sqrt(a² + b²)
  in abs (r * r - (a * a + b * b)) < 1e-8

-- QR tests

prop_qrReconstruction :: Property
prop_qrReconstruction = forAll (genMatrix @4 @3) $ \a ->
  let (q, r) = qr a
      qr_ = matMul q r
  in matApproxEq @4 @3 a qr_

prop_qOrthogonal :: Property
prop_qOrthogonal = forAll (genMatrix @3 @3) $ \a ->
  let (q, _) = qr a
      qt = transpose q
      qtq = matMul qt q
      eye = identityMatrix @3 @M.P :: Matrix 3 3 M.P Double
  in matApproxEq @3 @3 qtq eye

prop_rUpperTriangular :: Property
prop_rUpperTriangular = forAll (genMatrix @3 @3) $ \a ->
  let (_, r) = qr a
  in all (\(i, j) -> abs (r ! (i, j)) < 1e-8)
     [(i, j) | i <- [0..2], j <- [0..2], i > j]

test_qrKnown :: Assertion
test_qrKnown = do
  let a = makeMatrix @3 @3 @M.P $ \i j -> case (i,j) of
            (0,0) -> 12; (0,1) -> -51; (0,2) -> 4
            (1,0) -> 6;  (1,1) -> 167; (1,2) -> -68
            (2,0) -> -4; (2,1) -> 24;  (2,2) -> -41
            _ -> 0 :: Double
      (q, r) = qr a
      qr_ = matMul q r
  assertBool "QR reconstruction" $ matApproxEq @3 @3 a qr_

-- Least squares test

test_leastSquares :: Assertion
test_leastSquares = do
  -- Overdetermined system: 4 equations, 2 unknowns
  -- A = [[1,1],[1,2],[1,3],[1,4]], b = [6,5,7,10]
  -- Least squares fit: y = a + bx
  let a = makeMatrix @4 @2 @M.P $ \i j -> case (i,j) of
            (0,0) -> 1; (0,1) -> 1
            (1,0) -> 1; (1,1) -> 2
            (2,0) -> 1; (2,1) -> 3
            (3,0) -> 1; (3,1) -> 4
            _ -> 0 :: Double
      b = makeVector @4 @M.P $ \i -> case i of
            0 -> 6; 1 -> 5; 2 -> 7; _ -> 10 :: Double
      x = leastSquaresQR a b
  -- The solution should minimize ‖Ax - b‖₂
  -- With these values: x ≈ [3.5, 1.4]
  assertBool "intercept reasonable" $ abs (x !. 0 - 3.5) < 0.5
  assertBool "slope reasonable" $ abs (x !. 1 - 1.4) < 0.5

-- QR scaled residual tests

prop_qrScaledResidual5 :: Property
prop_qrScaledResidual5 = forAll (genMatrix @5 @5) $ \a ->
  let (q, r) = qr a
  in scaledResidualQR a q r < 100

prop_qOrthogonalityResidual5 :: Property
prop_qOrthogonalityResidual5 = forAll (genMatrix @5 @5) $ \a ->
  let (q, _) = qr a
  in orthogonalityResidual @5 q < 100

prop_qOrthogonalityResidual10 :: Property
prop_qOrthogonalityResidual10 = withMaxSuccess 20 $ forAll (genMatrix @10 @10) $ \a ->
  let (q, _) = qr a
  in orthogonalityResidual @10 q < 100

-- QR reconstruction (larger)

prop_qrReconstruction10x10 :: Property
prop_qrReconstruction10x10 = forAll (genMatrix @10 @10) $ \a ->
  let (q, r) = qr a
  in normFrob (mSub a (matMul q r)) / (normFrob a + 1e-15) < 1e-6

prop_qrReconstruction10x5 :: Property
prop_qrReconstruction10x5 = withMaxSuccess 20 $ forAll (genMatrix @10 @5) $ \a ->
  let (q, r) = qr a
  in normFrob (mSub a (matMul q r)) / (normFrob a + 1e-15) < 1e-6

-- Givens vs Householder

prop_givensMatchesHouseholder :: Property
prop_givensMatchesHouseholder = forAll (genMatrix @5 @5) $ \a ->
  let (q1, r1) = qr a
      (q2, r2) = qrGivens a
      a1 = matMul q1 r1
      a2 = matMul q2 r2
  in normFrob (mSub a1 a2) / (normFrob a1 + 1e-15) < 1e-6

-- Least squares (property)

prop_leastSquares10x5 :: Property
prop_leastSquares10x5 = forAll ((,) <$> genMatrix @10 @5 <*> genVector @10) $ \(a, b) ->
  let x = leastSquaresQR a b
      -- Compute Ax as a 10-vector
      ax = matvec a x
      -- Compute residual r = Ax - b
      r_vec = makeVector @10 @M.P $ \i -> (ax !. i) - (b !. i)
      -- Compute A^T r as a 5-vector
      at = transpose a
      atr = matvec at r_vec
      -- Normal equations: A^T(Ax - b) should be approximately zero
  in vnorm2 atr / (vnorm2 b + 1e-15) < 1e-4
