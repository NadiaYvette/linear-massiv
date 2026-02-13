{-# LANGUAGE AllowAmbiguousTypes #-}

module Test.Orthogonal (orthogonalTests) where

import Test.Tasty
import Test.Tasty.QuickCheck
import Test.Tasty.HUnit
import qualified Data.Massiv.Array as M

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.BLAS.Level3 (matMul, transpose)
import Numeric.LinearAlgebra.Massiv.Orthogonal.Householder
import Numeric.LinearAlgebra.Massiv.Orthogonal.Givens
import Numeric.LinearAlgebra.Massiv.Orthogonal.QR
import Numeric.LinearAlgebra.Massiv.Orthogonal.LeastSquares
import Numeric.LinearAlgebra.Massiv.BLAS.Level2 (matvec)
import Test.Types

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
