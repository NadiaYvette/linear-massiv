{-# LANGUAGE AllowAmbiguousTypes #-}

module Test.Solve (solveTests) where

import Test.Tasty
import Test.Tasty.QuickCheck
import Test.Tasty.HUnit
import qualified Data.Massiv.Array as M

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.BLAS.Level2 (matvec)
import Numeric.LinearAlgebra.Massiv.BLAS.Level3 (matMul, transpose)
import Numeric.LinearAlgebra.Massiv.Solve.Triangular
import Numeric.LinearAlgebra.Massiv.Solve.LU
import Numeric.LinearAlgebra.Massiv.Solve.Cholesky
import Numeric.LinearAlgebra.Massiv.Solve.Banded (tridiagSolve)
import Numeric.LinearAlgebra.Massiv.Eigen.Symmetric (symmetricEigen)
import Test.Types
import Test.Residuals

solveTests :: TestTree
solveTests = testGroup "Solve"
  [ testGroup "Triangular"
    [ testProperty "forward sub: Lx = b roundtrip" prop_forwardSubRoundtrip
    , testProperty "back sub: Ux = b roundtrip" prop_backSubRoundtrip
    ]
  , testGroup "LU"
    [ testProperty "PA = LU reconstruction" prop_luReconstruction
    , testProperty "luSolve: Ax = b roundtrip" prop_luSolveRoundtrip
    , testCase "3x3 LU known" test_luKnown
    ]
  , testGroup "Cholesky"
    [ testProperty "A = GGᵀ reconstruction" prop_choleskyReconstruction
    , testProperty "choleskySolve roundtrip" prop_choleskySolveRoundtrip
    ]
  , testGroup "Tridiagonal"
    [ testCase "tridiag solve known" test_tridiagKnown
    ]
  , testGroup "Ill-conditioned"
    [ testCase "LU on Hilbert 5×5" test_luHilbert
    , testProperty "Cholesky near-singular SPD" prop_choleskyNearSingular
    ]
  , testGroup "Randsvd"
    [ testProperty "luSolve with randsvd matrix" prop_luSolveRandsvd
    , testProperty "choleskySolve with randsvd SPD" prop_choleskySolveRandsvd
    ]
  , testGroup "Scaled Residuals"
    [ testProperty "LU residual bound 3×3" prop_luResidualBound
    , testProperty "Cholesky residual bound 3×3" prop_choleskyResidualBound
    ]
  , testGroup "Cross-module"
    [ testCase "det ≈ product of eigenvalues" test_detEqualsEigenProduct
    ]
  , testGroup "Larger Sizes"
    [ testProperty "LU reconstruction 5×5" prop_luReconstruction5
    , testProperty "Cholesky reconstruction 5×5" prop_choleskyReconstruction5
    , testProperty "luSolve 10×10" prop_luSolve10
    ]
  ]

-- Triangular tests

prop_forwardSubRoundtrip :: Property
prop_forwardSubRoundtrip = forAll ((,) <$> genLowerTriangular @3 <*> genVector @3) $ \(l, b) ->
  let x = forwardSub l b
      b' = matvec l x
  in vecApproxEq @3 b b'

prop_backSubRoundtrip :: Property
prop_backSubRoundtrip = forAll ((,) <$> genUpperTriangular @3 <*> genVector @3) $ \(u, b) ->
  let x = backSub u b
      b' = matvec u x
  in vecApproxEq @3 b b'

-- LU tests

prop_luReconstruction :: Property
prop_luReconstruction = forAll (genMatrix @3 @3) $ \a ->
  let (luMat, pivArr) = lu a
      nn = 3 :: Int
      l = makeMatrix @3 @3 @M.P $ \i j ->
        if i == j then 1
        else if i > j then luMat ! (i, j)
        else 0 :: Double
      u = makeMatrix @3 @3 @M.P $ \i j ->
        if i <= j then luMat ! (i, j)
        else 0 :: Double
      lu_ = matMul l u
      -- PA
      pa = makeMatrix @3 @3 @M.P $ \i j ->
        a ! (M.index' pivArr i, j)
  in matApproxEq @3 @3 pa lu_

prop_luSolveRoundtrip :: Property
prop_luSolveRoundtrip = forAll ((,) <$> genMatrix @3 @3 <*> genVector @3) $ \(a, b) ->
  -- Only test for non-singular matrices (skip near-singular)
  let (_, pivArr) = lu a
      luMat = fst (lu a)
      -- Check diagonal of U
      diagOk = all (\i -> abs (luMat ! (i, i)) > 1e-6) [0..2]
  in diagOk ==>
    let x = luSolve a b
        b' = matvec a x
    in vecApproxEq @3 b b'

test_luKnown :: Assertion
test_luKnown = do
  -- A = [[2,1,1],[4,3,3],[8,7,9]]
  let a = makeMatrix @3 @3 @M.P $ \i j -> case (i,j) of
            (0,0) -> 2; (0,1) -> 1; (0,2) -> 1
            (1,0) -> 4; (1,1) -> 3; (1,2) -> 3
            (2,0) -> 8; (2,1) -> 7; (2,2) -> 9
            _ -> 0 :: Double
      b = makeVector @3 @M.P $ \i -> case i of
            0 -> 1; 1 -> 1; _ -> 1 :: Double
      x = luSolve a b
      b' = matvec a x
  assertBool "LU solve roundtrip" $ vecApproxEq @3 b b'

-- Cholesky tests

prop_choleskyReconstruction :: Property
prop_choleskyReconstruction = forAll (genSPDMatrix @3) $ \a ->
  let g = cholesky a
      gt = transpose g
      ggt = matMul g gt
  in matApproxEq @3 @3 a ggt

prop_choleskySolveRoundtrip :: Property
prop_choleskySolveRoundtrip = forAll ((,) <$> genSPDMatrix @3 <*> genVector @3) $ \(a, b) ->
  let x = choleskySolve a b
      b' = matvec a x
  in vecApproxEq @3 b b'

-- Tridiagonal test

test_tridiagKnown :: Assertion
test_tridiagKnown = do
  -- Tridiagonal SPD: A = [[2,-1,0],[-1,2,-1],[0,-1,2]]
  -- b = [1, 0, 1]
  -- x should be [1, 1, 1]
  let diag_ = makeVector @3 @M.P $ \_ -> 2 :: Double
      supdiag = makeVector @3 @M.P $ \i -> if i < 2 then -1 else 0 :: Double
      b = makeVector @3 @M.P $ \i -> case i of { 0 -> 1; 1 -> 0; _ -> 1 } :: Double
      x = tridiagSolve diag_ supdiag b
      -- Verify Ax = b manually: 2*1 + (-1)*1 = 1, (-1)*1 + 2*1 + (-1)*1 = 0, (-1)*1 + 2*1 = 1
  assertBool "x(0) ~ 1" $ (x !. 0) ~= 1
  assertBool "x(1) ~ 1" $ (x !. 1) ~= 1
  assertBool "x(2) ~ 1" $ (x !. 2) ~= 1

------------------------------------------------------------------------
-- Ill-conditioned tests
------------------------------------------------------------------------

-- | Test 1: LU solve on the 5×5 Hilbert matrix (very ill-conditioned).
test_luHilbert :: Assertion
test_luHilbert = do
  let h = hilbertMatrix @5
      x = makeVector @5 @M.P $ \i -> fromIntegral (i + 1)
      b = matvec h x
      x' = luSolve h b
      residual = scaledResidualLinear @5 h x' b
  assertBool ("LU Hilbert residual too large: " ++ show residual)
    (residual < 1e-4)

-- | Test 2: Cholesky on near-singular SPD matrix (cond = 1e4).
prop_choleskyNearSingular :: Property
prop_choleskyNearSingular = withMaxSuccess 20 $
  forAll (genSPDMatrixWithCond @5 1e4) $ \a ->
    let g = cholesky a
        residual = scaledResidualCholesky @5 a g
    in counterexample ("Cholesky residual = " ++ show residual) $
       residual < 1000

------------------------------------------------------------------------
-- Randsvd tests
------------------------------------------------------------------------

-- | Test 3: LU solve with a randsvd matrix (cond = 100).
prop_luSolveRandsvd :: Property
prop_luSolveRandsvd = withMaxSuccess 30 $
  forAll ((,) <$> genMatrixWithCond @5 @5 100.0 <*> genVector @5) $ \(a, b) ->
    let x = luSolve a b
        residual = scaledResidualLinear @5 a x b
    in counterexample ("LU randsvd residual = " ++ show residual) $
       residual < 0.01

-- | Test 4: Cholesky solve with a randsvd SPD matrix (cond = 100).
prop_choleskySolveRandsvd :: Property
prop_choleskySolveRandsvd = withMaxSuccess 20 $
  forAll ((,) <$> genSPDMatrixWithCond @5 100.0 <*> genVector @5) $ \(a, b) ->
    let x = choleskySolve a b
        residual = scaledResidualLinear @5 a x b
    in counterexample ("Cholesky randsvd residual = " ++ show residual) $
       residual < 0.01

------------------------------------------------------------------------
-- Scaled residual bound tests
------------------------------------------------------------------------

-- | Test 5: LU residual bound for well-conditioned 3×3 random matrices.
prop_luResidualBound :: Property
prop_luResidualBound = forAll ((,) <$> genMatrix @3 @3 <*> genVector @3) $ \(a, b) ->
  let (luMat, _) = lu a
      diagOk = all (\i -> abs (luMat ! (i, i)) > 1e-6) [0..2]
  in diagOk ==>
    let x = luSolve a b
        residual = scaledResidualLinear @3 a x b
    in counterexample ("LU residual = " ++ show residual) $
       residual < 1e-6

-- | Test 6: Cholesky residual bound for well-conditioned 3×3 SPD matrices.
prop_choleskyResidualBound :: Property
prop_choleskyResidualBound = forAll ((,) <$> genSPDMatrix @3 <*> genVector @3) $ \(a, b) ->
  let x = choleskySolve a b
      residual = scaledResidualLinear @3 a x b
  in counterexample ("Cholesky residual = " ++ show residual) $
     residual < 1e-6

------------------------------------------------------------------------
-- Cross-module tests
------------------------------------------------------------------------

-- | Test 7: Determinant equals the product of eigenvalues for a known SPD matrix.
--
-- We use a diagonal-dominant matrix: diag(4,3,2,1) + 0.1*ones(4,4).
test_detEqualsEigenProduct :: Assertion
test_detEqualsEigenProduct = do
  let a = makeMatrix @4 @4 @M.P $ \i j ->
            let diag_ = case i of { 0 -> 4; 1 -> 3; 2 -> 2; _ -> 1 } :: Double
            in (if i == j then diag_ else 0) + 0.1
      d = det a
      (eigenvals, _) = symmetricEigen a 200 1e-12
      eigenProd = product [ eigenvals !. i | i <- [0..3] ]
      relErr = abs (d - eigenProd) / (abs d + 1e-15)
  assertBool ("det vs eigenproduct relative error too large: " ++ show relErr)
    (relErr < 1e-4)

------------------------------------------------------------------------
-- Larger-size tests
------------------------------------------------------------------------

-- | Test 8: LU reconstruction at 5×5 (like existing 3×3 test).
prop_luReconstruction5 :: Property
prop_luReconstruction5 = forAll (genMatrix @5 @5) $ \a ->
  let (luMat, pivArr) = lu a
      l = makeMatrix @5 @5 @M.P $ \i j ->
        if i == j then 1
        else if i > j then luMat ! (i, j)
        else 0 :: Double
      u = makeMatrix @5 @5 @M.P $ \i j ->
        if i <= j then luMat ! (i, j)
        else 0 :: Double
      lu_ = matMul l u
      pa = makeMatrix @5 @5 @M.P $ \i j ->
        a ! (M.index' pivArr i, j)
  in matApproxEq @5 @5 pa lu_

-- | Test 9: Cholesky reconstruction at 5×5 (like existing 3×3 test).
prop_choleskyReconstruction5 :: Property
prop_choleskyReconstruction5 = forAll (genSPDMatrix @5) $ \a ->
  let g = cholesky a
      gt = transpose g
      ggt = matMul g gt
  in matApproxEq @5 @5 a ggt

-- | Test 10: LU solve at 10×10 with diagonal guard.
prop_luSolve10 :: Property
prop_luSolve10 = withMaxSuccess 20 $
  forAll ((,) <$> genMatrix @10 @10 <*> genVector @10) $ \(a, b) ->
    let (luMat, _) = lu a
        diagOk = all (\i -> abs (luMat ! (i, i)) > 1e-6) [0..9]
    in diagOk ==>
      let x = luSolve a b
          residual = scaledResidualLinear @10 a x b
      in counterexample ("LU 10x10 residual = " ++ show residual) $
         residual < 1e-6
