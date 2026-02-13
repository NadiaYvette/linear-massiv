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
import Test.Types

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
