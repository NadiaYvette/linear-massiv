{-# LANGUAGE AllowAmbiguousTypes #-}

module Test.Eigen (eigenTests) where

import Test.Tasty
import Test.Tasty.QuickCheck
import Test.Tasty.HUnit
import qualified Data.Massiv.Array as M
import Data.List (sort)

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.BLAS.Level2 (matvec)
import Numeric.LinearAlgebra.Massiv.BLAS.Level3 (matMul, transpose)
import Numeric.LinearAlgebra.Massiv.BLAS.Level1 (nrm2, scal)
import Numeric.LinearAlgebra.Massiv.Eigen.Power
import Numeric.LinearAlgebra.Massiv.Eigen.Hessenberg
import Numeric.LinearAlgebra.Massiv.Eigen.Symmetric
import Numeric.LinearAlgebra.Massiv.Eigen.SVD
import Test.Types

eigenTests :: TestTree
eigenTests = testGroup "Eigenvalue"
  [ testGroup "Power Method"
    [ testCase "dominant eigenvalue of diagonal" test_powerDiagonal
    ]
  , testGroup "Hessenberg"
    [ testProperty "A = QHQᵀ reconstruction" prop_hessenbergReconstruction
    , testProperty "H is upper Hessenberg" prop_hessenbergForm
    ]
  , testGroup "Symmetric"
    [ testProperty "A = QΛQᵀ reconstruction" prop_symmetricEigenReconstruction
    , testProperty "Q orthogonal" prop_symmetricQOrthogonal
    , testCase "eigenvalues of diagonal" test_symmetricDiagonal
    ]
  , testGroup "Jacobi"
    [ testCase "jacobi eigenvalues of known matrix" test_jacobiKnown
    ]
  , testGroup "SVD"
    [ testProperty "A ≈ UΣVᵀ reconstruction" prop_svdReconstruction
    , testCase "singular values of diagonal" test_svdDiagonal
    ]
  ]

-- Power method

test_powerDiagonal :: Assertion
test_powerDiagonal = do
  -- A = diag(3, 2, 1) → dominant eigenvalue = 3
  let a = makeMatrix @3 @3 @M.P $ \i j ->
            if i == j then case i of { 0 -> 3; 1 -> 2; _ -> 1 } else 0 :: Double
      q0 = makeVector @3 @M.P $ \_ -> 1 / sqrt 3 :: Double
      (lam, _) = powerMethod a q0 100 1e-10
  assertBool "eigenvalue ~ 3" $ abs (lam - 3) < 0.01

-- Hessenberg

prop_hessenbergReconstruction :: Property
prop_hessenbergReconstruction = forAll (genMatrix @4 @4) $ \a ->
  let (q, h) = hessenberg a
      qt = transpose q
      qhqt = matMul q (matMul h qt)
  in matApproxEq @4 @4 a qhqt

prop_hessenbergForm :: Property
prop_hessenbergForm = forAll (genMatrix @4 @4) $ \a ->
  let (_, h) = hessenberg a
  in all (\(i, j) -> abs (h ! (i, j)) < 1e-8)
     [(i, j) | i <- [0..3], j <- [0..3], i > j + 1]

-- Symmetric eigenvalue

prop_symmetricEigenReconstruction :: Property
prop_symmetricEigenReconstruction = forAll (genSPDMatrix @3) $ \a ->
  let (eigvals, q) = symmetricEigen a 500 1e-12
      qt = transpose q
      -- Reconstruct: A ≈ Q diag(λ) Qᵀ
      lambda = makeMatrix @3 @3 @M.P $ \i j ->
        if i == j then eigvals !. i else 0
      qlqt = matMul q (matMul lambda qt)
      -- Use relaxed tolerance for iterative eigenvalue decomposition
  in all (\(i, j) -> abs (a ! (i,j) - qlqt ! (i,j)) < 1e-4 * (1 + abs (a ! (i,j))))
     [(i, j) | i <- [0..2], j <- [0..2]]

prop_symmetricQOrthogonal :: Property
prop_symmetricQOrthogonal = forAll (genSPDMatrix @3) $ \a ->
  let (_, q) = symmetricEigen a 500 1e-12
      qt = transpose q
      qtq = matMul qt q
      eye = identityMatrix @3 @M.P :: Matrix 3 3 M.P Double
  in matApproxEq @3 @3 qtq eye

test_symmetricDiagonal :: Assertion
test_symmetricDiagonal = do
  -- Eigenvalues of diag(5, 3, 1) should be {1, 3, 5}
  let a = makeMatrix @3 @3 @M.P $ \i j ->
            if i == j then case i of { 0 -> 5; 1 -> 3; _ -> 1 } else 0 :: Double
      (eigvals, _) = symmetricEigen a 100 1e-12
      evs = sort [eigvals !. 0, eigvals !. 1, eigvals !. 2]
  assertBool "eigenvalue 1" $ abs (evs !! 0 - 1) < 0.01
  assertBool "eigenvalue 3" $ abs (evs !! 1 - 3) < 0.01
  assertBool "eigenvalue 5" $ abs (evs !! 2 - 5) < 0.01

-- Jacobi

test_jacobiKnown :: Assertion
test_jacobiKnown = do
  let a = makeMatrix @3 @3 @M.P $ \i j -> case (i,j) of
            (0,0) -> 4; (0,1) -> 1; (0,2) -> 0
            (1,0) -> 1; (1,1) -> 3; (1,2) -> 1
            (2,0) -> 0; (2,1) -> 1; (2,2) -> 2
            _ -> 0 :: Double
      (eigvals, q) = jacobiEigen a 100 1e-12
      qt = transpose q
      lambda = makeMatrix @3 @3 @M.P $ \i j ->
        if i == j then eigvals !. i else 0
      qlqt = matMul q (matMul lambda qt)
  assertBool "Jacobi reconstruction" $ matApproxEq @3 @3 a qlqt

-- SVD

prop_svdReconstruction :: Property
prop_svdReconstruction = forAll (genMatrix @3 @3) $ \a ->
  let (u, sigma, v) = svd a
      vt = transpose v
      -- Reconstruct: U * diag(σ) * Vᵀ
      sigMat = makeMatrix @3 @3 @M.P $ \i j ->
        if i == j then sigma !. i else 0
      usv = matMul u (matMul sigMat vt)
  in matApproxEq @3 @3 a usv

test_svdDiagonal :: Assertion
test_svdDiagonal = do
  -- SVD of diag(5, 3, 1) should give singular values {5, 3, 1}
  let a = makeMatrix @3 @3 @M.P $ \i j ->
            if i == j then case i of { 0 -> 5; 1 -> 3; _ -> 1 } else 0 :: Double
      (_, sigma, _) = svd a
      svs = sort [sigma !. 0, sigma !. 1, sigma !. 2]
  assertBool "sv 1" $ abs (svs !! 0 - 1) < 0.1
  assertBool "sv 3" $ abs (svs !! 1 - 3) < 0.1
  assertBool "sv 5" $ abs (svs !! 2 - 5) < 0.1
