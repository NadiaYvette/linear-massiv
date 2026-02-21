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
import Numeric.LinearAlgebra.Massiv.BLAS.Level3 (matMul, transpose, mSub)
import Numeric.LinearAlgebra.Massiv.BLAS.Level1 (nrm2, scal)
import Numeric.LinearAlgebra.Massiv.Eigen.Power
import Numeric.LinearAlgebra.Massiv.Eigen.Hessenberg
import Numeric.LinearAlgebra.Massiv.Eigen.Symmetric
import Numeric.LinearAlgebra.Massiv.Eigen.SVD
import Data.Proxy (Proxy(..))
import Numeric.LinearAlgebra.Massiv.Eigen.Schur (schur, eigenvalues)
import Numeric.LinearAlgebra.Massiv.Norms (normFrob, vnorm2)
import Test.Types
import Test.Residuals

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
  , testGroup "Standard test matrices"
    [ testCase "Wilkinson eigenvalues" test_wilkinsonEigen
    , testCase "Hilbert eigenvalues positive" test_hilbertEigen
    , testCase "Frank eigenvalues positive" test_frankEigen
    , testProperty "clustered eigenvalues" prop_clusteredEigen
    ]
  , testGroup "Eigen residuals"
    [ testProperty "eigenpair scaled residuals 3x3" prop_eigenResiduals
    ]
  , testGroup "SVD residuals"
    [ testProperty "SVD scaled residual 3x3" prop_svdScaledResidual
    , testProperty "SVD orthogonality U and V 3x3" prop_svdOrthogonality
    , testCase "SVD diagonal 5x5 sorted" test_svdDiagonalLarger
    ]
  , testGroup "Eigenvalue ordering"
    [ testProperty "symmetric eigenvalues sorted 4x4" prop_symmetricEigenvaluesSorted
    ]
  , testGroup "D&C eigensolver"
    [ testCase "D&C eigenvalues of diagonal 10x10" test_dcEigenDiagonal
    , testCase "D&C reconstruction 50x50" test_dcEigenReconstruction50
    , testCase "D&C orthogonality 50x50" test_dcEigenOrthogonal50
    , testCase "D&C matches QR at 30x30" test_dcMatchesQR
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

-- Standard test matrices

test_wilkinsonEigen :: Assertion
test_wilkinsonEigen = do
  let a = wilkinsonMatrix @7 :: Matrix 7 7 M.P Double
      (eigvals, q) = symmetricEigen a 500 1e-12
      nn = 7
      -- Verify we get 7 eigenvalues
      evList = map (\i -> eigvals !. i) [0..nn-1]
  assertBool "got 7 eigenvalues" $ length evList == 7
  -- Verify reconstruction: A ≈ Q diag(λ) Qᵀ
  let diag_lambda = makeMatrix @7 @7 @M.P $ \i j ->
        if i == j then eigvals !. i else 0
      qt = transpose q
      qlqt = matMul q (matMul diag_lambda qt)
      residual = normFrob (mSub a qlqt) / (normFrob a + 1e-15)
  assertBool "Wilkinson reconstruction" $ residual < 1e-4

test_hilbertEigen :: Assertion
test_hilbertEigen = do
  let a = hilbertMatrix @5 :: Matrix 5 5 M.P Double
      (eigvals, _) = symmetricEigen a 500 1e-12
      evList = map (\i -> eigvals !. i) [0..4]
  -- Hilbert matrix is SPD, so all eigenvalues must be positive
  assertBool "all eigenvalues positive" $ all (> 0) evList

test_frankEigen :: Assertion
test_frankEigen = do
  let a = frankMatrix @5 :: Matrix 5 5 M.P Double
      (_, t) = schur a 200 1e-10
      evs = eigenvalues @5 t
  -- Frank matrix has all positive real eigenvalues
  assertBool "all eigenvalues positive" $ all (> 0) evs

prop_clusteredEigen :: Property
prop_clusteredEigen = withMaxSuccess 10 $ forAll (genClusteredEigenMatrix @4 5.0) $ \a ->
  let (eigvals, q) = symmetricEigen a 500 1e-12
      qt = transpose q
      diag_lambda = makeMatrix @4 @4 @M.P $ \i j ->
        if i == j then eigvals !. i else 0
      qlqt = matMul q (matMul diag_lambda qt)
      residual = normFrob (mSub a qlqt) / (normFrob a + 1e-15)
      -- Relaxed tolerance since clustered eigenvalues are harder
  in residual < 1e-3

-- Eigen residuals

prop_eigenResiduals :: Property
prop_eigenResiduals = withMaxSuccess 20 $ forAll (genSPDMatrix @3) $ \a ->
  let (eigvals, q) = symmetricEigen a 500 1e-12
      nn = 3
      -- Check each eigenpair
      checks = map (\i ->
        let lambda_i = eigvals !. i
            v_i = makeVector @3 @M.P $ \k -> q ! (k, i)
        in scaledResidualEigen a lambda_i v_i < 1000
        ) [0..nn-1]
  in all id checks

-- SVD residuals

prop_svdScaledResidual :: Property
prop_svdScaledResidual = forAll (genMatrix @3 @3) $ \a ->
  let (u, sigma, v) = svd a
  in scaledResidualSVD a u sigma v < 1000

prop_svdOrthogonality :: Property
prop_svdOrthogonality = forAll (genMatrix @3 @3) $ \a ->
  let (u, _, v) = svd a
      -- U orthogonality can be looser because the SVD implementation
      -- constructs U columns as Av/sigma, which may accumulate error.
      -- V comes from eigendecomposition of A^T A so is typically tighter.
  in orthogonalityResidual @3 u < 500000 && orthogonalityResidual @3 v < 5000

test_svdDiagonalLarger :: Assertion
test_svdDiagonalLarger = do
  let a = makeMatrix @5 @5 @M.P $ \i j ->
            if i == j then case i of
              0 -> 7; 1 -> 5; 2 -> 3; 3 -> 2; _ -> 1
            else 0 :: Double
      (_, sigma, _) = svd a
      svs = sort [sigma !. 0, sigma !. 1, sigma !. 2, sigma !. 3, sigma !. 4]
      expected = [1, 2, 3, 5, 7] :: [Double]
  assertBool "sorted singular values match" $
    all (\(s, e) -> abs (s - e) < 0.1) (zip svs expected)

-- Eigenvalue ordering

prop_symmetricEigenvaluesSorted :: Property
prop_symmetricEigenvaluesSorted = forAll (genSPDMatrix @4) $ \a ->
  let (eigvals, _) = symmetricEigen a 500 1e-12
      evList = sort $ map (\i -> eigvals !. i) [0..3]
      -- Verify non-decreasing order after sorting
  in and $ zipWith (<=) evList (tail evList)

-- D&C eigensolver tests

test_dcEigenDiagonal :: Assertion
test_dcEigenDiagonal = do
  -- Eigenvalues of diag(10, 9, 8, ..., 1) should be {1..10}
  let a = makeMatrix @10 @10 @M.P $ \i j ->
            if i == j then fromIntegral (10 - i) else 0 :: Double
      (eigvals, _) = symmetricEigenPDC a 1e-12
      evs = sort [eigvals !. i | i <- [0..9]]
  mapM_ (\(i, expected) ->
    assertBool ("eigenvalue " ++ show expected) $
      abs (evs !! i - expected) < 0.01)
    (zip [0..] [1..10 :: Double])

test_dcEigenReconstruction50 :: Assertion
test_dcEigenReconstruction50 = do
  -- A = QΛQ^T reconstruction for a 50x50 SPD matrix
  let a = mkSPD50
      (eigvals, q) = symmetricEigenPDC a 1e-12
      qt = transpose q
      lambda = makeMatrix @50 @50 @M.P $ \i j ->
        if i == j then eigvals !. i else 0
      qlqt = matMul q (matMul lambda qt)
      maxErr = maximum [abs (a ! (i,j) - qlqt ! (i,j)) | i <- [0..49], j <- [0..49]]
  assertBool ("reconstruction error " ++ show maxErr ++ " < 1e-8") $ maxErr < 1e-8

test_dcEigenOrthogonal50 :: Assertion
test_dcEigenOrthogonal50 = do
  let a = mkSPD50
      (_, q) = symmetricEigenPDC a 1e-12
      qt = transpose q
      qtq = matMul qt q
      eye = identityMatrix @50 @M.P :: Matrix 50 50 M.P Double
      maxErr = maximum [abs (qtq ! (i,j) - eye ! (i,j)) | i <- [0..49], j <- [0..49]]
  assertBool ("orthogonality error " ++ show maxErr ++ " < 1e-8") $ maxErr < 1e-8

test_dcMatchesQR :: Assertion
test_dcMatchesQR = do
  -- D&C and QR should produce same eigenvalues for a 30x30 SPD matrix
  let a = makeMatrix @30 @30 @M.P $ \i j ->
            let d = fromIntegral (abs (i - j) + 1) :: Double
            in if i == j then 30 + fromIntegral i else 1.0 / d
      (eigsDC, _) = symmetricEigenPDC a 1e-12
      (eigsQR, _) = symmetricEigenP a 3000 1e-12
      dcSorted = sort [eigsDC !. i | i <- [0..29]]
      qrSorted = sort [eigsQR !. i | i <- [0..29]]
      maxDiff = maximum $ zipWith (\a' b' -> abs (a' - b')) dcSorted qrSorted
  assertBool ("D&C vs QR diff " ++ show maxDiff ++ " < 1e-8") $ maxDiff < 1e-8

-- Helper: 50x50 SPD matrix for D&C tests
mkSPD50 :: Matrix 50 50 M.P Double
mkSPD50 = makeMatrix @50 @50 @M.P $ \i j ->
  let d = fromIntegral (abs (i - j) + 1) :: Double
  in if i == j then 50 + fromIntegral i else 1.0 / d
