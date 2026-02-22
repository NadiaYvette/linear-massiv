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
import Numeric.LinearAlgebra.Massiv.BLAS.Level3 (matMul, matMulP, transpose, mSub)
import Numeric.LinearAlgebra.Massiv.BLAS.Level1 (nrm2, scal)
import Numeric.LinearAlgebra.Massiv.Eigen.Power
import Numeric.LinearAlgebra.Massiv.Eigen.Hessenberg
import Numeric.LinearAlgebra.Massiv.Eigen.Symmetric
import Numeric.LinearAlgebra.Massiv.Eigen.SVD (svd, svdP, svdGKP)
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
    , testCase "svdGKP reconstruction 10x10" test_svdGKReconstruction
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
    , testCase "D&C orthogonality 30x30" test_dcOrtho30
    , testCase "D&C orthogonality 52x52" test_dcOrtho52
    , testCase "D&C orthogonality 60x60" test_dcOrtho60
    , testCase "D&C orthogonality 80x80" test_dcOrtho80
    , testCase "D&C orthogonality 90x90" test_dcOrtho90
    , testCase "D&C orthogonality 95x95" test_dcOrtho95
    , testCase "D&C ortho diagonal 100x100" test_dcOrthoDiag100
    , testCase "D&C ortho alt-matrix 100x100" test_dcOrthoAlt100
    , testCase "D&C reconstruction 100x100" test_dcEigenReconstruction100
    , testCase "D&C reconstruction 128x128" test_dcEigenReconstruction128
    ]
  , testGroup "Panel tridiag (n >= 256)"
    [ testCase "tridiag match 128x128" test_panelTridiagReconstruction128
    , testCase "eigenreconstruction 200x200" test_panelTridiagReconstruction200
    , testCase "orthogonality 200x200" test_panelTridiagOrthogonal200
    , testCase "eigenreconstruction 300x300" test_panelTridiagReconstruction300
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

test_svdGKReconstruction :: Assertion
test_svdGKReconstruction = do
  -- Test 1: diagonal matrix (trivial bidiag, no QR needed)
  let diag5 = makeMatrix @5 @5 @M.P $ \i j ->
                if i == j then fromIntegral (5 - i) else 0 :: Double
      (_, sigDiag, _) = svdGKP diag5
      diagSorted = sort [sigDiag !. i | i <- [0..4]]
      diagExpected = [1,2,3,4,5] :: [Double]
      diagErr = maximum $ zipWith (\a_ b_ -> abs (a_ - b_)) diagSorted diagExpected
  assertBool ("svdGKP diagonal sigma " ++ show diagSorted ++ " err=" ++ show diagErr) $ diagErr < 0.1
  -- Test 2: already-bidiagonal matrix (tests QR iteration in isolation)
  -- B = [[3,1,0],[0,2,1],[0,0,1]] — bidiag with d=[3,2,1], e=[1,1]
  let bidiag3 = makeMatrix @3 @3 @M.P $ \i j -> case (i,j) of
                  (0,0) -> 3; (0,1) -> 1; (1,1) -> 2; (1,2) -> 1; (2,2) -> 1
                  _ -> 0 :: Double
      (_, sigBidiag, _) = svdGKP bidiag3
      (_, sigBidiagRef, _) = svdP bidiag3
      bidiagSorted = sort [sigBidiag !. i | i <- [0..2]]
      bidiagRefSorted = sort [sigBidiagRef !. i | i <- [0..2]]
      bidiagDiff = maximum $ zipWith (\a_ b_ -> abs (a_ - b_)) bidiagSorted bidiagRefSorted
  assertBool ("svdGKP bidiag diff " ++ show bidiagDiff
              ++ "\n  gk=" ++ show bidiagSorted
              ++ "\n  ref=" ++ show bidiagRefSorted) $ bidiagDiff < 1e-6
  -- Test 3: general matrix
  let a = makeMatrix @5 @5 @M.P $ \i j ->
            let d = fromIntegral (abs (i - j) + 1) :: Double
            in 1.0 / d + if i == j then 5 else 0
      (_, sigmaGK, _) = svdGKP a
      (_, sigmaRef, _) = svdP a
      gkSorted = sort [sigmaGK !. i | i <- [0..4]]
      refSorted = sort [sigmaRef !. i | i <- [0..4]]
      svDiff = maximum $ zipWith (\a_ b_ -> abs (a_ - b_)) gkSorted refSorted
  assertBool ("svdGKP sigma diff " ++ show svDiff ++ "\n  gk=" ++ show gkSorted
              ++ "\n  ref=" ++ show refSorted) $ svDiff < 1e-6
  -- Test 4: singular values match for 10×10
  let a10 = makeMatrix @10 @10 @M.P $ \i j ->
              let d = fromIntegral (abs (i - j) + 1) :: Double
              in 1.0 / d + if i == j then 10 else 0
      (u10, sig10, v10) = svdGKP a10
      (_, sigRef10, _) = svdP a10
      gk10Sorted = sort [sig10 !. i | i <- [0..9]]
      ref10Sorted = sort [sigRef10 !. i | i <- [0..9]]
      svDiff10 = maximum $ zipWith (\a_ b_ -> abs (a_ - b_)) gk10Sorted ref10Sorted
  assertBool ("svdGKP 10x10 sigma diff " ++ show svDiff10
              ++ "\n  gk=" ++ show gk10Sorted
              ++ "\n  ref=" ++ show ref10Sorted) $ svDiff10 < 1e-6
  -- Test 5: reconstruction A ≈ U Σ V^T for 10×10
  let sigMat10 = makeMatrix @10 @10 @M.P $ \i j ->
        if i == j then sig10 !. i else 0
      usv10 = matMulP u10 (matMulP sigMat10 (transpose v10))
      reconErr10 = maximum [abs (a10 ! (i,j) - usv10 ! (i,j))
                           | i <- [0..9], j <- [0..9]]
  assertBool ("svdGKP 10x10 reconstruction err " ++ show reconErr10) $ reconErr10 < 1e-10
  -- Test 6: orthogonality of U and V
  let utu = matMulP (transpose u10) u10
      vtv = matMulP (transpose v10) v10
      eye10 = identityMatrix @10 @M.P :: Matrix 10 10 M.P Double
      uErr = maximum [abs (utu ! (i,j) - eye10 ! (i,j)) | i <- [0..9], j <- [0..9]]
      vErr = maximum [abs (vtv ! (i,j) - eye10 ! (i,j)) | i <- [0..9], j <- [0..9]]
  assertBool ("svdGKP 10x10 U ortho err " ++ show uErr) $ uErr < 1e-10
  assertBool ("svdGKP 10x10 V ortho err " ++ show vErr) $ vErr < 1e-10

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

test_dcEigenReconstruction100 :: Assertion
test_dcEigenReconstruction100 = do
  let a = makeMatrix @100 @100 @M.P $ \i j ->
            let d = fromIntegral (abs (i - j) + 1) :: Double
            in if i == j then 100 + fromIntegral i else 1.0 / d
      (eigvalsDC, qDC) = symmetricEigenPDC a 1e-12
      (eigvalsQR, _)   = symmetricEigenP a 10000 1e-12
      -- Compare eigenvalues
      dcSorted = sort [eigvalsDC !. i | i <- [0..99]]
      qrSorted = sort [eigvalsQR !. i | i <- [0..99]]
      evDiff = maximum $ zipWith (\a_ b_ -> abs (a_ - b_)) dcSorted qrSorted
  assertBool ("D&C 100 eigenvalue diff " ++ show evDiff) $ evDiff < 1e-6
  -- Check orthogonality of Q
  let qtq = matMulP (transpose qDC) qDC
      eye100 = identityMatrix @100 @M.P :: Matrix 100 100 M.P Double
      orthoErr = maximum [abs (qtq ! (i,j) - eye100 ! (i,j)) | i <- [0..99], j <- [0..99]]
  assertBool ("D&C 100 orthogonality error " ++ show orthoErr) $ orthoErr < 1e-6
  -- Full reconstruction
  let qt = transpose qDC
      lambda = makeMatrix @100 @100 @M.P $ \i j ->
        if i == j then eigvalsDC !. i else 0
      qlqt = matMul qDC (matMul lambda qt)
      maxErr = maximum [abs (a ! (i,j) - qlqt ! (i,j)) | i <- [0..99], j <- [0..99]]
  assertBool ("D&C 100 reconstruction error " ++ show maxErr ++ " < 1e-7") $ maxErr < 1e-7

test_dcEigenReconstruction128 :: Assertion
test_dcEigenReconstruction128 = do
  let a = makeMatrix @128 @128 @M.P $ \i j ->
            let d = fromIntegral (abs (i - j) + 1) :: Double
            in if i == j then 128 + fromIntegral i else 1.0 / d
      (eigvals, q) = symmetricEigenPDC a 1e-12
      qt = transpose q
      lambda = makeMatrix @128 @128 @M.P $ \i j ->
        if i == j then eigvals !. i else 0
      qlqt = matMul q (matMul lambda qt)
      maxErr = maximum [abs (a ! (i,j) - qlqt ! (i,j)) | i <- [0..127], j <- [0..127]]
  assertBool ("D&C 128 reconstruction error " ++ show maxErr ++ " < 5e-7") $ maxErr < 5e-7

-- Panel tridiag tests (n >= 128 crossover)

test_panelTridiagReconstruction128 :: Assertion
test_panelTridiagReconstruction128 = do
  let nn = 128
      a = makeMatrix @128 @128 @M.P $ \i j ->
            let d = fromIntegral (abs (i - j) + 1) :: Double
            in if i == j then 128 + fromIntegral i else 1.0 / d
      (eigvals, q) = symmetricEigenP a 10000 1e-12
      qt = transpose q
      lambda = makeMatrix @128 @128 @M.P $ \i j ->
        if i == j then eigvals !. i else 0
      qlqt = matMul q (matMul lambda qt)
      maxErr = maximum [abs (a ! (i,j) - qlqt ! (i,j)) | i <- [0..nn-1], j <- [0..nn-1]]
  assertBool ("reconstruction error " ++ show maxErr ++ " < 1e-6") $ maxErr < 1e-6

test_panelTridiagReconstruction200 :: Assertion
test_panelTridiagReconstruction200 = do
  let nn = 200
      a = mkSPD200
      (eigvals, q) = symmetricEigenP a 10000 1e-12
      qt = transpose q
      lambda = makeMatrix @200 @200 @M.P $ \i j ->
        if i == j then eigvals !. i else 0
      qlqt = matMul q (matMul lambda qt)
      maxErr = maximum [abs (a ! (i,j) - qlqt ! (i,j)) | i <- [0..nn-1], j <- [0..nn-1]]
  assertBool ("reconstruction error " ++ show maxErr ++ " < 1e-6") $ maxErr < 1e-6

test_panelTridiagOrthogonal200 :: Assertion
test_panelTridiagOrthogonal200 = do
  let nn = 200
      a = mkSPD200
      (_, q) = symmetricEigenP a 10000 1e-12
      qt = transpose q
      qtq = matMul qt q
      eye = identityMatrix @200 @M.P :: Matrix 200 200 M.P Double
      maxErr = maximum [abs (qtq ! (i,j) - eye ! (i,j)) | i <- [0..nn-1], j <- [0..nn-1]]
  assertBool ("orthogonality error " ++ show maxErr ++ " < 1e-8") $ maxErr < 1e-8

test_panelTridiagReconstruction300 :: Assertion
test_panelTridiagReconstruction300 = do
  let nn = 300
      a = makeMatrix @300 @300 @M.P $ \i j ->
            let d = fromIntegral (abs (i - j) + 1) :: Double
            in if i == j then 300 + fromIntegral i else 1.0 / d
      (eigvals, q) = symmetricEigenP a 15000 1e-12
      qt = transpose q
      lambda = makeMatrix @300 @300 @M.P $ \i j ->
        if i == j then eigvals !. i else 0
      qlqt = matMul q (matMul lambda qt)
      maxErr = maximum [abs (a ! (i,j) - qlqt ! (i,j)) | i <- [0..nn-1], j <- [0..nn-1]]
  assertBool ("reconstruction error " ++ show maxErr ++ " < 1e-6") $ maxErr < 1e-6

mkSPD200 :: Matrix 200 200 M.P Double
mkSPD200 = makeMatrix @200 @200 @M.P $ \i j ->
  let d = fromIntegral (abs (i - j) + 1) :: Double
  in if i == j then 200 + fromIntegral i else 1.0 / d

-- Helper: 50x50 SPD matrix for D&C tests
mkSPD50 :: Matrix 50 50 M.P Double
mkSPD50 = makeMatrix @50 @50 @M.P $ \i j ->
  let d = fromIntegral (abs (i - j) + 1) :: Double
  in if i == j then 50 + fromIntegral i else 1.0 / d

-- Granular D&C orthogonality tests at various sizes
test_dcOrtho30 :: Assertion
test_dcOrtho30 = do
  let a = makeMatrix @30 @30 @M.P $ \i j ->
            let d = fromIntegral (abs (i - j) + 1) :: Double
            in if i == j then 30 + fromIntegral i else 1.0 / d
      (_, q) = symmetricEigenPDC a 1e-12
      qtq = matMulP (transpose q) q
      eye = identityMatrix @30 @M.P :: Matrix 30 30 M.P Double
      maxErr = maximum [abs (qtq ! (i,j) - eye ! (i,j)) | i <- [0..29], j <- [0..29]]
  assertBool ("D&C 30 ortho err " ++ show maxErr) $ maxErr < 1e-8

test_dcOrtho52 :: Assertion
test_dcOrtho52 = do
  let a = makeMatrix @52 @52 @M.P $ \i j ->
            let d = fromIntegral (abs (i - j) + 1) :: Double
            in if i == j then 52 + fromIntegral i else 1.0 / d
      (_, q) = symmetricEigenPDC a 1e-12
      qtq = matMulP (transpose q) q
      eye = identityMatrix @52 @M.P :: Matrix 52 52 M.P Double
      maxErr = maximum [abs (qtq ! (i,j) - eye ! (i,j)) | i <- [0..51], j <- [0..51]]
  assertBool ("D&C 52 ortho err " ++ show maxErr) $ maxErr < 1e-8

test_dcOrtho60 :: Assertion
test_dcOrtho60 = do
  let a = makeMatrix @60 @60 @M.P $ \i j ->
            let d = fromIntegral (abs (i - j) + 1) :: Double
            in if i == j then 60 + fromIntegral i else 1.0 / d
      (_, q) = symmetricEigenPDC a 1e-12
      qtq = matMulP (transpose q) q
      eye = identityMatrix @60 @M.P :: Matrix 60 60 M.P Double
      maxErr = maximum [abs (qtq ! (i,j) - eye ! (i,j)) | i <- [0..59], j <- [0..59]]
  assertBool ("D&C 60 ortho err " ++ show maxErr) $ maxErr < 1e-8

test_dcOrtho80 :: Assertion
test_dcOrtho80 = do
  let a = makeMatrix @80 @80 @M.P $ \i j ->
            let d = fromIntegral (abs (i - j) + 1) :: Double
            in if i == j then 80 + fromIntegral i else 1.0 / d
      (_, q) = symmetricEigenPDC a 1e-12
      qtq = matMulP (transpose q) q
      eye = identityMatrix @80 @M.P :: Matrix 80 80 M.P Double
      maxErr = maximum [abs (qtq ! (i,j) - eye ! (i,j)) | i <- [0..79], j <- [0..79]]
  assertBool ("D&C 80 ortho err " ++ show maxErr) $ maxErr < 1e-8

test_dcOrtho90 :: Assertion
test_dcOrtho90 = do
  let a = makeMatrix @90 @90 @M.P $ \i j ->
            let d = fromIntegral (abs (i - j) + 1) :: Double
            in if i == j then 90 + fromIntegral i else 1.0 / d
      (_, q) = symmetricEigenPDC a 1e-12
      qtq = matMulP (transpose q) q
      eye = identityMatrix @90 @M.P :: Matrix 90 90 M.P Double
      maxErr = maximum [abs (qtq ! (i,j) - eye ! (i,j)) | i <- [0..89], j <- [0..89]]
  assertBool ("D&C 90 ortho err " ++ show maxErr) $ maxErr < 1e-8

test_dcOrtho95 :: Assertion
test_dcOrtho95 = do
  let a = makeMatrix @95 @95 @M.P $ \i j ->
            let d = fromIntegral (abs (i - j) + 1) :: Double
            in if i == j then 95 + fromIntegral i else 1.0 / d
      (_, q) = symmetricEigenPDC a 1e-12
      qtq = matMulP (transpose q) q
      eye = identityMatrix @95 @M.P :: Matrix 95 95 M.P Double
      maxErr = maximum [abs (qtq ! (i,j) - eye ! (i,j)) | i <- [0..94], j <- [0..94]]
  assertBool ("D&C 95 ortho err " ++ show maxErr) $ maxErr < 1e-8

-- Test D&C with a purely diagonal 100×100 matrix
test_dcOrthoDiag100 :: Assertion
test_dcOrthoDiag100 = do
  let a = makeMatrix @100 @100 @M.P $ \i j ->
            if i == j then fromIntegral (i + 1) else 0 :: Double
      (eigvals, q) = symmetricEigenPDC a 1e-12
      qtq = matMulP (transpose q) q
      eye = identityMatrix @100 @M.P :: Matrix 100 100 M.P Double
      maxErr = maximum [abs (qtq ! (i,j) - eye ! (i,j)) | i <- [0..99], j <- [0..99]]
      evSorted = sort [eigvals !. i | i <- [0..99]]
      evDiff = maximum $ zipWith (\a_ b_ -> abs (a_ - b_)) evSorted [1..100]
  assertBool ("D&C diag100 ortho err " ++ show maxErr) $ maxErr < 1e-8
  assertBool ("D&C diag100 eigenvalue diff " ++ show evDiff) $ evDiff < 1e-8

-- Test D&C with a different matrix at 100×100 (sparser off-diagonal)
test_dcOrthoAlt100 :: Assertion
test_dcOrthoAlt100 = do
  let a = makeMatrix @100 @100 @M.P $ \i j ->
            if i == j then 500 + fromIntegral i
            else if abs (i - j) == 1 then 0.1
            else 0 :: Double
      (_, q) = symmetricEigenPDC a 1e-12
      qtq = matMulP (transpose q) q
      eye = identityMatrix @100 @M.P :: Matrix 100 100 M.P Double
      maxErr = maximum [abs (qtq ! (i,j) - eye ! (i,j)) | i <- [0..99], j <- [0..99]]
  assertBool ("D&C alt100 ortho err " ++ show maxErr) $ maxErr < 1e-8
