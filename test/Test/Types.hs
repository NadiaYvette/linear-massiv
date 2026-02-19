{-# LANGUAGE AllowAmbiguousTypes #-}

-- | Test helpers: Arbitrary instances, approximate equality, matrix generators.
module Test.Types
  ( -- * Approximate equality
    (~=)
  , matApproxEq
  , vecApproxEq
    -- * Matrix generators
  , genMatrix
  , genVector
  , genSPDMatrix
  , genUpperTriangular
  , genLowerTriangular
  , genSymmetric
    -- * Tolerance
  , defaultTol
    -- * Standard test matrices
  , hilbertMatrix
  , wilkinsonMatrix
  , frankMatrix
  , hadamardMatrix
    -- * Generators with controlled properties
  , genMatrixWithCond
  , genNearSingularMatrix
  , genClusteredEigenMatrix
  , genSPDMatrixWithCond
  ) where

import qualified Data.Massiv.Array as M
import Data.Massiv.Array (Ix2(..), Sz(..), Comp(..))
import GHC.TypeNats (KnownNat, natVal)
import Data.Proxy (Proxy(..))
import Data.Bits (popCount)
import qualified Data.Bits as Bits
import Test.QuickCheck

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.BLAS.Level3 (matMul, transpose)
import Numeric.LinearAlgebra.Massiv.Orthogonal.QR (qr)

-- | Default tolerance for floating-point comparisons.
defaultTol :: Double
defaultTol = 1e-8

-- | Approximate equality for scalars.
(~=) :: Double -> Double -> Bool
x ~= y = abs (x - y) < defaultTol * (1 + abs x + abs y)

-- | Approximate equality for matrices.
matApproxEq :: forall m n. (KnownNat m, KnownNat n)
            => Matrix m n M.P Double -> Matrix m n M.P Double -> Bool
matApproxEq a b =
  let r = dimVal @m
      c = dimVal @n
  in all (\(i, j) -> (a ! (i, j)) ~= (b ! (i, j)))
     [(i, j) | i <- [0..r-1], j <- [0..c-1]]

-- | Approximate equality for vectors.
vecApproxEq :: forall n. KnownNat n
            => Vector n M.P Double -> Vector n M.P Double -> Bool
vecApproxEq a b =
  let nn = dimVal @n
  in all (\i -> (a !. i) ~= (b !. i)) [0..nn-1]

-- | Generate a random m×n matrix with entries in [-10, 10].
genMatrix :: forall m n. (KnownNat m, KnownNat n) => Gen (Matrix m n M.P Double)
genMatrix = do
  let r = fromIntegral (natVal (Proxy @m))
      c = fromIntegral (natVal (Proxy @n))
  entries <- vectorOf (r * c) (choose (-10, 10))
  pure $ makeMatrix @m @n @M.P $ \i j -> entries !! (i * c + j)

-- | Generate a random n-element vector with entries in [-10, 10].
genVector :: forall n. KnownNat n => Gen (Vector n M.P Double)
genVector = do
  let nn = fromIntegral (natVal (Proxy @n))
  entries <- vectorOf nn (choose (-10, 10))
  pure $ makeVector @n @M.P $ \i -> entries !! i

-- | Generate a symmetric positive definite matrix: A = BBᵀ + εI.
genSPDMatrix :: forall n. KnownNat n => Gen (Matrix n n M.P Double)
genSPDMatrix = do
  let nn = fromIntegral (natVal (Proxy @n))
  entries <- vectorOf (nn * nn) (choose (-5, 5))
  let b = makeMatrix @n @n @M.P $ \i j -> entries !! (i * nn + j)
      -- BBᵀ + εI
      epsilon = 0.1 :: Double
  pure $ makeMatrix @n @n @M.P $ \i j ->
    let bbT = foldl' (\acc k -> acc + (b ! (i, k)) * (b ! (j, k))) 0 [0..nn-1]
    in bbT + if i == j then epsilon else 0

-- | Generate an upper triangular matrix with nonzero diagonal.
genUpperTriangular :: forall n. KnownNat n => Gen (Matrix n n M.P Double)
genUpperTriangular = do
  let nn = fromIntegral (natVal (Proxy @n))
  entries <- vectorOf (nn * nn) (choose (-10, 10))
  diags <- vectorOf nn (choose (1, 10))  -- Ensure nonzero diagonal
  pure $ makeMatrix @n @n @M.P $ \i j ->
    if i == j then diags !! i
    else if i < j then entries !! (i * nn + j)
    else 0

-- | Generate a lower triangular matrix with nonzero diagonal.
genLowerTriangular :: forall n. KnownNat n => Gen (Matrix n n M.P Double)
genLowerTriangular = do
  let nn = fromIntegral (natVal (Proxy @n))
  entries <- vectorOf (nn * nn) (choose (-10, 10))
  diags <- vectorOf nn (choose (1, 10))
  pure $ makeMatrix @n @n @M.P $ \i j ->
    if i == j then diags !! i
    else if i > j then entries !! (i * nn + j)
    else 0

-- | Generate a symmetric matrix.
genSymmetric :: forall n. KnownNat n => Gen (Matrix n n M.P Double)
genSymmetric = do
  let nn = fromIntegral (natVal (Proxy @n))
  entries <- vectorOf (nn * nn) (choose (-10, 10))
  pure $ makeMatrix @n @n @M.P $ \i j ->
    if i <= j then entries !! (i * nn + j)
    else entries !! (j * nn + i)

------------------------------------------------------------------------
-- Standard test matrices
------------------------------------------------------------------------

-- | Hilbert matrix of order /n/: \(H_{ij} = 1/(i+j+1)\).
--
-- The Hilbert matrix is symmetric positive definite with an
-- exponentially growing condition number. Its inverse has integer
-- entries.
--
-- /Reference:/ Higham (2002), Section 28.1; GVL4 p. 128.
hilbertMatrix :: forall n. KnownNat n => Matrix n n M.P Double
hilbertMatrix = makeMatrix @n @n @M.P $ \i j ->
  1 / fromIntegral (i + j + 1)

-- | Wilkinson matrix \(W_n\): symmetric tridiagonal with
-- diagonal \(|i - \lfloor n/2 \rfloor|\) and unit sub/superdiagonal.
--
-- Has near-degenerate eigenvalue pairs that stress eigenvalue solvers.
--
-- /Reference:/ Higham (2002), Section 28.6; MATLAB @gallery('wilk', n)@.
wilkinsonMatrix :: forall n. KnownNat n => Matrix n n M.P Double
wilkinsonMatrix =
  let nn = dimVal @n
      mid = nn `div` 2
  in makeMatrix @n @n @M.P $ \i j ->
    if i == j then fromIntegral (abs (i - mid))
    else if abs (i - j) == 1 then 1
    else 0

-- | Frank matrix of order /n/: upper Hessenberg with
-- \(\det(F) = 1\) and known positive real eigenvalues.
--
-- \(F_{ij} = n - \max(i,j)\) for \(j \ge i-1\), zero otherwise
-- (0-indexed).
--
-- /Reference:/ Higham (2002), Section 28.5; Frank (1958).
frankMatrix :: forall n. KnownNat n => Matrix n n M.P Double
frankMatrix =
  let nn = dimVal @n
  in makeMatrix @n @n @M.P $ \i j ->
    if j >= i - 1 && i - 1 >= 0 || j >= i
    then fromIntegral (nn - max i j)
    else 0

-- | Hadamard matrix of order /n/ via the Sylvester\/Walsh construction.
--
-- \(H_{ij} = (-1)^{\mathrm{popcount}(i \mathbin{\&} j)}\)
--
-- Produces a proper Hadamard matrix when /n/ is a power of 2,
-- satisfying \(H^T H = n I\). All entries are \(\pm 1\).
--
-- /Reference:/ Higham (2002), Section 28.3.
hadamardMatrix :: forall n. KnownNat n => Matrix n n M.P Double
hadamardMatrix = makeMatrix @n @n @M.P $ \i j ->
  if even (popCount ((Bits..&.) i j)) then 1 else -1

------------------------------------------------------------------------
-- Generators with controlled properties
------------------------------------------------------------------------

-- | Generate a matrix with prescribed 2-norm condition number /κ/.
--
-- Uses the @randsvd@ construction: \(A = U \Sigma V^T\) where
-- \(U\) and \(V\) are random orthogonal matrices obtained from
-- QR factorization of random matrices, and the singular values
-- are geometrically spaced from 1 to \(1/\kappa\).
--
-- /Reference:/ Higham (2002), Section 28.3; Fasi & Higham,
-- "Generating Extreme-Scale Matrices With Specified Singular
-- Values or Condition Number," SIAM J. Sci. Comput. 43(5), 2021.
genMatrixWithCond :: forall m n. (KnownNat m, KnownNat n)
  => Double -> Gen (Matrix m n M.P Double)
genMatrixWithCond kappa = do
  uRaw <- genMatrix @m @m
  vRaw <- genMatrix @n @n
  let (u, _) = qr uRaw
      (v, _) = qr vRaw
      mm = dimVal @m
      nn = dimVal @n
      minDim = min mm nn
      -- Singular values geometrically spaced: sigma_0 = 1, sigma_{minDim-1} = 1/kappa
      sigma = makeMatrix @m @n @M.P $ \i j ->
        if i == j && i < minDim
        then if minDim <= 1 then 1
             else let t = fromIntegral i / fromIntegral (minDim - 1)
                  in kappa ** (-t)
        else 0
  pure $ matMul u (matMul sigma (transpose v))

-- | Generate a near-singular matrix with smallest singular value /ε/.
--
-- Wrapper around 'genMatrixWithCond' with \(\kappa = 1/\varepsilon\).
genNearSingularMatrix :: forall n. KnownNat n
  => Double -> Gen (Matrix n n M.P Double)
genNearSingularMatrix smallSV = genMatrixWithCond @n @n (1 / smallSV)

-- | Generate a symmetric matrix with clustered eigenvalues near /c/.
--
-- Constructs \(A = Q \Lambda Q^T\) where \(\Lambda\) is diagonal
-- with entries \(c + i \cdot 10^{-6}\) and \(Q\) is a random
-- orthogonal matrix. This stresses iterative eigenvalue solvers
-- that must resolve near-degenerate eigenvalue pairs.
genClusteredEigenMatrix :: forall n. KnownNat n
  => Double -> Gen (Matrix n n M.P Double)
genClusteredEigenMatrix cluster = do
  qRaw <- genMatrix @n @n
  let (q, _) = qr qRaw
      qt = transpose q
      lambda = makeMatrix @n @n @M.P $ \i j ->
        if i == j then cluster + fromIntegral i * 1e-6 else 0
  pure $ matMul q (matMul lambda qt)

-- | Generate an SPD matrix with prescribed condition number /κ/.
--
-- Uses the @randsvd@ construction with \(U = V\) (ensuring symmetry):
-- \(A = Q \Lambda Q^T\) where eigenvalues are geometrically spaced
-- from 1 to \(1/\kappa\).
genSPDMatrixWithCond :: forall n. KnownNat n
  => Double -> Gen (Matrix n n M.P Double)
genSPDMatrixWithCond kappa = do
  qRaw <- genMatrix @n @n
  let (q, _) = qr qRaw
      qt = transpose q
      nn = dimVal @n
      lambda = makeMatrix @n @n @M.P $ \i j ->
        if i == j
        then if nn <= 1 then 1
             else let t = fromIntegral i / fromIntegral (nn - 1)
                  in kappa ** (-t)
        else 0
  pure $ matMul q (matMul lambda qt)
