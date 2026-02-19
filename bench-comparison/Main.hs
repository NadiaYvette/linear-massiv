{-# LANGUAGE AllowAmbiguousTypes #-}

-- | Cross-library benchmark: linear-massiv vs hmatrix vs linear.
--
-- Compares performance of numerical linear algebra operations across
-- three Haskell libraries:
--
--   * __linear-massiv__ — pure Haskell, massiv-backed, type-safe dimensions
--   * __hmatrix__       — FFI to BLAS\/LAPACK (OpenBLAS on this system)
--   * __linear__        — pure Haskell, optimised for small fixed-size (V2–V4)
--
-- Run with @+RTS -N1@ for fair single-threaded comparison.
module Main (main) where

import Criterion.Main
import qualified Data.Massiv.Array as M
import GHC.TypeNats (KnownNat)

-- linear-massiv
import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.BLAS.Level1 (dotP)
import Numeric.LinearAlgebra.Massiv.BLAS.Level2 (matvecP)
import Numeric.LinearAlgebra.Massiv.BLAS.Level3 (matMulP)
import Numeric.LinearAlgebra.Massiv.Solve.LU (luSolve, luSolveP)
import Numeric.LinearAlgebra.Massiv.Solve.Cholesky (choleskySolve, choleskySolveP)
import Numeric.LinearAlgebra.Massiv.Orthogonal.QR (qr, qrP)
import Numeric.LinearAlgebra.Massiv.Eigen.Symmetric (symmetricEigen, symmetricEigenP)
import Numeric.LinearAlgebra.Massiv.Eigen.SVD (svd)

-- hmatrix
import qualified Numeric.LinearAlgebra as H

-- linear (small fixed-size)
import Linear.V4 (V4(..))
import qualified Linear.Matrix as LM
import qualified Linear.Metric as LMet

------------------------------------------------------------------------
-- linear-massiv matrix generators
------------------------------------------------------------------------

mkMatLM :: forall m n. (KnownNat m, KnownNat n) => Matrix m n M.P Double
mkMatLM = makeMatrix @m @n @M.P $ \i j ->
  fromIntegral (i * 7 + j * 3 + 1) / 100.0

mkVecLM :: forall n. KnownNat n => Vector n M.P Double
mkVecLM = makeVector @n @M.P $ \i -> fromIntegral (i + 1) / 10.0

-- Diagonally dominant for LU
mkDDLM :: forall n. KnownNat n => Matrix n n M.P Double
mkDDLM = makeMatrix @n @n @M.P $ \i j ->
  fromIntegral (i * 7 + j * 3 + 1) / 100.0
    + if i == j then fromIntegral (dimVal @n) else 0

-- SPD: B^T B + nI
mkSPDLM :: forall n. KnownNat n => Matrix n n M.P Double
mkSPDLM =
  let nn = dimVal @n
      b = makeMatrix @n @n @M.P $ \i j ->
            fromIntegral (i * nn + j + 1) / fromIntegral (nn * nn)
  in makeMatrix @n @n @M.P $ \i j ->
       foldl' (\acc k -> acc + (b ! (i, k)) * (b ! (j, k)))
              (if i == j then 1 else 0)
              [0..nn-1]

------------------------------------------------------------------------
-- hmatrix matrix generators (same numerical entries)
------------------------------------------------------------------------

mkMatHM :: Int -> H.Matrix Double
mkMatHM n = (n H.>< n) [ fromIntegral (i * 7 + j * 3 + 1) / 100.0
                        | i <- [0..n-1], j <- [0..n-1] ]

mkVecHM :: Int -> H.Vector Double
mkVecHM n = H.fromList [ fromIntegral (i + 1) / 10.0 | i <- [0..n-1] ]

mkDDHM :: Int -> H.Matrix Double
mkDDHM n = (n H.>< n) [ fromIntegral (i * 7 + j * 3 + 1) / 100.0
                           + if i == j then fromIntegral n else 0
                       | i <- [0..n-1], j <- [0..n-1] ]

mkSPDHM :: Int -> H.Matrix Double
mkSPDHM n =
  let b = (n H.>< n) [ fromIntegral (i * n + j + 1) / fromIntegral (n * n)
                      | i <- [0..n-1], j <- [0..n-1] ]
  in H.tr b H.<> b + H.scale (fromIntegral n) (H.ident n)

------------------------------------------------------------------------
-- linear (V4) data
------------------------------------------------------------------------

linM44a :: V4 (V4 Double)
linM44a = V4 (V4 0.01 0.04 0.07 0.10)
              (V4 0.08 0.11 0.14 0.17)
              (V4 0.15 0.18 0.21 0.24)
              (V4 0.22 0.25 0.28 0.31)

linM44b :: V4 (V4 Double)
linM44b = V4 (V4 0.34 0.37 0.40 0.43)
              (V4 0.41 0.44 0.47 0.50)
              (V4 0.48 0.51 0.54 0.57)
              (V4 0.55 0.58 0.61 0.64)

linV4a :: V4 Double
linV4a = V4 0.1 0.2 0.3 0.4

linV4b :: V4 Double
linV4b = V4 0.5 0.6 0.7 0.8

------------------------------------------------------------------------
-- hmatrix helpers (avoid operator section issues)
------------------------------------------------------------------------

hmGemm :: H.Matrix Double -> H.Matrix Double -> H.Matrix Double
hmGemm = (H.<>)

hmMatvec :: H.Matrix Double -> H.Vector Double -> H.Vector Double
hmMatvec = (H.#>)

hmDot :: H.Vector Double -> H.Vector Double -> Double
hmDot = H.dot

hmLinearSolve :: H.Matrix Double -> H.Vector Double -> H.Matrix Double
hmLinearSolve a b = case H.linearSolve a (H.asColumn b) of
  Just x  -> x
  Nothing -> error "hmLinearSolve: singular matrix"

hmCholSolve :: H.Matrix Double -> H.Vector Double -> H.Matrix Double
hmCholSolve a b =
  let r = H.chol (H.trustSym a)
  in H.cholSolve r (H.asColumn b)

------------------------------------------------------------------------
-- Benchmarks
------------------------------------------------------------------------

main :: IO ()
main = do
  -- Pre-compute all matrices to avoid construction overhead in benchmarks
  let hm4   = mkMatHM 4;   hm10  = mkMatHM 10;  hm50  = mkMatHM 50
      hm100 = mkMatHM 100; hm200 = mkMatHM 200
      hv4   = mkVecHM 4;   hv10  = mkVecHM 10;   hv50  = mkVecHM 50
      hv100 = mkVecHM 100; hv1000 = mkVecHM 1000
      dd10 = mkDDHM 10; dd50 = mkDDHM 50; dd100 = mkDDHM 100
      spd10 = mkSPDHM 10; spd50 = mkSPDHM 50; spd100 = mkSPDHM 100

  defaultMain
    [ bgroup "GEMM"
      [ bgroup "4x4"
        [ bench "linear"        $ nf (linM44a LM.!*!) linM44b
        , bench "hmatrix"       $ nf (hmGemm hm4) hm4
        , bench "linear-massiv" $ nf (matMulP (mkMatLM @4 @4)) (mkMatLM @4 @4)
        ]
      , bgroup "10x10"
        [ bench "hmatrix"       $ nf (hmGemm hm10) hm10
        , bench "linear-massiv" $ nf (matMulP (mkMatLM @10 @10)) (mkMatLM @10 @10)
        ]
      , bgroup "50x50"
        [ bench "hmatrix"       $ nf (hmGemm hm50) hm50
        , bench "linear-massiv" $ nf (matMulP (mkMatLM @50 @50)) (mkMatLM @50 @50)
        ]
      , bgroup "100x100"
        [ bench "hmatrix"       $ nf (hmGemm hm100) hm100
        , bench "linear-massiv" $ nf (matMulP (mkMatLM @100 @100)) (mkMatLM @100 @100)
        ]
      , bgroup "200x200"
        [ bench "hmatrix"       $ nf (hmGemm hm200) hm200
        , bench "linear-massiv" $ nf (matMulP (mkMatLM @200 @200)) (mkMatLM @200 @200)
        ]
      ]
    , bgroup "dot"
      [ bgroup "4"
        [ bench "linear"        $ nf (LMet.dot linV4a) linV4b
        , bench "hmatrix"       $ nf (hmDot hv4) hv4
        , bench "linear-massiv" $ nf (dotP (mkVecLM @4)) (mkVecLM @4)
        ]
      , bgroup "100"
        [ bench "hmatrix"       $ nf (hmDot hv100) hv100
        , bench "linear-massiv" $ nf (dotP (mkVecLM @100)) (mkVecLM @100)
        ]
      , bgroup "1000"
        [ bench "hmatrix"       $ nf (hmDot hv1000) hv1000
        , bench "linear-massiv" $ nf (dotP (mkVecLM @1000)) (mkVecLM @1000)
        ]
      ]
    , bgroup "matvec"
      [ bgroup "4"
        [ bench "linear"        $ nf (linM44a LM.!*) linV4a
        , bench "hmatrix"       $ nf (hmMatvec hm4) hv4
        , bench "linear-massiv" $ nf (matvecP (mkMatLM @4 @4)) (mkVecLM @4)
        ]
      , bgroup "50"
        [ bench "hmatrix"       $ nf (hmMatvec hm50) hv50
        , bench "linear-massiv" $ nf (matvecP (mkMatLM @50 @50)) (mkVecLM @50)
        ]
      , bgroup "100"
        [ bench "hmatrix"       $ nf (hmMatvec hm100) hv100
        , bench "linear-massiv" $ nf (matvecP (mkMatLM @100 @100)) (mkVecLM @100)
        ]
      ]
    , bgroup "luSolve"
      [ bgroup "10x10"
        [ bench "hmatrix"       $ nf (hmLinearSolve dd10) hv10
        , bench "linear-massiv" $ nf (luSolveP (mkDDLM @10)) (mkVecLM @10)
        , bench "lm-generic"    $ nf (luSolve (mkDDLM @10)) (mkVecLM @10)
        ]
      , bgroup "50x50"
        [ bench "hmatrix"       $ nf (hmLinearSolve dd50) hv50
        , bench "linear-massiv" $ nf (luSolveP (mkDDLM @50)) (mkVecLM @50)
        , bench "lm-generic"    $ nf (luSolve (mkDDLM @50)) (mkVecLM @50)
        ]
      , bgroup "100x100"
        [ bench "hmatrix"       $ nf (hmLinearSolve dd100) hv100
        , bench "linear-massiv" $ nf (luSolveP (mkDDLM @100)) (mkVecLM @100)
        , bench "lm-generic"    $ nf (luSolve (mkDDLM @100)) (mkVecLM @100)
        ]
      ]
    , bgroup "choleskySolve"
      [ bgroup "10x10"
        [ bench "hmatrix"       $ nf (hmCholSolve spd10) hv10
        , bench "linear-massiv" $ nf (choleskySolveP (mkSPDLM @10)) (mkVecLM @10)
        , bench "lm-generic"    $ nf (choleskySolve (mkSPDLM @10)) (mkVecLM @10)
        ]
      , bgroup "50x50"
        [ bench "hmatrix"       $ nf (hmCholSolve spd50) hv50
        , bench "linear-massiv" $ nf (choleskySolveP (mkSPDLM @50)) (mkVecLM @50)
        , bench "lm-generic"    $ nf (choleskySolve (mkSPDLM @50)) (mkVecLM @50)
        ]
      , bgroup "100x100"
        [ bench "hmatrix"       $ nf (hmCholSolve spd100) hv100
        , bench "linear-massiv" $ nf (choleskySolveP (mkSPDLM @100)) (mkVecLM @100)
        , bench "lm-generic"    $ nf (choleskySolve (mkSPDLM @100)) (mkVecLM @100)
        ]
      ]
    , bgroup "QR"
      [ bgroup "10x10"
        [ bench "hmatrix"       $ nf H.qr hm10
        , bench "linear-massiv" $ nf qrP (mkMatLM @10 @10)
        , bench "lm-generic"    $ nf qr (mkMatLM @10 @10)
        ]
      , bgroup "50x50"
        [ bench "hmatrix"       $ nf H.qr hm50
        , bench "linear-massiv" $ nf qrP (mkMatLM @50 @50)
        , bench "lm-generic"    $ nf qr (mkMatLM @50 @50)
        ]
      , bgroup "100x100"
        [ bench "hmatrix"       $ nf H.qr hm100
        , bench "linear-massiv" $ nf qrP (mkMatLM @100 @100)
        , bench "lm-generic"    $ nf qr (mkMatLM @100 @100)
        ]
      ]
    , bgroup "eigenSH"
      [ bgroup "10x10"
        [ bench "hmatrix"       $ nf (H.eigSH . H.trustSym) spd10
        , bench "linear-massiv" $ nf (\a -> symmetricEigenP a 200 1e-12) (mkSPDLM @10)
        , bench "lm-generic"    $ nf (\a -> symmetricEigen a 200 1e-12) (mkSPDLM @10)
        ]
      , bgroup "50x50"
        [ bench "hmatrix"       $ nf (H.eigSH . H.trustSym) spd50
        , bench "linear-massiv" $ nf (\a -> symmetricEigenP a 500 1e-12) (mkSPDLM @50)
        , bench "lm-generic"    $ nf (\a -> symmetricEigen a 500 1e-12) (mkSPDLM @50)
        ]
      ]
    , bgroup "SVD"
      [ bgroup "10x10"
        [ bench "hmatrix"       $ nf H.svd hm10
        , bench "linear-massiv" $ nf svd (mkMatLM @10 @10)
        ]
      , bgroup "50x50"
        [ bench "hmatrix"       $ nf H.svd hm50
        , bench "linear-massiv" $ nf svd (mkMatLM @50 @50)
        ]
      ]
    ]
