{-# LANGUAGE AllowAmbiguousTypes #-}

module Bench.Eigen (eigenBenchmarks) where

import Criterion.Main
import qualified Data.Massiv.Array as M
import GHC.TypeNats (KnownNat)

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.Eigen.Symmetric (symmetricEigen, jacobiEigen)
import Numeric.LinearAlgebra.Massiv.Eigen.SVD (svd)

mkSPDP :: forall n. KnownNat n => Matrix n n M.P Double
mkSPDP =
  let nn = dimVal @n
      b = makeMatrix @n @n @M.P $ \i j -> fromIntegral (i * nn + j + 1) / fromIntegral (nn * nn)
  in makeMatrix @n @n @M.P $ \i j ->
    foldl' (\acc k -> acc + (b ! (i, k)) * (b ! (j, k))) (if i == j then 1 else 0) [0..nn-1]

mkMatP :: forall m n. (KnownNat m, KnownNat n) => Matrix m n M.P Double
mkMatP = makeMatrix @m @n @M.P $ \i j -> fromIntegral (i * 7 + j * 3 + 1) / 100.0

eigenBenchmarks :: [Benchmark]
eigenBenchmarks =
  [ bgroup "symmetricEigen"
    [ bench "10x10/P" $ nf (\a -> symmetricEigen a 200 1e-12) (mkSPDP @10)
    , bench "20x20/P" $ nf (\a -> symmetricEigen a 500 1e-12) (mkSPDP @20)
    , bench "50x50/P" $ nf (\a -> symmetricEigen a 500 1e-12) (mkSPDP @50)
    ]
  , bgroup "jacobiEigen"
    [ bench "10x10/P" $ nf (\a -> jacobiEigen a 50 1e-12) (mkSPDP @10)
    , bench "20x20/P" $ nf (\a -> jacobiEigen a 50 1e-12) (mkSPDP @20)
    ]
  , bgroup "svd"
    [ bench "10x10/P" $ nf svd (mkMatP @10 @10)
    , bench "20x20/P" $ nf svd (mkMatP @20 @20)
    , bench "50x50/P" $ nf svd (mkMatP @50 @50)
    ]
  ]
