{-# LANGUAGE AllowAmbiguousTypes #-}

module Bench.Orthogonal (orthogonalBenchmarks) where

import Criterion.Main
import qualified Data.Massiv.Array as M
import GHC.TypeNats (KnownNat)

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.Orthogonal.QR (qr, qrGivens)

mkMatP :: forall m n. (KnownNat m, KnownNat n) => Matrix m n M.P Double
mkMatP = makeMatrix @m @n @M.P $ \i j -> fromIntegral (i * 7 + j * 3 + 1) / 100.0

orthogonalBenchmarks :: [Benchmark]
orthogonalBenchmarks =
  [ bgroup "qr-householder"
    [ bench "10x10/P"  $ nf qr (mkMatP @10 @10)
    , bench "50x50/P"  $ nf qr (mkMatP @50 @50)
    , bench "100x100/P" $ nf qr (mkMatP @100 @100)
    ]
  , bgroup "qr-givens"
    [ bench "10x10/P"  $ nf qrGivens (mkMatP @10 @10)
    , bench "50x50/P"  $ nf qrGivens (mkMatP @50 @50)
    ]
  ]
