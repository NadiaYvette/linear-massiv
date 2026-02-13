{-# LANGUAGE AllowAmbiguousTypes #-}

-- | Benchmarks specifically for measuring parallelism scalability.
--
-- Runs matrix multiplication at various sizes with explicit thread counts
-- using massiv's ParN constructor.
module Bench.Parallel (parallelBenchmarks) where

import Criterion.Main
import qualified Data.Massiv.Array as M
import Data.Massiv.Array (Comp(..))
import GHC.TypeNats (KnownNat)

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.BLAS.Level3 (matMulComp)

mkMatP :: forall m n. (KnownNat m, KnownNat n) => Matrix m n M.P Double
mkMatP = makeMatrix @m @n @M.P $ \i j -> fromIntegral (i * 7 + j * 3 + 1) / 100.0

-- | Benchmarks measuring how performance scales with thread count.
parallelBenchmarks :: [Benchmark]
parallelBenchmarks =
  [ bgroup "scaling-gemm"
    [ bgroup "100x100"
        [ bench "Seq"     $ nf (uncurry (matMulComp @100 @100 @100 Seq))     (mkMatP @100 @100, mkMatP @100 @100)
        , bench "Par"     $ nf (uncurry (matMulComp @100 @100 @100 Par))     (mkMatP @100 @100, mkMatP @100 @100)
        , bench "ParN-1"  $ nf (uncurry (matMulComp @100 @100 @100 (ParN 1)))  (mkMatP @100 @100, mkMatP @100 @100)
        , bench "ParN-2"  $ nf (uncurry (matMulComp @100 @100 @100 (ParN 2)))  (mkMatP @100 @100, mkMatP @100 @100)
        , bench "ParN-4"  $ nf (uncurry (matMulComp @100 @100 @100 (ParN 4)))  (mkMatP @100 @100, mkMatP @100 @100)
        , bench "ParN-8"  $ nf (uncurry (matMulComp @100 @100 @100 (ParN 8)))  (mkMatP @100 @100, mkMatP @100 @100)
        , bench "ParN-16" $ nf (uncurry (matMulComp @100 @100 @100 (ParN 16))) (mkMatP @100 @100, mkMatP @100 @100)
        , bench "ParN-20" $ nf (uncurry (matMulComp @100 @100 @100 (ParN 20))) (mkMatP @100 @100, mkMatP @100 @100)
        ]
    , bgroup "200x200"
        [ bench "Seq"     $ nf (uncurry (matMulComp @200 @200 @200 Seq))     (mkMatP @200 @200, mkMatP @200 @200)
        , bench "Par"     $ nf (uncurry (matMulComp @200 @200 @200 Par))     (mkMatP @200 @200, mkMatP @200 @200)
        , bench "ParN-1"  $ nf (uncurry (matMulComp @200 @200 @200 (ParN 1)))  (mkMatP @200 @200, mkMatP @200 @200)
        , bench "ParN-2"  $ nf (uncurry (matMulComp @200 @200 @200 (ParN 2)))  (mkMatP @200 @200, mkMatP @200 @200)
        , bench "ParN-4"  $ nf (uncurry (matMulComp @200 @200 @200 (ParN 4)))  (mkMatP @200 @200, mkMatP @200 @200)
        , bench "ParN-8"  $ nf (uncurry (matMulComp @200 @200 @200 (ParN 8)))  (mkMatP @200 @200, mkMatP @200 @200)
        , bench "ParN-16" $ nf (uncurry (matMulComp @200 @200 @200 (ParN 16))) (mkMatP @200 @200, mkMatP @200 @200)
        , bench "ParN-20" $ nf (uncurry (matMulComp @200 @200 @200 (ParN 20))) (mkMatP @200 @200, mkMatP @200 @200)
        ]
    ]
  ]
