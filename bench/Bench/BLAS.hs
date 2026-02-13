{-# LANGUAGE AllowAmbiguousTypes #-}

module Bench.BLAS (blasBenchmarks) where

import Criterion.Main
import qualified Data.Massiv.Array as M
import Data.Massiv.Array (Comp(..))
import Control.DeepSeq ()
import GHC.TypeNats (KnownNat)

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.BLAS.Level1 (dot)
import Numeric.LinearAlgebra.Massiv.BLAS.Level2 (matvec)
import Numeric.LinearAlgebra.Massiv.BLAS.Level3 (matMul, matMulComp)

-- linear library imports for comparison
import Linear.V4 (V4(..))
import Linear.Matrix ((!*!), (!*))
import qualified Linear.Metric as LM
import Linear.V4 ()

-- Massiv matrix generators
mkMatP :: forall m n. (KnownNat m, KnownNat n) => Matrix m n M.P Double
mkMatP = makeMatrix @m @n @M.P $ \i j -> fromIntegral (i * 7 + j * 3 + 1) / 100.0

mkVecP :: forall n. KnownNat n => Vector n M.P Double
mkVecP = makeVector @n @M.P $ \i -> fromIntegral (i + 1) / 10.0

-- linear library 4x4 matrices for comparison
linearM44 :: V4 (V4 Double)
linearM44 = V4 (V4 1 2 3 4)
                (V4 5 6 7 8)
                (V4 9 10 11 12)
                (V4 13 14 15 16)

linearM44b :: V4 (V4 Double)
linearM44b = V4 (V4 17 18 19 20)
                 (V4 21 22 23 24)
                 (V4 25 26 27 28)
                 (V4 29 30 31 32)

linearV4 :: V4 Double
linearV4 = V4 1 2 3 4

linearV4b :: V4 Double
linearV4b = V4 5 6 7 8

blasBenchmarks :: [Benchmark]
blasBenchmarks =
  [ bgroup "gemm"
    [ -- Small sizes: compare linear vs massiv
      bgroup "4x4"
        [ bench "linear-V4" $ nf (linearM44 !*!) linearM44b
        , bench "massiv-P"  $ nf (uncurry (matMul @4 @4 @4)) (mkMatP @4 @4, mkMatP @4 @4)
        ]
    , bgroup "10x10"
        [ bench "P/Seq" $ nf (uncurry (matMul @10 @10 @10)) (mkMatP @10 @10, mkMatP @10 @10)
        ]
    , bgroup "50x50"
        [ bench "P/Seq" $ nf (uncurry (matMul @50 @50 @50)) (mkMatP @50 @50, mkMatP @50 @50)
        , bench "P/Par" $ nf (uncurry (matMulComp @50 @50 @50 Par)) (mkMatP @50 @50, mkMatP @50 @50)
        ]
    , bgroup "100x100"
        [ bench "P/Seq" $ nf (uncurry (matMul @100 @100 @100)) (mkMatP @100 @100, mkMatP @100 @100)
        , bench "P/Par" $ nf (uncurry (matMulComp @100 @100 @100 Par)) (mkMatP @100 @100, mkMatP @100 @100)
        ]
    , bgroup "200x200"
        [ bench "P/Seq" $ nf (uncurry (matMul @200 @200 @200)) (mkMatP @200 @200, mkMatP @200 @200)
        , bench "P/Par" $ nf (uncurry (matMulComp @200 @200 @200 Par)) (mkMatP @200 @200, mkMatP @200 @200)
        ]
    -- Representation comparison at 100x100
    , bgroup "repr-100x100"
        [ bench "P" $ nf (uncurry (matMul @100 @100 @100)) (mkMatP @100 @100, mkMatP @100 @100)
        ]
    ]
  , bgroup "dot"
    [ bgroup "4"
        [ bench "linear-V4" $ nf (LM.dot linearV4) linearV4b
        , bench "massiv-P"  $ nf (uncurry dot) (mkVecP @4, mkVecP @4)
        ]
    , bgroup "100"
        [ bench "P" $ nf (uncurry dot) (mkVecP @100, mkVecP @100)
        ]
    , bgroup "1000"
        [ bench "P" $ nf (uncurry dot) (mkVecP @1000, mkVecP @1000)
        ]
    , bgroup "10000"
        [ bench "P" $ nf (uncurry dot) (mkVecP @10000, mkVecP @10000)
        ]
    ]
  , bgroup "matvec"
    [ bgroup "4"
        [ bench "linear-V4" $ nf (linearM44 !*) linearV4
        , bench "massiv-P"  $ nf (uncurry (matvec @4 @4)) (mkMatP @4 @4, mkVecP @4)
        ]
    , bgroup "50"
        [ bench "P" $ nf (uncurry (matvec @50 @50)) (mkMatP @50 @50, mkVecP @50)
        ]
    , bgroup "100"
        [ bench "P" $ nf (uncurry (matvec @100 @100)) (mkMatP @100 @100, mkVecP @100)
        ]
    ]
  ]
