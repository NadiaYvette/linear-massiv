{-# LANGUAGE AllowAmbiguousTypes #-}

module Bench.Solve (solveBenchmarks) where

import Criterion.Main
import qualified Data.Massiv.Array as M
import GHC.TypeNats (KnownNat)

import Numeric.LinearAlgebra.Massiv.Types
import Numeric.LinearAlgebra.Massiv.Internal
import Numeric.LinearAlgebra.Massiv.Solve.LU (lu, luSolve)
import Numeric.LinearAlgebra.Massiv.Solve.Cholesky (cholesky, choleskySolve)

-- Diagonally dominant matrix for LU
mkMatP :: forall n. KnownNat n => Matrix n n M.P Double
mkMatP = makeMatrix @n @n @M.P $ \i j ->
  fromIntegral (i * 7 + j * 3 + 1) / 100.0 + if i == j then fromIntegral (dimVal @n) else 0

-- SPD matrix: A = B^T B + nI
mkSPDP :: forall n. KnownNat n => Matrix n n M.P Double
mkSPDP =
  let nn = dimVal @n
      b = makeMatrix @n @n @M.P $ \i j -> fromIntegral (i * nn + j + 1) / fromIntegral (nn * nn)
  in makeMatrix @n @n @M.P $ \i j ->
    foldl' (\acc k -> acc + (b ! (i, k)) * (b ! (j, k))) (if i == j then 1 else 0) [0..nn-1]

mkVecP :: forall n. KnownNat n => Vector n M.P Double
mkVecP = makeVector @n @M.P $ \i -> fromIntegral (i + 1)

solveBenchmarks :: [Benchmark]
solveBenchmarks =
  [ bgroup "lu"
    [ bench "10x10/P"  $ nf lu (mkMatP @10)
    , bench "50x50/P"  $ nf lu (mkMatP @50)
    , bench "100x100/P" $ nf lu (mkMatP @100)
    ]
  , bgroup "luSolve"
    [ bench "10x10/P"  $ nf (uncurry luSolve) (mkMatP @10, mkVecP @10)
    , bench "50x50/P"  $ nf (uncurry luSolve) (mkMatP @50, mkVecP @50)
    , bench "100x100/P" $ nf (uncurry luSolve) (mkMatP @100, mkVecP @100)
    ]
  , bgroup "cholesky"
    [ bench "10x10/P"  $ nf cholesky (mkSPDP @10)
    , bench "50x50/P"  $ nf cholesky (mkSPDP @50)
    , bench "100x100/P" $ nf cholesky (mkSPDP @100)
    ]
  , bgroup "choleskySolve"
    [ bench "10x10/P"  $ nf (uncurry choleskySolve) (mkSPDP @10, mkVecP @10)
    , bench "50x50/P"  $ nf (uncurry choleskySolve) (mkSPDP @50, mkVecP @50)
    , bench "100x100/P" $ nf (uncurry choleskySolve) (mkSPDP @100, mkVecP @100)
    ]
  ]
