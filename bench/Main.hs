module Main (main) where

import Criterion.Main

import Bench.BLAS (blasBenchmarks)
import Bench.Solve (solveBenchmarks)
import Bench.Orthogonal (orthogonalBenchmarks)
import Bench.Eigen (eigenBenchmarks)
import Bench.Parallel (parallelBenchmarks)

main :: IO ()
main = defaultMain
  [ bgroup "BLAS" blasBenchmarks
  , bgroup "Solve" solveBenchmarks
  , bgroup "Orthogonal" orthogonalBenchmarks
  , bgroup "Eigen" eigenBenchmarks
  , bgroup "Parallel" parallelBenchmarks
  ]
