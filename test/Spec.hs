module Main (main) where

import Test.Tasty

import Test.BLAS (blasTests)
import Test.Solve (solveTests)
import Test.Orthogonal (orthogonalTests)
import Test.Eigen (eigenTests)
import Test.Norms (normTests)

main :: IO ()
main = defaultMain $ testGroup "linear-massiv"
  [ blasTests
  , solveTests
  , orthogonalTests
  , eigenTests
  , normTests
  ]
