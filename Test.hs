module Main where

import LevMar (dlevmar_dif, defaultOpts)

ys :: Fractional a => [a]
ys = [ 2.90, 3.40, 3.10, 2.79, 3.15
     , 3.37, 2.87, 3.67, 2.57, 2.78
     , 2.96, 3.03, 3.34, 2.87, 2.20
     ]

params :: Num a => [a]
params = [0]

model :: [a] -> Int -> a
model [p] _ = p
model _   _ = error "FOUT!"

main :: IO ()
main = do result <- dlevmar_dif model params (zip [0..] ys) 1000 defaultOpts
          print result
