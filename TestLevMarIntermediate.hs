module Main where

import LevMar.Intermediate
import System.Random

model_1 :: Model a r
model_1 _ [a] = a
model_1 _ _ = error "params"

model_2 :: Num r => Model r r
model_2 x [a, b] = a*x + b
model_2 _ _ = error "params"

model_3 :: Num r => Model r r
model_3 x [a, b, c] = a*x*x + b*x + c
model_3 _ _ = error "params"

model_4 :: Num r => Model r r
model_4 x [a, b, c, d] = a*x*x*x + b*x*x + c*x + d
model_4 _ _ = error "params"

rndGenSeed :: Int
rndGenSeed = 123456

test :: Show a
     => Model a Double
     -> [Double] -- solution
     -> [a]      -- x-values of samples
     -> Double   -- noise level
     -> ([Double], Info Double, CovarMatrix Double)
test f ps xs noise = levmar f
                            Nothing
                            (replicate (length ps) 0) -- all params 0
                            samples'
                            1000
                            defaultOpts
                            Nothing
                            Nothing
    where
      ns = take (length xs) $ randoms $ mkStdGen rndGenSeed
      samples = zip xs $ map (\x -> f x ps) xs
      samples' = zipWith (\(x, y) n -> (x, y + (n - 0.5) * 2 * noise)) samples ns

main :: IO ()
main = print $ test model_3 [5, 3, 7] [1..50] 2
