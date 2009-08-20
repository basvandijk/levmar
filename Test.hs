module Main where

import LevMar
import System.Random

model_1 :: [r] -> a -> r
model_1 [a] _ = a
model_1 _ _ = error "params"

model_2 :: Num r => [r] -> r -> r
model_2 [a, b] x = a*x + b
model_2 _ _ = error "params"

model_3 :: Num r => [r] -> r -> r
model_3 [a, b, c] x = a*x*x + b*x + c
model_3 _ _ = error "params"

model_4 :: Num r => [r] -> r -> r
model_4 [a, b, c, d] x = a*x*x*x + b*x*x + c*x + d
model_4 _ _ = error "params"

rndGenSeed :: Int
rndGenSeed = 123456

test :: Show a
     => Model Double a
     -> [Double] -- solution
     -> [a]      -- x-values of samples
     -> Double   -- noise level
     -> IO ()
test f ps xs noise = do let ns = take (length xs) $ randoms $ mkStdGen rndGenSeed
                        let samples = zip xs $ map (f ps) xs
                        let samples' = zipWith (\(x, y) n -> (x, y + (n - 0.5) * 2 * noise)) samples ns
                        result <- dlevmar_dif f
                                              (replicate (length ps) 0) -- all params 0
                                              samples'
                                              1000
                                              defaultOpts
                        print result


main :: IO ()
main = test model_3 [5, 3, 7] [1..50] 2
