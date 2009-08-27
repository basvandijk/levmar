module Main where

import LevMar.Intermediate
import System.Random

model_1 :: Model a r
model_1 [a] _ = a
model_1 _ _ = error "params"

model_2 :: Num r => Model r r
model_2 [a, b] x = a*x + b
model_2 _ _ = error "params"

model_3 :: Num r => Model r r
model_3 [a, b, c] x = a*x*x + b*x + c
model_3 _ _ = error "params"

model_4 :: Num r => Model r r
model_4 [a, b, c, d] x = a*x*x*x + b*x*x + c*x + d
model_4 _ _ = error "params"

expfunc :: Floating r => [r] -> r -> r
expfunc [a, b, c] x = a * exp ((-b) * x) + c
expfunc _ _ = error "params"

rndGenSeed :: Int
rndGenSeed = 123456

test :: Show a
     => Model Double a
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
      samples = zip xs $ map (\x -> f ps x) xs
      samples' = zipWith (\(x, y) n -> (x, y + (n - 0.5) * 2 * noise)) samples ns

-- main :: IO ()
-- main = print $ test model_3 [5, 3, 7] [1..50] 2

main :: IO ()
main = print $ test_expfunc

test_expfunc = levmar expfunc
                      Nothing
                      ([1, 0, 0] :: [Double])
                      samples'
                      1000
                      ( defaultOpts { opt_epsilon1 = 1e-15
                                    , opt_epsilon2 = 1e-15
                                    , opt_epsilon3 = 1e-20
                                    }
                      )
                      Nothing
                      Nothing
    where
      ns = take (length xs) $ randoms $ mkStdGen rndGenSeed
      xs = [1..40]
      samples = zip xs $ map (\x -> expfunc [5.0, 0.1, 1.0] x) xs
      samples' = zipWith (\(x, y) n -> (x, y + (n - 0.5) * 2 * 0.1)) samples ns
