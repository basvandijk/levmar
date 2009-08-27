module Main where

import LevMar.Intermediate
import System.Random

model1 :: Model a r
model1 [a] _ = a
model1 _ _ = error "params"

model2 :: Num r => Model r r
model2 [a, b] x = a*x + b
model2 _ _ = error "params"

model3 :: Num r => Model r r
model3 [a, b, c] x = a*x*x + b*x + c
model3 _ _ = error "params"

model4 :: Num r => Model r r
model4 [a, b, c, d] x = a*x*x*x + b*x*x + c*x + d
model4 _ _ = error "params"

expFunc :: Floating r => [r] -> r -> r
expFunc [a, b, c] x = a * exp ((-b) * x) + c
expFunc _ _ = error "params"

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
                            Nothing
                            Nothing
    where
      ns = take (length xs) $ randoms $ mkStdGen rndGenSeed
      samples = zip xs $ map (f ps) xs
      samples' = zipWith (\(x, y) n -> (x, y + (n - 0.5) * 2 * noise)) samples ns

-- main :: IO ()
-- main = print $ test model3 [5, 3, 7] [1..50] 2

main :: IO ()
main = print testExpFunc

testExpFunc = levmar expFunc
                     Nothing
                     ([1, 0, 0] :: [Double])
                     samples'
                     1000
                     ( defaultOpts { optEpsilon1 = 1e-15
                                   , optEpsilon2 = 1e-15
                                   , optEpsilon3 = 1e-20
                                   }
                     )
                     Nothing
                     Nothing
                     Nothing
                     Nothing
    where
      ns = take (length xs) $ randoms $ mkStdGen rndGenSeed
      xs = [1..40]
      samples = zip xs $ map (expFunc [5.0, 0.1, 1.0]) xs
      samples' = zipWith (\(x, y) n -> (x, y + (n - 0.5) * 2 * 0.1)) samples ns
