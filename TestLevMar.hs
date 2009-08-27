module Main where

import LevMar
import SizedList (lengthSL, replicateSL)
import NFunction(($*))
import System.Random

model_1 :: r -> a -> r
model_1 a _ = a

model_2 :: Num r => r -> r -> r -> r
model_2 a b x = a*x + b

model_3 :: Num r => r -> r -> r -> r -> r
model_3 a b c x = a*x*x + b*x + c

model_4 :: Num r => r -> r -> r -> r -> r -> r
model_4 a b c d x = a*x*x*x + b*x*x + c*x + d

rndGenSeed :: Int
rndGenSeed = 123456

test :: (Show a, Nat n)
     => (Model n Double a)
     -> SizedList n Double
     -> [a]
     -> Double
     -> (SizedList n Double, Info Double, CovarMatrix n Double)
test f ps xs noise = levmar f
                            Nothing
                            (replicateSL (lengthSL ps) 0)
                            samples'
                            1000
                            defaultOpts
                            Nothing
                            Nothing
    where
      ns = take (length xs) $ randoms $ mkStdGen rndGenSeed
      samples = zip xs $ map (\x -> (f $* ps) x) xs
      samples' = zipWith (\(x, y) n -> (x, y + (n - 0.5) * 2 * noise)) samples ns

main :: IO ()
main = print $ test model_3 (5 ::: 3 ::: 7 ::: Nil) [1..50] 2

