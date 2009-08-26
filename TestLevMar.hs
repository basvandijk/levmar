module Main where

import LevMar

import TypeLevelNat (N(..))
import SizedList (lengthSL, replicateSL)

import System.Random

model_1 :: a -> r -> r
model_1 _ a = a

model_2 :: Num r => r -> r -> r -> r
model_2 x a b = a*x + b

model_3 :: Num r => r -> r -> r -> r -> r
model_3 x a b c = a*x*x + b*x + c

model_4 :: Num r => r -> r -> r -> r -> r -> r
model_4 x a b c d = a*x*x*x + b*x*x + c*x + d

rndGenSeed :: Int
rndGenSeed = 123456

test :: (Show a, Nat n)
     => (Model a Double n)
     -> SizedList Double n
     -> [a]
     -> Double
     -> (SizedList Double n, Info Double, CovarMatrix Double n)
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
      samples = zip xs $ map (\x -> f x $* ps) xs
      samples' = zipWith (\(x, y) n -> (x, y + (n - 0.5) * 2 * noise)) samples ns

main :: IO ()
main = print $ test model_3 (5 ::: 3 ::: 7 ::: Nil) [1..50] 2

