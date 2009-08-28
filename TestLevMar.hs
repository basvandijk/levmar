module Main where

import LevMar
import qualified SizedList as SL
import NFunction(($*))
import System.Random

type N0 = Z
type N1 = S N0
type N2 = S N1
type N3 = S N2
type N4 = S N3
type N5 = S N4


model1 :: Model N1 r a
model1 a _ = a

model2 :: Num r => Model N2 r r
model2 a b x = a*x + b

model3 :: Num r => Model N3 r r
model3 a b c x = a*x*x + b*x + c

model4 :: Num r => Model N4 r r
model4 a b c d x = a*x*x*x + b*x*x + c*x + d

rndGenSeed :: Int
rndGenSeed = 123456

test :: (Show a, Nat n)
     => (Model n Double a)
     -> SizedList n Double
     -> [a]
     -> Double
     -> Maybe (SizedList n Double, Info Double, CovarMatrix n Double)
test f ps xs noise = levmar f
                            Nothing
                            (SL.replicate (SL.length ps) 0)
                            samples'
                            1000
                            defaultOpts
                            Nothing
                            Nothing
                            noLinearConstraints
                            Nothing
    where
      ns = take (length xs) $ randoms $ mkStdGen rndGenSeed
      samples = zip xs $ map (f $* ps) xs
      samples' = zipWith (\(x, y) n -> (x, y + (n - 0.5) * 2 * noise)) samples ns

-- main :: IO ()
-- main = print $ test expfunc (1 ::: 0 ::: 0 ::: Nil) [1..40] 2

main :: IO ()
main = print testExpFunc

expFunc :: Floating r => Model N3 r r
expFunc a b c x = a * exp ((-b) * x) + c

jacExpFunc :: Floating r => Jacobian N3 r r
jacExpFunc a b _ x = let u = exp ((-b) * x)
                     in u ::: (-a) * x * u ::: 1.0 ::: Nil

testExpFunc :: Maybe (SizedList N3 Double, Info Double, CovarMatrix N3 Double)
testExpFunc = levmar expFunc
                     (Just jacExpFunc)
                     initParams
                     samples'
                     1000
                     ( defaultOpts { optEpsilon1 = 1e-15
                                   , optEpsilon2 = 1e-15
                                   , optEpsilon3 = 1e-20
                                   }
                     )
                     Nothing
                     Nothing
                     noLinearConstraints
                     Nothing
    where
      ns = take (length xs) $ randoms $ mkStdGen rndGenSeed
      xs = [1..40]

      initParams, actualParams :: SizedList N3 Double
      initParams   = 1.0 ::: 0.0 ::: 0.0 ::: Nil
      actualParams = 5.0 ::: 0.1 ::: 1.0 ::: Nil

      samples  = zip xs $ map (expFunc $* actualParams) xs
      samples' = zipWith (\(x, y) n -> (x, y + (n - 0.5) * 2 * 0.1)) samples ns
