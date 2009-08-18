module Main where

import LevMar

import Foreign.Marshal.Array
import Foreign.C.Types
import Foreign.Ptr

samples :: [CDouble]
samples = [ 2.90, 3.40, 3.10, 2.79, 3.15
          , 3.37, 2.87, 3.67, 2.57, 2.78
          , 2.96, 3.03, 3.34, 2.87, 2.20
          ]

params :: [CDouble]
params = [0]

model :: Model CDouble
model parPtr m pd md _ = do [p] <- peekArray (fromIntegral pd) parPtr
                            pokeArray m $ replicate (fromIntegral md) p

main :: IO ()
main = withArray params  $ \paramsPtr  ->
         withArray samples $ \samplesPtr -> do
           modelPtr <- mkModel model
           i <- c_dlevmar_dif modelPtr
                              paramsPtr
                              samplesPtr
                              (fromIntegral $ length params)
                              (fromIntegral $ length samples)
                              1000
                              nullPtr
                              nullPtr
                              nullPtr
                              nullPtr
                              nullPtr
           freeHaskellFunPtr modelPtr
           putStrLn $ "Stopped after " ++ (show i) ++ " iterations"
           result <- peekArray (length params) paramsPtr
           putStrLn $ "Result: " ++ (show result)
