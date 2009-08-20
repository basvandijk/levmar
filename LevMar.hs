module LevMar ( Model
              , Options(..)
              , StopReason(..)
              , Info(..)
              , CoVarMatrix
              , LevMarDif
              , defaultOpts
              , dlevmar_dif
              , slevmar_dif
              ) where

import Foreign.Marshal.Array ( allocaArray
                             , peekArray
                             , pokeArray
                             , withArray
                             )
import Foreign.Ptr           (nullPtr)
import Foreign.Storable      (Storable)

import qualified Bindings.LevMar as C_LMA


type Model r a = [r] -> a -> r

data Options r = Opts { opt_mu       :: r
                      , opt_epsilon1 :: r
                      , opt_epsilon2 :: r
                      , opt_epsilon3 :: r
                      , opt_delta    :: r
                      } deriving Show

data StopReason = SmallGradient
                | SmallDp
                | MaxIterations
                | SingularMatrix
                | SmallestError
                | SmallE_2
                | InvalidValues
                  deriving (Show, Enum)

data Info r = Info { inf_values          :: [r]
                   , inf_numIter         :: Integer
                   , inf_stopReason      :: StopReason
                   , inf_numFuncEvals    :: Integer
                   , inf_numJacobEvals   :: Integer
                   , inf_numLinSysSolved :: Integer
                   } deriving Show

type CoVarMatrix r = [r]

type LevMarDif r a =  Model r a
                   -> [r]       -- initial parameters
                   -> [(a, r)]  -- samples
                   -> Integer   -- itmax
                   -> Options r -- opts
                   -> IO ([r], Info r, CoVarMatrix r)

defaultOpts :: Fractional r => Options r
defaultOpts = Opts { opt_mu       = 1e-3
                   , opt_epsilon1 = 1e-17
                   , opt_epsilon2 = 1e-17
                   , opt_epsilon3 = 1e-17
                   , opt_delta    = 1e-6
                   }

optsToList :: Options r -> [r]
optsToList (Opts mu eps1 eps2 eps3 delta) = [mu, eps1, eps2, eps3, delta]

listToInfo :: RealFrac r => [r] -> Info r
listToInfo [a,b,c,d,e,f,g,h,i,j] =
    Info { inf_values          = [a,b,c,d,e]
         , inf_numIter         = floor f
         , inf_stopReason      = toEnum $ floor g - 1
         , inf_numFuncEvals    = floor h
         , inf_numJacobEvals   = floor i
         , inf_numLinSysSolved = floor j
         }
listToInfo _ = error "liftToInfo: wrong list length"

convertModel :: (Real h, Fractional h, Storable c, Real c, Fractional c)
             => [a] -> Model h a -> C_LMA.Model c
convertModel xs f = \parPtr hxPtr numPar _ _ -> do
                      params <- peekArray (fromIntegral numPar) parPtr
                      pokeArray hxPtr $
                        map (realToFrac . f (map realToFrac params)) xs

gen_levmar_dif :: (Storable cr, Real cr, Fractional cr, RealFrac r)
               => C_LMA.LevMarDif cr -> LevMarDif r a
gen_levmar_dif minimise f ps samples itMax opts =
    let lenPs    = length ps
        (xs, ys) = unzip samples
    in withArray (map realToFrac ps) $ \psPtr ->
         withArray (map realToFrac ys) $ \ysPtr ->
           withArray (map realToFrac $ optsToList opts) $ \optsPtr ->
             allocaArray 10 $ \infoPtr ->
               allocaArray (lenPs * lenPs) $ \coVarPtr ->
                 C_LMA.withModel (convertModel xs f) $ \modelPtr -> do

                   _ <- minimise modelPtr
                                 psPtr
                                 ysPtr
                                 (fromIntegral $ lenPs)
                                 (fromIntegral $ length samples)
                                 (fromIntegral itMax)
                                 optsPtr
                                 infoPtr
                                 nullPtr -- work
                                 coVarPtr
                                 nullPtr -- adata

                   result <- peekArray lenPs psPtr
                   info   <- peekArray 10 infoPtr
                   coVar  <- peekArray (lenPs * lenPs) coVarPtr

                   return $ ( map realToFrac result
                            , listToInfo $ map realToFrac info
                            , map realToFrac coVar
                            )

dlevmar_dif :: LevMarDif Double a
dlevmar_dif = gen_levmar_dif C_LMA.dlevmar_dif

slevmar_dif :: LevMarDif Float a
slevmar_dif = gen_levmar_dif C_LMA.slevmar_dif
