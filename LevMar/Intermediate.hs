module LevMar.Intermediate
    ( Model
    , Options(..)
    , StopReason(..)
    , Info(..)
    , CovarMatrix
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
import Foreign.Ptr           (nullPtr, plusPtr)
import Foreign.Storable      (Storable)
import System.IO.Unsafe      (unsafePerformIO)

import qualified LevMar.Binding as C_LMA

type Model a r = a -> [r] -> r

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

type CovarMatrix r = [[r]]

type LevMarDif r a =  Model a r
                   -> [r]       -- initial parameters
                   -> [(a, r)]  -- samples
                   -> Integer   -- itmax
                   -> Options r -- opts
                   -> ([r], Info r, CovarMatrix r)

defaultOpts :: Fractional r => Options r
defaultOpts = Opts { opt_mu       = C_LMA._LM_INIT_MU
                   , opt_epsilon1 = C_LMA._LM_STOP_THRESH
                   , opt_epsilon2 = C_LMA._LM_STOP_THRESH
                   , opt_epsilon3 = C_LMA._LM_STOP_THRESH
                   , opt_delta    = C_LMA._LM_DIFF_DELTA
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

convertModel :: (Real r, Fractional r, Storable c, Real c, Fractional c)
             => [a] -> Model a r -> C_LMA.Model c
convertModel xs f = \parPtr hxPtr numPar _ _ -> do
                      params <- peekArray (fromIntegral numPar) parPtr
                      pokeArray hxPtr $
                        map (\x -> realToFrac $ f x $ map realToFrac params) xs

gen_levmar_dif :: (Storable cr, Real cr, Fractional cr, RealFrac r)
               => C_LMA.LevMarDif cr -> LevMarDif r a
gen_levmar_dif lma f ps samples itMax opts = unsafePerformIO $
    let lenPs    = length ps
        coVarLen = lenPs * lenPs
        (xs, ys) = unzip samples
    in withArray (map realToFrac ps) $ \psPtr ->
         withArray (map realToFrac ys) $ \ysPtr ->
           withArray (map realToFrac $ optsToList opts) $ \optsPtr ->
             allocaArray C_LMA._LM_INFO_SZ $ \infoPtr ->
               allocaArray coVarLen $ \coVarPtr ->
                 C_LMA.withModel (convertModel xs f) $ \modelPtr -> do

                   r <- lma modelPtr
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
                   info   <- peekArray C_LMA._LM_INFO_SZ infoPtr

                   let coVarPtrEnd = plusPtr coVarPtr coVarLen
                   let mkCoVarMatrix ptr
                           | ptr == coVarPtrEnd = return []
                           | otherwise = do row <- peekArray lenPs ptr
                                            rows <- mkCoVarMatrix $ plusPtr ptr lenPs
                                            return (row : rows)

                   coVar <- mkCoVarMatrix coVarPtr

                   return $ ( map realToFrac result
                            , listToInfo $ map realToFrac info
                            , map (map realToFrac) coVar
                            )

dlevmar_dif :: LevMarDif Double a
dlevmar_dif = gen_levmar_dif C_LMA.dlevmar_dif

slevmar_dif :: LevMarDif Float a
slevmar_dif = gen_levmar_dif C_LMA.slevmar_dif
