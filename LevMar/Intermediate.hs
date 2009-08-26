module LevMar.Intermediate
    ( ParamFunc
    , Model
    , Jacobian
    , Options(..)
    , StopReason(..)
    , Info(..)
    , CovarMatrix

    , defaultOpts

    , LevMarable
    , levmar
    , LevMar
    ) where

import Foreign.Marshal.Array (allocaArray, peekArray, pokeArray, withArray)
import Foreign.Ptr           (Ptr, nullPtr, plusPtr)
import Foreign.Storable      (Storable)
import System.IO.Unsafe      (unsafePerformIO)
import Data.Maybe            (fromJust, isJust)

import qualified LevMar.Binding.CurryFriendly as LMA_C

type ParamFunc  r = [r] -> r
type Model    a r = a -> ParamFunc r
type Jacobian a r = Model a r

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

defaultOpts :: Fractional r => Options r
defaultOpts = Opts { opt_mu       = LMA_C._LM_INIT_MU
                   , opt_epsilon1 = LMA_C._LM_STOP_THRESH
                   , opt_epsilon2 = LMA_C._LM_STOP_THRESH
                   , opt_epsilon3 = LMA_C._LM_STOP_THRESH
                   , opt_delta    = LMA_C._LM_DIFF_DELTA
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
             => [a] -> Model a r -> LMA_C.Model c
convertModel xs f = \parPtr hxPtr numPar _ _ -> do
                      params <- peekArray (fromIntegral numPar) parPtr
                      pokeArray hxPtr $
                        map (\x -> realToFrac $ f x $ map realToFrac params) xs

convertJacobian :: (Real r, Fractional r, Storable c, Real c, Fractional c)
                => [a] -> Jacobian a r -> LMA_C.Jacobian c
convertJacobian = convertModel

{-
-- All smaller than
(<*) :: Ord a => [a] -> [a] -> Bool
xs <* ys = and $ zipWith (<) xs ys

type LevMarBLecDer r a =  Model a r
                       -> Jacobian a r
                       -> [r]       -- initial parameters
                       -> [(a, r)]  -- samples
                       -> Integer   -- itmax
                       -> Options r -- opts
                       -> [[r]]     -- constraints matrix
                       -> [r]       -- right hand constraints vector
                       -> ([r], Info r, CovarMatrix r)

gen_levmar_blec_der :: (Storable cr, Real cr, Fractional cr, RealFrac r)
                   => LMA_C.LevMarBLecDer cr -> LevMarBLecDer r a
gen_levmar_blec_der lma f j ps samples itMax opts cMat rhcVec
    | lenSamples < lenPs + lenCMat = notEnoughSamplesError "gen_levmar_lec_der" -- TODO: mention linear constraints
    | otherwise = unsafePerformIO $
    withArray (map realToFrac ps) $ \psPtr ->
      withArray (map realToFrac ys) $ \ysPtr ->
        withArray (map realToFrac $ optsToList opts) $ \optsPtr ->
          withArray (map realToFrac $ concat cMat) $ \cMatPtr ->
            withArray (map realToFrac $ rhcVec) $ \rhcVecPtr ->
              allocaArray (workSize lenPs lenSamples) $ \workPtr ->
                allocaArray LMA_C._LM_INFO_SZ $ \infoPtr ->
                  allocaArray covarLen $ \covarPtr ->
                    LMA_C.withModel (convertModel xs f) $ \modelPtr ->
                      LMA_C.withJacobian (convertJacobian xs j) $ \jacobPtr -> do

                        _ <- lma modelPtr
                                 jacobPtr
                                 psPtr
                                 ysPtr
                                 (fromIntegral lenPs)
                                 (fromIntegral lenSamples)
                                 cMatPtr
                                 rhcVecPtr
                                 (fromIntegral lenCMat)
                                 (fromIntegral itMax)
                                 optsPtr
                                 infoPtr
                                 workPtr
                                 covarPtr
                                 nullPtr -- adata

                        result <- peekArray lenPs psPtr
                        info   <- peekArray LMA_C._LM_INFO_SZ infoPtr

                        let covarPtrEnd = plusPtr covarPtr covarLen
                        let mkCovarMatrix ptr
                                | ptr == covarPtrEnd = return []
                                | otherwise = do row <- peekArray lenPs ptr
                                                 rows <- mkCovarMatrix $ plusPtr ptr lenPs
                                                 return (row : rows)

                        covar <- mkCovarMatrix covarPtr

                        return $ ( map realToFrac result
                                 , listToInfo $ map realToFrac info
                                 , map (map realToFrac) covar
                                 )
  where lenSamples   = length samples
        lenPs        = length ps
        lenCMat      = length cMat
        covarLen     = lenPs * lenPs
        (xs, ys)     = unzip samples
        workSize n m = 2*m + 4*n + m*n + n*n
-}

-------------------------------------------------------------------------------
-- All-in-one LMA function
-------------------------------------------------------------------------------

maybeWithArray :: Storable a => Maybe [a] -> (Ptr a -> IO b) -> IO b
maybeWithArray Nothing   f = f nullPtr
maybeWithArray (Just xs) f = withArray xs f

type LevMar r a =  Model a r
                -> Maybe (Jacobian a r)
                -> [r]       -- initial parameters
                -> [(a, r)]  -- samples
                -> Integer   -- itmax
                -> Options r -- opts
                -> Maybe [r] -- lower bounds
                -> Maybe [r] -- upper bounds
                -> ([r], Info r, CovarMatrix r)

gen_levmar :: (Storable cr, Real cr, Fractional cr, RealFrac r)
       => LMA_C.LevMarDer cr
       -> LMA_C.LevMarDif cr
       -> LMA_C.LevMarBCDer cr
       -> LMA_C.LevMarBCDif cr
       -> LevMar r a
gen_levmar lma_der lma_dif lma_bc_der lma_bc_dif
           model mJac ps samples itMax opts mLowBs mUpBs
    | lenSamples < lenPs = error "levmar: not enough samples"
    | bcError            = error "levmar: wrong number of box constraints"
    | otherwise = unsafePerformIO $
        withArray (map realToFrac ps) $ \psPtr ->
        withArray (map realToFrac ys) $ \ysPtr ->
        withArray (map realToFrac $ optsToList opts) $ \optsPtr ->
        allocaArray LMA_C._LM_INFO_SZ $ \infoPtr ->
        allocaArray workSize $ \workPtr ->
        allocaArray covarLen $ \covarPtr ->
        LMA_C.withModel (convertModel xs model) $ \modelPtr -> do
          let lma f = f modelPtr
                        psPtr
                        ysPtr
                        (fromIntegral lenPs)
                        (fromIntegral lenSamples)
                        (fromIntegral itMax)
                        optsPtr
                        infoPtr
                        workPtr
                        covarPtr
                        nullPtr

          _ <- ( case mJac of
                   Just jac -> LMA_C.withJacobian (convertJacobian xs jac) $ \jacobPtr ->
                                 if boxConstrained
                                 then withBoxConstraints $ \lBsPtr uBsPtr -> lma $ lma_bc_der lBsPtr uBsPtr jacobPtr
                                 else lma $ lma_der jacobPtr

                   Nothing -> if boxConstrained
                              then withBoxConstraints $ \lBsPtr uBsPtr -> lma $ lma_bc_dif lBsPtr uBsPtr
                              else lma lma_dif
               )

          result <- peekArray lenPs psPtr
          info   <- peekArray LMA_C._LM_INFO_SZ infoPtr

          let covarPtrEnd = plusPtr covarPtr covarLen
          let mkCovarMatrix ptr | ptr == covarPtrEnd = return []
                                | otherwise = do row <- peekArray lenPs ptr
                                                 rows <- mkCovarMatrix $ plusPtr ptr lenPs
                                                 return (row : rows)
          covar <- mkCovarMatrix covarPtr

          return $ ( map realToFrac result
                   , listToInfo $ map realToFrac info
                   , map (map realToFrac) covar
                   )
    where
      lenPs      = length ps
      lenSamples = length samples
      lBs        = fromJust mLowBs
      uBs        = fromJust mUpBs
      covarLen   = lenPs * lenPs
      (xs, ys)   = unzip samples

      workSize = let n = lenPs
                     m = lenSamples
                     a | boxConstrained = 2
                       | isJust mJac    = 2
                       | otherwise      = 4
                 in a*m + 4*n + m*n + n*n

      bcError | isJust mLowBs && length lBs /= lenPs = True
              | isJust mUpBs  && length uBs /= lenPs = True
              | otherwise = False

      boxConstrained | isJust mLowBs = True
                     | isJust mUpBs  = True
                     | otherwise     = False

      withBoxConstraints f = maybeWithArray ((fmap . fmap) realToFrac mLowBs) $ \lBsPtr ->
                               maybeWithArray ((fmap . fmap) realToFrac mUpBs) $ \uBsPtr ->
                                 f lBsPtr uBsPtr

slevmar :: LevMar Float a
slevmar = gen_levmar LMA_C.slevmar_der
                     LMA_C.slevmar_dif
                     LMA_C.slevmar_bc_der
                     LMA_C.slevmar_bc_dif

dlevmar :: LevMar Double a
dlevmar = gen_levmar LMA_C.dlevmar_der
                     LMA_C.dlevmar_dif
                     LMA_C.dlevmar_bc_der
                     LMA_C.dlevmar_bc_dif

class LevMarable r where
    levmar :: LevMar r a

instance LevMarable Float where
    levmar = slevmar

instance LevMarable Double where
    levmar = dlevmar
