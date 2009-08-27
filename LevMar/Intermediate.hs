{-# LANGUAGE ScopedTypeVariables #-}

module LevMar.Intermediate
    ( Model
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
import Foreign.C.Types       (CInt)
import System.IO.Unsafe      (unsafePerformIO)
import Data.Maybe            (fromJust, isJust)

import qualified Bindings.LevMar.CurryFriendly as LMA_C

type Model    r a = [r] -> a -> r
type Jacobian r a = [r] -> a -> [r]

data Options r = Opts { optMu       :: r
                      , optEpsilon1 :: r
                      , optEpsilon2 :: r
                      , optEpsilon3 :: r
                      , optDelta    :: r
                      } deriving Show

data StopReason = SmallGradient
                | SmallDp
                | MaxIterations
                | SingularMatrix
                | SmallestError
                | SmallE_2
                | InvalidValues
                  deriving (Show, Enum)

data Info r = Info { infValues          :: [r]
                   , infNumIter         :: Integer
                   , infStopReason      :: StopReason
                   , infNumFuncEvals    :: Integer
                   , infNumJacobEvals   :: Integer
                   , infNumLinSysSolved :: Integer
                   } deriving Show

type CovarMatrix r = [[r]]

defaultOpts :: Fractional r => Options r
defaultOpts = Opts { optMu       = LMA_C._LM_INIT_MU
                   , optEpsilon1 = LMA_C._LM_STOP_THRESH
                   , optEpsilon2 = LMA_C._LM_STOP_THRESH
                   , optEpsilon3 = LMA_C._LM_STOP_THRESH
                   , optDelta    = LMA_C._LM_DIFF_DELTA
                   }

optsToList :: Options r -> [r]
optsToList (Opts mu eps1 eps2 eps3 delta) = [mu, eps1, eps2, eps3, delta]

listToInfo :: RealFrac r => [r] -> Info r
listToInfo [a,b,c,d,e,f,g,h,i,j] =
    Info { infValues          = [a,b,c,d,e]
         , infNumIter         = floor f
         , infStopReason      = toEnum $ floor g - 1
         , infNumFuncEvals    = floor h
         , infNumJacobEvals   = floor i
         , infNumLinSysSolved = floor j
         }
listToInfo _ = error "liftToInfo: wrong list length"

convertModel :: (Real r, Fractional r, Storable c, Real c, Fractional c)
             => [a] -> Model r a -> LMA_C.Model c
convertModel xs f = \parPtr hxPtr numPar _ _ -> do
                      params <- peekArray (fromIntegral numPar) parPtr
                      pokeArray hxPtr $
                        map (realToFrac . f (map realToFrac params)) xs

convertJacobian :: (Real r, Fractional r, Storable c, Real c, Fractional c)
                => [a] -> Jacobian r a -> LMA_C.Jacobian c
convertJacobian xs f = \parPtr jPtr numPar _ _ -> do
                         params <- peekArray (fromIntegral numPar) parPtr
                         let params' = map realToFrac params
                         pokeArray jPtr $
                           concatMap (map realToFrac . f params') xs

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

type LevMar r a =  Model r a
                -> Maybe (Jacobian r a)
                -> [r]                -- initial parameters
                -> [(a, r)]           -- samples
                -> Integer            -- itmax
                -> Options r          -- opts
                -> Maybe [r]          -- lower bounds
                -> Maybe [r]          -- upper bounds
                -> Maybe ([[r]], [r]) -- linear constraints
                -> Maybe [r]          -- weights
                -> ([r], Info r, CovarMatrix r)


{-

Preconditions:
  length samples >= length ps

     isJust mLowBs && length (fromJust mLowBs) == length ps
  && isJust mUpBs  && length (fromJust mUpBs)  == length ps
-}
gen_levmar :: forall cr r a. (Storable cr, Real cr, Fractional cr, RealFrac r)
           => LMA_C.LevMarDer cr
           -> LMA_C.LevMarDif cr
           -> LMA_C.LevMarBCDer cr
           -> LMA_C.LevMarBCDif cr
           -> LMA_C.LevMarLecDer cr
           -> LMA_C.LevMarLecDif cr
           -> LMA_C.LevMarBLecDer cr
           -> LMA_C.LevMarBLecDif cr
           -> LevMar r a
gen_levmar f_der f_dif f_bc_der f_bc_dif f_lec_der f_lec_dif f_blec_der f_blec_dif
           model mJac ps samples itMax opts mLowBs mUpBs mLinC mWghts
    | lenSamples < lenPs = error "gen_levmar: not enough samples"
    | otherwise = unsafePerformIO $
        withArray (map realToFrac ps) $ \psPtr ->
        withArray (map realToFrac ys) $ \ysPtr ->
        withArray (map realToFrac $ optsToList opts) $ \optsPtr ->
        allocaArray LMA_C._LM_INFO_SZ $ \infoPtr ->
        allocaArray workSize $ \workPtr ->
        allocaArray covarLen $ \covarPtr ->
        LMA_C.withModel (convertModel xs model) $ \modelPtr -> do
          let runDif :: LMA_C.LevMarDif cr -> IO CInt
              runDif f = f modelPtr
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
                                 let runDer :: LMA_C.LevMarDer cr -> IO CInt
                                     runDer f = runDif $ f jacobPtr
                                 in if boxConstrained
                                    then runDer `withBoxConstraints` f_bc_der
                                    else runDer f_der

                   Nothing -> if boxConstrained
                              then runDif `withBoxConstraints` f_bc_dif
                              else runDif f_dif
               )

          result <- peekArray lenPs psPtr
          info   <- peekArray LMA_C._LM_INFO_SZ infoPtr

          let covarPtrEnd = plusPtr covarPtr covarLen
          let mkCovarMatrix ptr | ptr == covarPtrEnd = return []
                                | otherwise = do row <- peekArray lenPs ptr
                                                 rows <- mkCovarMatrix $ plusPtr ptr lenPs
                                                 return (row : rows)
          covar <- mkCovarMatrix covarPtr

          return ( map realToFrac result
                 , listToInfo $ map realToFrac info
                 , map (map realToFrac) covar
                 )
    where
      lenPs      = length ps
      lenSamples = length samples
      covarLen   = lenPs * lenPs
      (xs, ys)   = unzip samples
      (cMat, rhsVec) = fromJust mLinC
      numConstraints = length cMat

      workSize = let p | linConstrained = lenPs - numConstraints
                       | otherwise      = lenPs
                     s |  linConstrained
                       && boxConstrained = lenSamples + lenPs
                       | otherwise       = lenSamples
                     a | isJust mJac         = 2
                       |  not linConstrained
                       && boxConstrained     = 2
                       | otherwise           = 4
                 in a*s + 4*p + s*p + p*p

      linConstrained = isJust mLinC
      boxConstrained = isJust mLowBs || isJust mUpBs

      withBoxConstraints f g = maybeWithArray ((fmap . fmap) realToFrac mLowBs) $ \lBsPtr ->
                                 maybeWithArray ((fmap . fmap) realToFrac mUpBs) $ \uBsPtr ->
                                   f $ g lBsPtr uBsPtr


class LevMarable r where
    levmar :: LevMar r a

instance LevMarable Float where
    levmar = gen_levmar LMA_C.slevmar_der
                        LMA_C.slevmar_dif
                        LMA_C.slevmar_bc_der
                        LMA_C.slevmar_bc_dif
                        LMA_C.slevmar_lec_der
                        LMA_C.slevmar_lec_dif
                        LMA_C.slevmar_blec_der
                        LMA_C.slevmar_blec_dif

instance LevMarable Double where
    levmar = gen_levmar LMA_C.dlevmar_der
                        LMA_C.dlevmar_dif
                        LMA_C.dlevmar_bc_der
                        LMA_C.dlevmar_bc_dif
                        LMA_C.dlevmar_lec_der
                        LMA_C.dlevmar_lec_dif
                        LMA_C.dlevmar_blec_der
                        LMA_C.dlevmar_blec_dif
