{-# LANGUAGE ScopedTypeVariables #-}

module LevMar.Intermediate
    ( levmar
    , Model
    , Jacobian

    , Options(..)
    , defaultOpts

    , Info(..)
    , StopReason(..)
    , CovarMatrix

    , LevMarable
    , levmar'
    , LevMar
    , Model'
    , Jacobian'
    ) where

import Foreign.Marshal.Array (allocaArray, peekArray, pokeArray, withArray)
import Foreign.Ptr           (Ptr, nullPtr, plusPtr)
import Foreign.Storable      (Storable)
import Foreign.C.Types       (CInt)
import System.IO.Unsafe      (unsafePerformIO)
import Data.Maybe            (fromJust, isJust)

import qualified Bindings.LevMar.CurryFriendly as LMA_C

levmar :: LevMarable r
       => Model r a
       -> Maybe (Jacobian r a)
       -> [r]                -- initial parameters
       -> [(a, r)]           -- samples
       -> Integer            -- itmax
       -> Options r          -- opts
       -> Maybe [r]          -- lower bounds
       -> Maybe [r]          -- upper bounds
       -> Maybe ([[r]], [r]) -- linear constraints
       -> Maybe [r]          -- weights
       -> Maybe ([r], Info r, CovarMatrix r)
levmar model mJac ps samples =
    levmar' (\ps -> map (model ps) xs)
            (fmap (\jac -> \ps -> map (jac ps) xs) mJac)
            ps
            ys
        where
          (xs, ys) = unzip samples

type Model  r a = [r] -> a -> r

type Jacobian r a = [r] -> a -> [r]

data Options r = Opts { optMu       :: r
                      , optEpsilon1 :: r
                      , optEpsilon2 :: r
                      , optEpsilon3 :: r
                      , optDelta    :: r
                      } deriving Show

defaultOpts :: Fractional r => Options r
defaultOpts = Opts { optMu       = LMA_C._LM_INIT_MU
                   , optEpsilon1 = LMA_C._LM_STOP_THRESH
                   , optEpsilon2 = LMA_C._LM_STOP_THRESH
                   , optEpsilon3 = LMA_C._LM_STOP_THRESH
                   , optDelta    = LMA_C._LM_DIFF_DELTA
                   }

optsToList :: Options r -> [r]
optsToList (Opts mu  eps1  eps2  eps3  delta) =
                [mu, eps1, eps2, eps3, delta]

data Info r = Info { infValues          :: [r]
                   , infNumIter         :: Integer
                   , infStopReason      :: StopReason
                   , infNumFuncEvals    :: Integer
                   , infNumJacobEvals   :: Integer
                   , infNumLinSysSolved :: Integer
                   } deriving Show

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

data StopReason = SmallGradient
                | SmallDp
                | MaxIterations
                | SingularMatrix
                | SmallestError
                | SmallE_2
                | InvalidValues
                  deriving (Show, Enum)

type CovarMatrix r = [[r]]

type LevMar r =  Model' r
              -> Maybe (Jacobian' r)
              -> [r]                -- initial parameters
              -> [r]                -- samples
              -> Integer            -- itmax
              -> Options r          -- opts
              -> Maybe [r]          -- lower bounds
              -> Maybe [r]          -- upper bounds
              -> Maybe ([[r]], [r]) -- linear constraints
              -> Maybe [r]          -- weights
              -> Maybe ([r], Info r, CovarMatrix r)

type Model' r = [r] -> [r]

type Jacobian' r = [r] -> [[r]]

class LevMarable r where
    levmar' :: LevMar r

instance LevMarable Float where
    levmar' = gen_levmar LMA_C.slevmar_der
                         LMA_C.slevmar_dif
                         LMA_C.slevmar_bc_der
                         LMA_C.slevmar_bc_dif
                         LMA_C.slevmar_lec_der
                         LMA_C.slevmar_lec_dif
                         LMA_C.slevmar_blec_der
                         LMA_C.slevmar_blec_dif

instance LevMarable Double where
    levmar' = gen_levmar LMA_C.dlevmar_der
                         LMA_C.dlevmar_dif
                         LMA_C.dlevmar_bc_der
                         LMA_C.dlevmar_bc_dif
                         LMA_C.dlevmar_lec_der
                         LMA_C.dlevmar_lec_dif
                         LMA_C.dlevmar_blec_der
                         LMA_C.dlevmar_blec_dif

{-|
Preconditions:
  length ys >= length ps

     isJust mLowBs && length (fromJust mLowBs) == length ps
  && isJust mUpBs  && length (fromJust mUpBs)  == length ps

  boxConstrained && (all $ zipWith (<=) (fromJust mLowBs) (fromJust mUpBs))
-}
gen_levmar :: forall cr r. (Storable cr, Real cr, Fractional cr, RealFrac r)
           => LMA_C.LevMarDer cr
           -> LMA_C.LevMarDif cr
           -> LMA_C.LevMarBCDer cr
           -> LMA_C.LevMarBCDif cr
           -> LMA_C.LevMarLecDer cr
           -> LMA_C.LevMarLecDif cr
           -> LMA_C.LevMarBLecDer cr
           -> LMA_C.LevMarBLecDif cr
           -> LevMar r
gen_levmar f_der
           f_dif
           f_bc_der
           f_bc_dif
           f_lec_der
           f_lec_dif
           f_blec_der
           f_blec_dif
           model mJac ps ys itMax opts mLowBs mUpBs mLinC mWeights
    | lenYs < lenPs = Nothing
    | otherwise = unsafePerformIO $
        withArray (map realToFrac ps) $ \psPtr ->
        withArray (map realToFrac ys) $ \ysPtr ->
        withArray (map realToFrac $ optsToList opts) $ \optsPtr ->
        allocaArray LMA_C._LM_INFO_SZ $ \infoPtr ->
        allocaArray covarLen $ \covarPtr ->
        LMA_C.withModel (convertModel model) $ \modelPtr -> do
          let runDif :: LMA_C.LevMarDif cr -> IO CInt
              runDif f = f modelPtr
                           psPtr
                           ysPtr
                           (fromIntegral lenPs)
                           (fromIntegral lenYs)
                           (fromIntegral itMax)
                           optsPtr
                           infoPtr
                           nullPtr
                           covarPtr
                           nullPtr

          r <- case mJac of
                 Just jac -> LMA_C.withJacobian (convertJacobian jac) $ \jacobPtr ->
                               let runDer :: LMA_C.LevMarDer cr -> IO CInt
                                   runDer f = runDif $ f jacobPtr
                               in if boxConstrained
                                  then if linConstrained
                                       then withBoxConstraints (withLinConstraints $ withWeights runDer) f_blec_der
                                       else withBoxConstraints runDer f_bc_der
                                  else if linConstrained
                                       then withLinConstraints runDer f_lec_der
                                       else runDer f_der

                 Nothing -> if boxConstrained
                            then if linConstrained
                                 then withBoxConstraints (withLinConstraints $ withWeights runDif) f_blec_dif
                                 else withBoxConstraints runDif f_bc_dif
                            else if linConstrained
                                 then withLinConstraints runDif f_lec_dif
                                 else runDif f_dif

          if r == LMA_C._LM_ERROR
            then return Nothing
            else do result <- peekArray lenPs psPtr
                    info   <- peekArray LMA_C._LM_INFO_SZ infoPtr

                    let covarPtrEnd = plusPtr covarPtr covarLen
                    let convertCovarMatrix ptr
                            | ptr == covarPtrEnd = return []
                            | otherwise = do row <- peekArray lenPs ptr
                                             rows <- convertCovarMatrix $ plusPtr ptr lenPs
                                             return $ row : rows

                    covar  <- convertCovarMatrix covarPtr

                    return $ Just ( map realToFrac result
                                  , listToInfo $ map realToFrac info
                                  , map (map realToFrac) covar
                                  )
    where
      lenPs          = length ps
      lenYs          = length ys
      covarLen       = lenPs * lenPs
      (cMat, rhcVec) = fromJust mLinC

      -- Whether the parameters are constrained by a linear equation.
      linConstrained = isJust mLinC

      -- Whether the parameters are constrained by a bounding box.
      boxConstrained = isJust mLowBs || isJust mUpBs

      withBoxConstraints f g = maybeWithArray ((fmap . fmap) realToFrac mLowBs) $ \lBsPtr ->
                                 maybeWithArray ((fmap . fmap) realToFrac mUpBs) $ \uBsPtr ->
                                   f $ g lBsPtr uBsPtr

      withLinConstraints f g = withArray (map realToFrac $ concat cMat) $ \cMatPtr ->
                                 withArray (map realToFrac rhcVec) $ \rhcVecPtr ->
                                   f $ g cMatPtr rhcVecPtr $ fromIntegral $ length cMat

      withWeights f g = maybeWithArray ((fmap . fmap) realToFrac mWeights) $ \weightsPtr ->
                          f $ g weightsPtr

convertModel :: (Real r, Fractional r, Storable c, Real c, Fractional c)
             =>  Model' r -> LMA_C.Model c
convertModel model = \parPtr hxPtr numPar _ _ -> do
                       params <- peekArray (fromIntegral numPar) parPtr
                       pokeArray hxPtr $ map realToFrac $ model $ map realToFrac params

convertJacobian :: (Real r, Fractional r, Storable c, Real c, Fractional c)
                => Jacobian' r -> LMA_C.Jacobian c
convertJacobian jac = \parPtr jPtr numPar _ _ -> do
                        params <- peekArray (fromIntegral numPar) parPtr
                        pokeArray jPtr $ concatMap (map realToFrac) $ jac $ map realToFrac params

maybeWithArray :: Storable a => Maybe [a] -> (Ptr a -> IO b) -> IO b
maybeWithArray Nothing   f = f nullPtr
maybeWithArray (Just xs) f = withArray xs f
