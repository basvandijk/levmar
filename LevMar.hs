{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ScopedTypeVariables #-}

module LevMar
    ( Model
    , Jacobian
    , levmar

      -- *Supporting types and default values
    , LMA_I.Options(..)
    , LMA_I.defaultOpts
    , LMA_I.StopReason(..)
    , LMA_I.Info(..)
    , LMA_I.LevMarable
    , LMA_I.LevMarError(..)
    , CovarMatrix
    , LinearConstraints
    , noLinearConstraints
    , Matrix

      -- *Type-level stuff
    , Z, S, Nat
    , SizedList(..)
    , NFunction
    )
    where

import qualified LevMar.Intermediate as LMA_I

import TypeLevelNat (Z, S, Nat)
import SizedList    (SizedList(..), toList, unsafeFromList)
import NFunction    (NFunction, ($*))

import Data.Either

type Model    n r = NFunction n r [r]

-- | See: <http://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant>
type Jacobian n r = NFunction n r [SizedList n r]

levmar :: forall n k r. (Nat n, Nat k, LMA_I.LevMarable r)
       => (Model n r)                     -- ^ Model
       -> Maybe (Jacobian n r)            -- ^ Optional jacobian
       -> SizedList n r                   -- ^ Initial parameters
       -> [r]                             -- ^ Samples
       -> Integer                         -- ^ Maximum number of iterations
       -> LMA_I.Options r                 -- ^ Options
       -> Maybe (SizedList n r)           -- ^ Optional lower bounds
       -> Maybe (SizedList n r)           -- ^ Optional upper bounds
       -> Maybe (LinearConstraints k n r) -- ^ Optional linear constraints
       -> Maybe (SizedList n r)           -- ^ Optional weights
       -> Either LMA_I.LevMarError (SizedList n r, LMA_I.Info r, CovarMatrix n r)

levmar model mJac params ys itMax opts mLowBs mUpBs mLinC mWghts =
    fmap convertResult $ LMA_I.levmar (convertModel model)
                                      (fmap convertJacob mJac)
                                      (toList params)
                                      ys
                                      itMax
                                      opts
                                      (fmap toList mLowBs)
                                      (fmap toList mUpBs)
                                      (fmap convertLinC mLinC)
                                      (fmap toList mWghts)
    where
      convertModel f = \ps ->              f $* (unsafeFromList ps :: SizedList n r)
      convertJacob f = \ps -> map toList ((f $* (unsafeFromList ps :: SizedList n r)) :: [SizedList n r])
      convertLinC (cMat, rhcVec) = ( map toList $ toList cMat
                                   , toList rhcVec
                                   )
      convertResult (psResult, info, covar) = ( unsafeFromList psResult
                                              , info
                                              , unsafeFromList $ map unsafeFromList covar
                                              )

type CovarMatrix n r = Matrix n n r

type LinearConstraints k n r = (Matrix k n r, SizedList k r)

type Matrix n m r = SizedList n (SizedList m r)

-- |Value to denote the absense of any linear constraints over the
-- parameters of the model function. Use this instead of 'Nothing'
-- because the type parameter which contains the number of constraints
-- can't be inferred.
noLinearConstraints :: Nat n => Maybe (LinearConstraints Z n r)
noLinearConstraints = Nothing
