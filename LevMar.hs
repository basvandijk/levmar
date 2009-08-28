{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ScopedTypeVariables #-}

module LevMar
    ( levmar

    , LevMar
    , NFunction
    , Model
    , Jacobian
    , CovarMatrix
    , LecMatrix

    , noLinearConstraints

    , Z, S, Nat
    , SizedList(..)

    , LMA_I.Options(..)
    , LMA_I.defaultOpts
    , LMA_I.StopReason(..)
    , LMA_I.Info(..)
    )
    where

import qualified LevMar.Intermediate as LMA_I

import TypeLevelNat (Z, S, Nat)
import SizedList    (SizedList(..), toList, unsafeFromList)
import NFunction    (NFunction, ($*))

type Model    n r a = NFunction n r (a -> r)
type Jacobian n r a = NFunction n r (a -> SizedList n r)
type Matrix n m r = SizedList n (SizedList m r)
type CovarMatrix n r = Matrix n n r
type LecMatrix k n r = Matrix k n r

-- |Value to denote the absense of any linear constraints over the
-- parameters of the model function. This is necessary because the
-- type parameter which contains the number of constraints can't be
-- inferred.
noLinearConstraints :: Nat n => Maybe (LecMatrix Z n r, SizedList Z r)
noLinearConstraints = Nothing

type LevMar n k r a =  (Model n r a)          -- model
                    -> Maybe (Jacobian n r a) -- jacobian
                    -> SizedList n r          -- init params
                    -> [(a, r)]               -- samples
                    -> Integer                -- max iterations
                    -> LMA_I.Options r
                    -> Maybe (SizedList n r)  -- lower bounds
                    -> Maybe (SizedList n r)  -- upper bounds
                    -> Maybe (LecMatrix k n r, SizedList k r) -- linear constraints
                    -> Maybe (SizedList n r) -- weights
                    -> Maybe (SizedList n r, LMA_I.Info r, CovarMatrix n r)

levmar :: forall n k r a. (Nat n, Nat k, LMA_I.LevMarable r) => LevMar n k r a
levmar model mJac params samples itMax opts mLowBs mUpBs mLinC mWghts =
    fmap convertResult $ LMA_I.levmar (convertModel model)
                                      (fmap convertJacob mJac)
                                      (toList params)
                                      samples
                                      itMax
                                      opts
                                      (fmap toList mLowBs)
                                      (fmap toList mUpBs)
                                      (fmap convertLinC mLinC)
                                      (fmap toList mWghts)
    where
      convertModel f = \ps x ->         (f $* (unsafeFromList ps :: SizedList n r)) x
      convertJacob f = \ps x -> toList ((f $* (unsafeFromList ps :: SizedList n r)) x :: SizedList n r)
      convertLinC (cMat, rhcVec) = ( map toList $ toList cMat
                                   , toList rhcVec
                                   )
      convertResult (psResult, info, covar) = ( unsafeFromList psResult
                                              , info
                                              , unsafeFromList $ map unsafeFromList covar
                                              )
