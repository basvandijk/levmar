{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ScopedTypeVariables #-}

module LevMar
    ( levmar

    , LevMar
    , NFunction
    , Model
    , Jacobian
    , CovarMatrix

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

type Model n r a = NFunction n r (a -> r)

type Jacobian n r a = NFunction n r (a -> SizedList n r)

type CovarMatrix n r = SizedList n (SizedList n r)

type LevMar n r a =  (Model n r a)          -- model
                  -> Maybe (Jacobian n r a) -- jacobian
                  -> SizedList n r          -- init params
                  -> [(a, r)]               -- samples
                  -> Integer                -- max iterations
                  -> LMA_I.Options r
                  -> Maybe (SizedList n r)  -- lower bounds
                  -> Maybe (SizedList n r)  -- upper bounds
                  -> (SizedList n r, LMA_I.Info r, CovarMatrix n r)

levmar :: forall n r a. (Nat n, LMA_I.LevMarable r) => LevMar n r a
levmar model mJac params samples itMax opts mLowBs mUpBs =
    let mkModel f = \ps x ->           (f $* (unsafeFromList ps :: SizedList n r)) x
        mkJacob f = \ps x -> toList $ ((f $* (unsafeFromList ps :: SizedList n r)) x :: SizedList n r)
        (psResult, info, covar) = LMA_I.levmar (mkModel model)
                                               (fmap mkJacob mJac)
                                               (toList params)
                                               samples
                                               itMax
                                               opts
                                               (fmap toList mLowBs)
                                               (fmap toList mUpBs)
    in ( unsafeFromList psResult
       , info
       , unsafeFromList $ map unsafeFromList covar
       )
