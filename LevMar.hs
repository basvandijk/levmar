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

type Model a r n = NFunction n r (a -> r)

type Jacobian a r n = NFunction n r (a -> SizedList r n)

type CovarMatrix r n = SizedList (SizedList r n) n

type LevMar a r n =  (Model a r n)          -- model
                  -> Maybe (Jacobian a r n) -- jacobian
                  -> SizedList r n          -- init params
                  -> [(a, r)]               -- samples
                  -> Integer                -- max iterations
                  -> LMA_I.Options r
                  -> Maybe (SizedList r n)  -- lower bounds
                  -> Maybe (SizedList r n)  -- upper bounds
                  -> (SizedList r n, LMA_I.Info r, CovarMatrix r n)

levmar :: forall n r a. (Nat n, LMA_I.LevMarable r) => LevMar a r n
levmar model mJac params samples itMax opts mLowBs mUpBs =
    let mkModel f = \ps x ->           (f $* (unsafeFromList ps :: SizedList r n)) x
        mkJacob f = \ps x -> toList $ ((f $* (unsafeFromList ps :: SizedList r n)) x :: SizedList r n)
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
