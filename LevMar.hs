{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ScopedTypeVariables #-}

module LevMar
    ( levmar
    , Model
    , Jacobian
    , CovarMatrix
    , LecMatrix

    , levmar'
    , Model'
    , Jacobian'

    , noLinearConstraints

    , Z, S, Nat
    , SizedList(..)
    , NFunction

    , LMA_I.Options(..)
    , LMA_I.defaultOpts
    , LMA_I.StopReason(..)
    , LMA_I.Info(..)
    )
    where

import qualified LevMar.Intermediate as LMA_I

import TypeLevelNat (Z, S, Nat, witnessNat)
import SizedList    (SizedList(..), toList, unsafeFromList)
import NFunction    (NFunction, ($*), ComposeN, compose)

levmar :: forall n k r a. (Nat n, ComposeN n, Nat k, LMA_I.LevMarable r)
       => (Model n r a)                          -- model
       -> Maybe (Jacobian n r a)                 -- jacobian
       -> SizedList n r                          -- init params
       -> [(a,r)]                                -- samples
       -> Integer                                -- max iterations
       -> LMA_I.Options r                        -- options
       -> Maybe (SizedList n r)                  -- lower bounds
       -> Maybe (SizedList n r)                  -- upper bounds
       -> Maybe (LecMatrix k n r, SizedList k r) -- linear constraints
       -> Maybe (SizedList n r)                  -- weights
       -> Maybe (SizedList n r, LMA_I.Info r, CovarMatrix n r)
levmar model mJac params samples = levmar' (convertModel model)
                                           (fmap convertJacob mJac)
                                           params
                                           ys
    where
      (xs, ys) = unzip samples

      convertModel :: Model n r a -> Model' n r
      convertModel = compose (witnessNat :: n) (undefined :: r)
                             (\(f :: a -> r) -> map f xs)

      convertJacob :: Jacobian n r a -> Jacobian' n r
      convertJacob = compose (witnessNat :: n) (undefined :: r)
                             (\(f :: a -> SizedList n r) -> map f xs)

type Model    n r a = NFunction n r (a -> r)
type Jacobian n r a = NFunction n r (a -> SizedList n r)

type CovarMatrix n r = Matrix n n r
type LecMatrix k n r = Matrix k n r

type Matrix n m r = SizedList n (SizedList m r)

type Model'    n r = NFunction n r [r]
type Jacobian' n r = NFunction n r [SizedList n r]

-- |Value to denote the absense of any linear constraints over the
-- parameters of the model function. This is necessary because the
-- type parameter which contains the number of constraints can't be
-- inferred.
noLinearConstraints :: Nat n => Maybe (LecMatrix Z n r, SizedList Z r)
noLinearConstraints = Nothing

levmar' :: forall n k r. (Nat n, Nat k, LMA_I.LevMarable r)
        => (Model' n r)                            -- model
        -> Maybe (Jacobian' n r)                   -- jacobian
        -> SizedList n r                          -- init params
        -> [r]                                    -- samples
        -> Integer                                -- max iterations
        -> LMA_I.Options r                        -- options
        -> Maybe (SizedList n r)                  -- lower bounds
        -> Maybe (SizedList n r)                  -- upper bounds
        -> Maybe (LecMatrix k n r, SizedList k r) -- linear constraints
        -> Maybe (SizedList n r)                  -- weights
        -> Maybe (SizedList n r, LMA_I.Info r, CovarMatrix n r)

levmar' model mJac params ys itMax opts mLowBs mUpBs mLinC mWghts =
    fmap convertResult $ LMA_I.levmar' (convertModel model)
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
