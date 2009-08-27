{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ScopedTypeVariables #-}

module LevMar
    ( levmar

    , LevMar
    , ParamFunc
    , Model
    , Jacobian
    , CovarMatrix

    , Z, S, Nat
    , SizedList(..)

    , LMA_I.Options(..)
    , LMA_I.defaultOpts
    , LMA_I.StopReason(..)
    , LMA_I.Info(..)

    -- TODO: The following is only needed for the test:
    , ($*)
    )
    where

import qualified LevMar.Intermediate as LMA_I

import TypeLevelNat (Z, S, Nat)
import SizedList (SizedList(..), toList, unsafeFromList)

-- | @ParamFunc r n@ represents a function from @n@ @r@'s to a @r@.
-- For example: @ParamFunc Double (S (S (S Z))) ~ Double -> Double -> Double -> Double@
type family ParamFunc r n :: *

type instance ParamFunc r Z     = r
type instance ParamFunc r (S n) = r -> ParamFunc r n

-- | @f $* ps@ applies the /n/-arity function @f@ to each of the parameters in
-- the /n/-sized list @ps@.
($*) :: ParamFunc r n -> SizedList r n -> r
f $* Nil        = f
f $* (p ::: ps) = f p $* ps

infixr 0 $* -- same as $

type Model a r n = a -> ParamFunc r n

type Jacobian a r n = Model a r n

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
    let mkModel f = \x ps -> f x $* (unsafeFromList ps :: SizedList r n)
        (psResult, info, covar) = LMA_I.levmar (mkModel model)
                                               (fmap mkModel mJac)
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
