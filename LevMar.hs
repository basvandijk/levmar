{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ScopedTypeVariables #-}

module LevMar
    ( -- *LevMar with a per-sample model
      Model
    , Jacobian
    , levmar

      -- *LevMar with a model for all samples
    , Model'
    , Jacobian'
    , levmar'

      -- *Supporting types and default values
    , LMA_I.Options(..)
    , LMA_I.defaultOpts
    , LMA_I.StopReason(..)
    , LMA_I.Info(..)
    , noLinearConstraints
    , CovarMatrix
    , LinearConstraints
    , LecMatrix

      -- *Type-level stuff
    , Z, S, Nat
    , SizedList(..)
    , NFunction
    , ComposeN
    )
    where

import qualified LevMar.Intermediate as LMA_I

import TypeLevelNat (Z, S, Nat, witnessNat)
import SizedList    (SizedList(..), toList, unsafeFromList)
import NFunction    (NFunction, ($*), ComposeN, compose)

-------------------------------------------------------------------------------

-- |A function from @n@ parameters of type @r@ and an x-value of type
-- @a@ to a value of type @r@.
--
-- For example, the quadratic function @f(x) = a*x^2 + b*x + c@ can be
-- written as:
--
-- @
--   type N3 = 'S' ('S' ('S' 'Z'))
--   quad :: 'Num' r => 'Model' N3 r r
--   quad a b c x = a*x^2 + b*x + c
-- @
type Model n r a = NFunction n r (a -> r)

-- |The jacobian of the 'Model' function. Expressed as a function from
-- @n@ parameters of type @r@ and an x-value of type @a@ to a vector
-- of @n@ values of type @r@.
--
-- For example, the jacobian of the quadratic function @f(x) = a*x^2 +
-- b*x + c@ can be written as:
--
-- @
--   type N3 = 'S' ('S' ('S' 'Z'))
--   quadJacob :: 'Num' r => 'Jacobian' N3 r r
--   quadJacob _ _ _ x =   x^2   -- with respect to a
--                     ::: x     -- with respect to b
--                     ::: 1     -- with respect to c
--                     ::: 'Nil'
-- @
--
-- Notice you don't have to differentiate for @x@.
type Jacobian n r a = NFunction n r (a -> SizedList n r)

levmar :: forall n k r a. (Nat n, ComposeN n, Nat k, LMA_I.LevMarable r)
       => (Model n r a)                          -- ^Model
       -> Maybe (Jacobian n r a)                 -- ^Jacobian
       -> SizedList n r                          -- ^Initial parameters
       -> [(a,r)]                                -- ^Samples
       -> Integer                                -- ^Maximum number of iterations
       -> LMA_I.Options r                        -- ^Options
       -> Maybe (SizedList n r)                  -- ^Lower bounds
       -> Maybe (SizedList n r)                  -- ^Upper bounds
       -> Maybe (LinearConstraints k n r) -- ^Linear constraints
       -> Maybe (SizedList n r)                  -- ^Weights
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

-------------------------------------------------------------------------------

type Model'    n r = NFunction n r [r]
type Jacobian' n r = NFunction n r [SizedList n r]

levmar' :: forall n k r. (Nat n, Nat k, LMA_I.LevMarable r)
        => (Model' n r)                           -- ^Model
        -> Maybe (Jacobian' n r)                  -- ^Jacobian
        -> SizedList n r                          -- ^Initial parameters
        -> [r]                                    -- ^Samples
        -> Integer                                -- ^Maximum number of iterations
        -> LMA_I.Options r                        -- ^Options
        -> Maybe (SizedList n r)                  -- ^Lower bounds
        -> Maybe (SizedList n r)                  -- ^Upper bounds
        -> Maybe (LecMatrix k n r, SizedList k r) -- ^Linear constraints
        -> Maybe (SizedList n r)                  -- ^Weights
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

-------------------------------------------------------------------------------

type CovarMatrix n r = Matrix n n r

type LinearConstraints k n r = (LecMatrix k n r, SizedList k r)
type LecMatrix         k n r = Matrix k n r

type Matrix n m r = SizedList n (SizedList m r)

-- |Value to denote the absense of any linear constraints over the
-- parameters of the model function. This is necessary because the
-- type parameter which contains the number of constraints can't be
-- inferred.
noLinearConstraints :: Nat n => Maybe (LecMatrix Z n r, SizedList Z r)
noLinearConstraints = Nothing

