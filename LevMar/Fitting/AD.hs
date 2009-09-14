{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleContexts #-}

--------------------------------------------------------------------------------
-- |
-- Module      :  LevMar.Fitting.AD
-- Copyright   :  (c) 2009 Roel van Dijk & Bas van Dijk
-- License     :  BSD-style (see the file LICENSE)
--
-- Maintainer  :  vandijk.roel@gmail.com, v.dijk.bas@gmail.com
-- Stability   :  Experimental
--
-- This module provides the Levenberg-Marquardt algorithm specialised
-- for curve-fitting that uses Automatic Differentiation to
-- automatically compute the Jacobian.
--
-- For additional documentation see the documentation of the levmar C
-- library which this library is based on:
-- <http://www.ics.forth.gr/~lourakis/levmar/>
--
--------------------------------------------------------------------------------

module LevMar.Fitting.AD
    ( -- * Model.
      Model
    , SimpleModel

      -- * Levenberg-Marquardt algorithm.
    , LMA_I.LevMarable
    , levmar

    , LinearConstraints
    , noLinearConstraints
    , Matrix

    -- * Minimization options.
    , LMA_I.Options(..)
    , LMA_I.defaultOpts

      -- * Output
    , LMA_I.Info(..)
    , LMA_I.StopReason(..)
    , CovarMatrix

    , LMA_I.LevMarError(..)

      -- *Type-level machinery
    , Z, S, Nat
    , SizedList(..)
    , NFunction
    ) where


import qualified LevMar.Intermediate.Fitting as LMA_I

import LevMar.Utils ( LinearConstraints
                    , noLinearConstraints
                    , convertLinearConstraints
                    , Matrix
                    , CovarMatrix
                    , convertResult
                    )

import TypeLevelNat ( Z, S, Nat )
import SizedList    ( SizedList(..), toList, unsafeFromList )
import NFunction    ( NFunction, ($*) )

import LevMar.Utils.AD  ( firstDeriv, constant, idDAt )

-- From vector-space:
import Data.Derivative  ( (:~>), (:>), powVal )
import Data.VectorSpace ( VectorSpace, Scalar, AdditiveGroup )
import Data.Basis       ( HasBasis, Basis )
import Data.MemoTrie    ( HasTrie )


--------------------------------------------------------------------------------
-- Model & Jacobian.
--------------------------------------------------------------------------------

{- | A functional relation describing measurements represented as a function
from @n@ parameters of type @r@ and an x-value of type @a@ to a value of type @r@.

For example, the quadratic function @f(x) = a*x^2 + b*x + c@ can be
written as:

@
type N3 = 'S' ('S' ('S' 'Z'))

quad :: 'Num' r => 'Model' N3 r (r :~> r)
quad a b c x = a*x^2 + b*x + c
@
-}
type Model n r a = NFunction n (r :~> r) (a -> (r :~> r))

-- | This type synonym expresses that usually the @a@ in @'Model' n r a@
-- equals the type of the parameters.
type SimpleModel n r = Model n r (r :~> r)


--------------------------------------------------------------------------------
-- Levenberg-Marquardt algorithm.
--------------------------------------------------------------------------------

-- | The Levenberg-Marquardt algorithm specialised for curve-fitting
-- that automatically computes the 'Jacobian' using automatic
-- differentiation of the model function.
levmar :: forall n k r a.
          ( Nat n
          , Nat k
          , HasBasis r
          , Basis r ~ ()
          , VectorSpace (Scalar r)
          , LMA_I.LevMarable r
          )
       => (Model n r a)                       -- ^ Model
       -> SizedList n r                       -- ^ Initial parameters
       -> [(a, r)]                            -- ^ Samples
       -> Integer                             -- ^ Maximum number of iterations
       -> LMA_I.Options r                       -- ^ Minimization options
       -> Maybe (SizedList n r)               -- ^ Optional lower bounds
       -> Maybe (SizedList n r)               -- ^ Optional upper bounds
       -> Maybe (LinearConstraints k n r)     -- ^ Optional linear constraints
       -> Maybe (SizedList n r)               -- ^ Optional weights
       -> Either LMA_I.LevMarError (SizedList n r, LMA_I.Info r, CovarMatrix n r)

levmar model params ys itMax opts mLowBs mUpBs mLinC mWghts =
    fmap convertResult $ LMA_I.levmar (convertModel model)
                                      (Just $ jacobianOf model)
                                      (toList params)
                                      ys
                                      itMax
                                      opts
                                      (fmap toList mLowBs)
                                      (fmap toList mUpBs)
                                      (fmap convertLinearConstraints mLinC)
                                      (fmap toList mWghts)
    where
      convertModel :: (HasBasis r, HasTrie (Basis r))
                   => Model n r a -> LMA_I.Model r a
      (convertModel f) ps x =
          let pDs = unsafeFromList (fmap constant ps) :: SizedList n (r :~> r)
          in powVal $ (f $* pDs :: a -> r :~> r) x undefined

      jacobianOf :: ( HasBasis r
                    , Basis r ~ ()
                    , HasTrie (Basis r)
                    , AdditiveGroup r
                    , VectorSpace (Scalar r)
                    ) => Model n r a -> LMA_I.Jacobian r a
      (jacobianOf f) ps x = fmap combine $ zip [0..] ps
          where
            combine (ix, p) = let pDs = unsafeFromList (idDAt ix ps) :: SizedList n (r :~> r)
                              in firstDeriv $ (f $* pDs :: a -> r :~> r) x p


-- The End ---------------------------------------------------------------------
