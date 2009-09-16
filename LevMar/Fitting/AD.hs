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
      LMA.Model
    , LMA.SimpleModel

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


import qualified LevMar.Fitting              as LMA
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

import LevMar.Utils.AD  ( value, firstDeriv, constant, idDAt )

-- From vector-space:
import Data.Derivative  ( (:~>) )
import Data.VectorSpace ( VectorSpace, Scalar )
import Data.Basis       ( HasBasis, Basis )


--------------------------------------------------------------------------------
-- Levenberg-Marquardt algorithm.
--------------------------------------------------------------------------------

-- | The Levenberg-Marquardt algorithm specialised for curve-fitting
-- that automatically computes the 'Jacobian' using automatic
-- differentiation of the model function.
--
-- /Warning/: Don't apply 'levmar' to 'LMA_I.Model's that apply methods of
-- the 'Eq' and 'Ord' classes to the parameters. These methods are
-- undefined for ':~>'!!!
levmar :: forall m k r a.
          ( Nat m
          , Nat k
          , HasBasis r
          , Basis r ~ ()
          , VectorSpace (Scalar r)
          , LMA_I.LevMarable r
          )
       => LMA.Model m (r :~> r) a             -- ^ Model. Note that
                                              --   ':~>' is overloaded
                                              --   for all the numeric
                                              --   classes.
       -> SizedList m r                       -- ^ Initial parameters
       -> [(a, r)]                            -- ^ Samples
       -> Integer                             -- ^ Maximum number of iterations
       -> LMA_I.Options r                       -- ^ Minimization options
       -> Maybe (SizedList m r)               -- ^ Optional lower bounds
       -> Maybe (SizedList m r)               -- ^ Optional upper bounds
       -> Maybe (LinearConstraints k m r)     -- ^ Optional linear constraints
       -> Maybe (SizedList m r)               -- ^ Optional weights
       -> Either LMA_I.LevMarError (SizedList m r, LMA_I.Info r, CovarMatrix m r)

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
      convertModel :: LMA.Model m (r :~> r) a -> LMA_I.Model r a
      (convertModel f) ps x = value $ (f $* pDs :: a -> r :~> r) x
          where
            pDs :: SizedList m (r :~> r)
            pDs = unsafeFromList $ fmap constant ps

      jacobianOf :: LMA.Model m (r :~> r) a -> LMA_I.Jacobian r a
      (jacobianOf f) ps x = fmap combine $ zip [0..] ps
          where
            combine (ix, p) = firstDeriv $ (f $* pDs :: a -> r :~> r) x p
                where
                  pDs :: SizedList m (r :~> r)
                  pDs = unsafeFromList $ idDAt ix ps


-- The End ---------------------------------------------------------------------
