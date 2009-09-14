{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleContexts #-}

--------------------------------------------------------------------------------
-- |
-- Module      :  LevMar.Intermediate.Fitting.AD
-- Copyright   :  (c) 2009 Roel van Dijk & Bas van Dijk
-- License     :  BSD-style (see the file LICENSE)
--
-- Maintainer  :  vandijk.roel@gmail.com, v.dijk.bas@gmail.com
-- Stability   :  Experimental
--
-- A levmar variant specialised for curve-fitting that uses Automatic
-- Differentiation to automatically compute the Jacobian.
--
-- For additional documentation see the documentation of the levmar C
-- library which this library is based on:
-- <http://www.ics.forth.gr/~lourakis/levmar/>
--
--------------------------------------------------------------------------------

module LevMar.Intermediate.Fitting.AD
    ( -- * Model.
      Model
    , SimpleModel

      -- * Levenberg-Marquardt algorithm.
    , LMA_I.LevMarable
    , levmar

    , LMA_I.LinearConstraints

      -- * Minimization options.
    , LMA_I.Options(..)
    , LMA_I.defaultOpts

      -- * Output
    , LMA_I.Info(..)
    , LMA_I.StopReason(..)
    , LMA_I.CovarMatrix

    , LMA_I.LevMarError(..)
    ) where


import qualified LevMar.Intermediate.Fitting as LMA_I

import LevMar.Utils.AD  ( firstDeriv, constant, idDAt )

-- From vector-space:
import Data.Derivative  ( (:~>), (:>), powVal )
import Data.VectorSpace ( VectorSpace, Scalar )
import Data.Basis       ( HasBasis, Basis )


--------------------------------------------------------------------------------
-- Model
--------------------------------------------------------------------------------

{- | A functional relation describing measurements.

 * Ensure that the length of the parameters list equals the length of the initial
   parameters list in 'levmar'.

For example:

@
quad :: 'Num' r => 'Model' r r
quad [a, b, c] x = a*x^2 + b*x + c
@
-}
type Model r a = [r :~> r] -> a -> (r :~> r)

-- | This type synonym expresses that usually the @a@ in @'Model' r a@
-- equals the type of the parameters.
type SimpleModel r = Model r (r :~> r)


--------------------------------------------------------------------------------
-- Levenberg-Marquardt algorithm.
--------------------------------------------------------------------------------

-- | The Levenberg-Marquardt algorithm specialised for curve-fitting
-- that automatically computes the 'Jacobian' using automatic
-- differentiation of the model function.
levmar :: forall r a.
          ( HasBasis r
          , Basis r ~ ()
          , VectorSpace (Scalar r)
          , LMA_I.LevMarable r
          )
       => Model r a                         -- ^ Model
       -> [r]                               -- ^ Initial parameters
       -> [(a, r)]                          -- ^ Samples
       -> Integer                           -- ^ Maximum iterations
       -> LMA_I.Options r                   -- ^ Minimization options
       -> Maybe [r]                         -- ^ Optional lower bounds
       -> Maybe [r]                         -- ^ Optional upper bounds
       -> Maybe (LMA_I.LinearConstraints r) -- ^ Optional linear constraints
       -> Maybe [r]                         -- ^ Optional weights
       -> Either LMA_I.LevMarError ([r], LMA_I.Info r, LMA_I.CovarMatrix r)

levmar model = LMA_I.levmar (convertModel model) $ Just $ jacobianOf model
    where
      convertModel :: Model r a -> LMA_I.Model r a
      convertModel f = \ps x -> powVal $ f (map constant ps) x undefined

      jacobianOf :: Model r a -> LMA_I.Jacobian r a
      jacobianOf f =
          \ps x -> map (\(ix, p) -> firstDeriv $ f (idDAt ix ps) x p) $ zip [0..] ps


-- The End ---------------------------------------------------------------------
