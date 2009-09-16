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
      LMA_I.Model
    , LMA_I.SimpleModel
    , LMA_I.Jacobian
    , LMA_I.SimpleJacobian
    , jacobianOf

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
levmar :: forall r a.
          ( HasBasis r
          , Basis r ~ ()
          , VectorSpace (Scalar r)
          , LMA_I.LevMarable r
          )
       => LMA_I.Model (r :~> r) a           -- ^ Model. Note that
                                            --   ':~>' is overloaded
                                            --   for all the numeric
                                            --   classes.
       -> [r]                               -- ^ Initial parameters
       -> [(a, r)]                          -- ^ Samples
       -> Integer                           -- ^ Maximum iterations
       -> LMA_I.Options r                   -- ^ Minimization options
       -> Maybe [r]                         -- ^ Optional lower bounds
       -> Maybe [r]                         -- ^ Optional upper bounds
       -> Maybe (LMA_I.LinearConstraints r) -- ^ Optional linear constraints
       -> Maybe [r]                         -- ^ Optional weights
       -> Either LMA_I.LevMarError ([r], LMA_I.Info r, LMA_I.CovarMatrix r)

levmar model = LMA_I.levmar (convertModel model) . Just $ jacobianOf model
    where
      convertModel :: LMA_I.Model (r :~> r) a -> LMA_I.Model r a
      convertModel mdl = \ps -> value . mdl (map constant ps)

-- | Compute the 'LMA_I.Jacobian' of the 'LMA_I.Model' using Automatic
-- Differentiation.
jacobianOf :: (HasBasis r, Basis r ~ (), VectorSpace (Scalar r))
           => LMA_I.Model (r :~> r) a -> LMA_I.Jacobian r a
jacobianOf mdl =
    \ps x -> map (\(ix, p) -> firstDeriv $ mdl (idDAt ix ps) x p) $
                 zip [0..] ps


-- The End ---------------------------------------------------------------------
