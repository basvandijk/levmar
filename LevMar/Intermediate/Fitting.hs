--------------------------------------------------------------------------------
-- |
-- Module      :  LevMar.Intermediate.Fitting
-- Copyright   :  (c) 2009 Roel van Dijk & Bas van Dijk
-- License     :  BSD-style (see the file LICENSE)
--
-- Maintainer  :  vandijk.roel@gmail.com, v.dijk.bas@gmail.com
-- Stability   :  Experimental
--
-- This module provides the Levenberg-Marquardt algorithm specialised
-- for curve-fitting.
--
-- For additional documentation see the documentation of the levmar C
-- library which this library is based on:
-- <http://www.ics.forth.gr/~lourakis/levmar/>
--
--------------------------------------------------------------------------------

module LevMar.Intermediate.Fitting
    ( -- * Model & Jacobian.
      Model
    , Jacobian

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


import qualified LevMar.Intermediate as LMA_I


--------------------------------------------------------------------------------
-- Model & Jacobian.
--------------------------------------------------------------------------------

type Model r a = [r] -> a -> r

-- | See: <http://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant>
type Jacobian r a = [r] -> a -> [r]


--------------------------------------------------------------------------------
-- Levenberg-Marquardt algorithm.
--------------------------------------------------------------------------------

-- | The Levenberg-Marquardt algorithm specialised for curve-fitting.
levmar :: LMA_I.LevMarable r
       => Model r a                         -- ^ Model
       -> Maybe (Jacobian r a)              -- ^ Optional jacobian
       -> [r]                               -- ^ Initial parameters
       -> [(a, r)]                          -- ^ Samples
       -> Integer                           -- ^ Maximum iterations
       -> LMA_I.Options r                   -- ^ Minimization options
       -> Maybe [r]                         -- ^ Optional lower bounds
       -> Maybe [r]                         -- ^ Optional upper bounds
       -> Maybe (LMA_I.LinearConstraints r) -- ^ Optional linear constraints
       -> Maybe [r]                         -- ^ Optional weights
       -> Either LMA_I.LevMarError ([r], LMA_I.Info r, LMA_I.CovarMatrix r)
levmar model mJac ps samples =
    LMA_I.levmar (\ps' -> map (model ps') xs)
                 (fmap (\jac -> \ps' -> map (jac ps') xs) mJac)
                 ps
                 ys
        where
          (xs, ys) = unzip samples


-- The End ---------------------------------------------------------------------
