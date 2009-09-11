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

{- | A functional relation describing measurements represented as a function
from a list of parameters of type @r@ and an x-value of type @a@ to a value of
type @r@.

 * Ensure that the length of the parameters list equals the lenght of the initial
   parameters list in 'levmar'.

For example, the quadratic function @f(x) = a*x^2 + b*x + c@ can be
written as:

@
quad :: 'Num' r => 'Model' r r
quad [a, b, c] x = a*x^2 + b*x + c
@
-}
type Model r a = [r] -> a -> r

{- | The jacobian of the 'Model' function. Expressed as a function from a list
of parameters of type @r@ and an x-value of type @a@ to a vector of @n@ values
of type @r@.

See: <http://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant>

 * Ensure that the length of the parameters list equals the lenght of the initial
   parameters list in 'levmar'.

 * Ensure that the length of the output parameter derivatives list equals the
   length of the input parameters list.

For example, the jacobian of the above @quad@ model can be written as:

@
quadJacob :: 'Num' r => 'Jacobian' N3 r r
quadJacob [_, _, _] x = [ x^2   -- with respect to a
                        , x     -- with respect to b
                        , 1     -- with respect to c
                        ]
@

Notice you don't have to differentiate for @x@.
-}
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
