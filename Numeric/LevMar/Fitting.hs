{-# LANGUAGE NoImplicitPrelude, UnicodeSyntax #-}

--------------------------------------------------------------------------------
-- |
-- Module:     Numeric.LevMar.Fitting
-- Copyright:  (c) 2009 - 2010 Roel van Dijk & Bas van Dijk
-- License:    BSD-style (see the file LICENSE)
-- Maintainer: Roel van Dijk <vandijk.roel@gmail.com>
--             Bas van Dijk <v.dijk.bas@gmail.com>
-- Stability:  Experimental
--
-- This module provides the Levenberg-Marquardt algorithm specialised
-- for curve-fitting.
--
-- For additional documentation see the documentation of the levmar C
-- library which this library is based on:
-- <http://www.ics.forth.gr/~lourakis/levmar/>
--
--------------------------------------------------------------------------------

module Numeric.LevMar.Fitting
    ( -- * Model & Jacobian.
      Model
    , SimpleModel
    , Jacobian
    , SimpleJacobian

      -- * Levenberg-Marquardt algorithm.
    , LevMar.LevMarable
    , levmar

    , LevMar.LinearConstraints

      -- * Minimization options.
    , LevMar.Options(..)
    , LevMar.defaultOpts

      -- * Output
    , LevMar.Info(..)
    , LevMar.StopReason(..)
    , LevMar.CovarMatrix

    , LevMar.LevMarError(..)
    ) where


--------------------------------------------------------------------------------
-- Imports
--------------------------------------------------------------------------------

-- from base:
import Data.Functor  ( fmap )
import Data.Either   ( Either )
import Data.List     ( map, unzip )
import Data.Maybe    ( Maybe )
import Prelude       ( Integer )

-- from levmar:
import qualified Numeric.LevMar as LevMar


--------------------------------------------------------------------------------
-- Model & Jacobian.
--------------------------------------------------------------------------------

{-| A functional relation describing measurements represented as a function
from a list of parameters and an x-value to an expected measurement.

 * Ensure that the length of the parameters list equals the lenght of the initial
   parameters list in 'levmar'.

For example, the quadratic function @f(x) = a*x^2 + b*x + c@ can be
written as:

@
quad :: 'Num' r => 'Model' r r
quad [a, b, c] x = a*x^2 + b*x + c
@
-}
type Model r a = [r] → (a → r)

-- | This type synonym expresses that usually the @a@ in @'Model' r a@
-- equals the type of the parameters.
type SimpleModel r = Model r r

{-| The jacobian of the 'Model' function. Expressed as a function from a list
of parameters and an x-value to the partial derivatives of the parameters.

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

(Notice you don't have to differentiate for @x@.)
-}
type Jacobian r a = [r] → (a → [r])

-- | This type synonym expresses that usually the @a@ in @'Jacobian' r a@
-- equals the type of the parameters.
type SimpleJacobian r = Jacobian r r


--------------------------------------------------------------------------------
-- Levenberg-Marquardt algorithm.
--------------------------------------------------------------------------------

-- | The Levenberg-Marquardt algorithm specialised for curve-fitting.
levmar ∷ LevMar.LevMarable r
       ⇒ Model r a                          -- ^ Model
       → Maybe (Jacobian r a)               -- ^ Optional jacobian
       → [r]                                -- ^ Initial parameters
       → [(a, r)]                           -- ^ Samples
       → Integer                            -- ^ Maximum iterations
       → LevMar.Options r                   -- ^ Minimization options
       → LevMar.Constraints r               -- ^ Constraints
       → Either LevMar.LevMarError ([r], LevMar.Info r, LevMar.CovarMatrix r)
levmar model mJac params samples =
    LevMar.levmar (convertModel model)
                  (fmap convertJacob mJac)
                  params
                  ys
        where
          (xs, ys) = unzip samples

          convertModel mdl = \ps → map (mdl ps) xs
          convertJacob jac = \ps → map (jac ps) xs


-- The End ---------------------------------------------------------------------
