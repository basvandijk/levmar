{-# LANGUAGE ScopedTypeVariables #-}

--------------------------------------------------------------------------------
-- |
-- Module      :  LevMar.Fitting
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

module LevMar.Fitting
    ( -- * Model & Jacobian.
      Model
    , SimpleModel
    , Jacobian
    , SimpleJacobian

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


--------------------------------------------------------------------------------
-- Model & Jacobian.
--------------------------------------------------------------------------------

{- | A functional relation describing measurements represented as a function
from @n@ parameters of type @r@ and an x-value of type @a@ to a value of type @r@.

For example, the quadratic function @f(x) = a*x^2 + b*x + c@ can be
written as:

@
type N3 = 'S' ('S' ('S' 'Z'))

quad :: 'Num' r => 'Model' N3 r r
quad a b c x = a*x^2 + b*x + c
@
-}
type Model n r a = NFunction n r (a -> r)

-- | This type synonym expresses that usually the @a@ in @'Model' n r a@
-- equals the type of the parameters.
type SimpleModel n r = Model n r r

{- | The jacobian of the 'Model' function. Expressed as a function from
@n@ parameters of type @r@ and an x-value of type @a@ to a vector
of @n@ values of type @r@.

See: <http://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant>

For example, the jacobian of the quadratic function @f(x) = a*x^2 +
b*x + c@ can be written as:

@
type N3 = 'S' ('S' ('S' 'Z'))

quadJacob :: 'Num' r => 'Jacobian' N3 r r
quadJacob _ _ _ x =   x^2   -- with respect to a
                  ::: x     -- with respect to b
                  ::: 1     -- with respect to c
                  ::: 'Nil'
@

Notice you don't have to differentiate for @x@.
-}
type Jacobian n r a = NFunction n r (a -> SizedList n r)

-- | This type synonym expresses that usually the @a@ in @'Jacobian' n r a@
-- equals the type of the parameters.
type SimpleJacobian n r = Jacobian n r r


--------------------------------------------------------------------------------
-- Levenberg-Marquardt algorithm.
--------------------------------------------------------------------------------

-- | The Levenberg-Marquardt algorithm specialised for curve-fitting.
levmar :: forall n k r a. (Nat n, Nat k, LMA_I.LevMarable r)
       => (Model n r a)                          -- ^ Model
       -> Maybe (Jacobian n r a)                 -- ^ Optional jacobian
       -> SizedList n r                          -- ^ Initial parameters
       -> [(a, r)]                               -- ^ Samples
       -> Integer                                -- ^ Maximum number of iterations
       -> LMA_I.Options r                          -- ^ Minimization options
       -> Maybe (SizedList n r)                  -- ^ Optional lower bounds
       -> Maybe (SizedList n r)                  -- ^ Optional upper bounds
       -> Maybe (LinearConstraints k n r)    -- ^ Optional linear constraints
       -> Maybe (SizedList n r)                  -- ^ Optional weights
       -> Either LMA_I.LevMarError (SizedList n r, LMA_I.Info r, CovarMatrix n r)
levmar model mJac params ys itMax opts mLowBs mUpBs mLinC mWghts =
    fmap convertResult $ LMA_I.levmar (convertModel model)
                                      (fmap convertJacob mJac)
                                      (toList params)
                                      ys
                                      itMax
                                      opts
                                      (fmap toList mLowBs)
                                      (fmap toList mUpBs)
                                      (fmap convertLinearConstraints mLinC)
                                      (fmap toList mWghts)
    where
      convertModel f = \ps   ->           f $* (unsafeFromList ps :: SizedList n r)
      convertJacob f = \ps x -> toList (((f $* (unsafeFromList ps :: SizedList n r)) x) :: SizedList n r)


-- The End ---------------------------------------------------------------------
