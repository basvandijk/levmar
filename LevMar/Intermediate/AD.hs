{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleContexts #-}

--------------------------------------------------------------------------------
-- |
-- Module      :  LevMar.Intermediate.AD
-- Copyright   :  (c) 2009 Roel van Dijk & Bas van Dijk
-- License     :  BSD-style (see the file LICENSE)
--
-- Maintainer  :  vandijk.roel@gmail.com, v.dijk.bas@gmail.com
-- Stability   :  Experimental
--
-- A levmar variant that uses Automatic Differentiation to
-- automatically compute the Jacobian.
--
-- For additional documentation see the documentation of the levmar C
-- library which this library is based on:
-- <http://www.ics.forth.gr/~lourakis/levmar/>
--
--------------------------------------------------------------------------------

module LevMar.Intermediate.AD
    ( -- * Model.
      LMA_I.Model
    , LMA_I.Jacobian
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


import qualified LevMar.Intermediate as LMA_I

import LevMar.Utils.AD  ( value, firstDeriv, constant, idDAt )

-- From vector-space:
import Data.Derivative  ( (:~>) )
import Data.VectorSpace ( VectorSpace, Scalar )
import Data.Basis       ( HasBasis, Basis )

import Data.List        ( transpose )


--------------------------------------------------------------------------------
-- Levenberg-Marquardt algorithm.
--------------------------------------------------------------------------------

-- | The Levenberg-Marquardt algorithm that automatically computes the
-- 'Jacobian' using automatic differentiation of the model function.
--
-- /Warning/: Don't apply 'levmar' to 'LMA_I.Model's that apply methods of
-- the 'Eq' and 'Ord' classes to the parameters. These methods are
-- undefined for ':~>'!!!
levmar :: forall r.
          ( HasBasis r
          , Basis r ~ ()
          , VectorSpace (Scalar r)
          , LMA_I.LevMarable r
          )
       => LMA_I.Model (r :~> r)             -- ^ Model. Note that
                                            --   ':~>' is overloaded
                                            --   for all the numeric
                                            --   classes.
       -> [r]                               -- ^ Initial parameters
       -> [r]                               -- ^ Samples
       -> Integer                           -- ^ Maximum iterations
       -> LMA_I.Options r                   -- ^ Minimization options
       -> Maybe [r]                         -- ^ Optional lower bounds
       -> Maybe [r]                         -- ^ Optional upper bounds
       -> Maybe (LMA_I.LinearConstraints r) -- ^ Optional linear constraints
       -> Maybe [r]                         -- ^ Optional weights
       -> Either LMA_I.LevMarError ([r], LMA_I.Info r, LMA_I.CovarMatrix r)

levmar model = LMA_I.levmar (convertModel model) . Just $ jacobianOf model
    where
      convertModel :: LMA_I.Model (r :~> r) -> LMA_I.Model r
      convertModel mdl = map value . mdl . map constant

-- | Compute the 'LMA_I.Jacobian' of the 'LMA_I.Model' using Automatic
-- Differentiation.
jacobianOf :: (HasBasis r, Basis r ~ (), VectorSpace (Scalar r))
           => LMA_I.Model (r :~> r) -> LMA_I.Jacobian r
(jacobianOf mdl) ps = map (\fs -> zipWith (firstDeriv .) fs ps)
                    . transpose $ map mdl pDs
    where
      pDs = [idDAt n ps | n <- [0 .. length ps - 1]]


-- The End ---------------------------------------------------------------------
