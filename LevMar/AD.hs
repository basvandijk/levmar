{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleContexts #-}

--------------------------------------------------------------------------------
-- |
-- Module      :  LevMar.AD
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


module LevMar.AD
    ( -- * Model
      Model

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
    )
    where


import qualified LevMar.Intermediate as LMA_I

import LevMar.Utils ( LinearConstraints
                    , noLinearConstraints
                    , Matrix
                    , CovarMatrix
                    , convertLinearConstraints
                    , convertResult
                    )

import TypeLevelNat ( Z, S, Nat )
import SizedList    ( SizedList(..), toList, unsafeFromList )
import NFunction    ( NFunction, ($*) )

import LevMar.Utils.AD  ( firstDeriv, constant, idDAt )

-- From vector-space:
import Data.Derivative  ( (:~>), (:>), powVal )
import Data.VectorSpace ( VectorSpace, Scalar )
import Data.Basis       ( HasBasis, Basis )

import Data.List        ( transpose )


--------------------------------------------------------------------------------
-- Model & Jacobian.
--------------------------------------------------------------------------------

{- | A functional relation describing measurements represented as a
function from @n@ parameters of type @r :~> r@ to a list of @r :~> r@.

 * Ensure that the length of the ouput list equals the length of the sample list
   in 'levmar'.

An example from /Demo.hs/:

@
type N4 = 'S' ('S' ('S' ('S' 'Z')))

hatfldc :: Model N4 Double
hatfldc p0 p1 p2 p3 = [ p0 - 1.0
                      , p0 - sqrt p1
                      , p1 - sqrt p2
                      , p3 - 1.0
                      ]
@
-}
type Model n r = NFunction n (r :~> r) [r :~> r]


--------------------------------------------------------------------------------
-- Levenberg-Marquardt algorithm.
--------------------------------------------------------------------------------

-- | The Levenberg-Marquardt algorithm that automatically computes the
-- 'Jacobian' using automatic differentiation of the model function.
levmar :: forall n k r.
          ( Nat n
          , Nat k
          , HasBasis r
          , Basis r ~ ()
          , VectorSpace (Scalar r)
          , LMA_I.LevMarable r
          )
       => (Model n r)                     -- ^ Model
       -> SizedList n r                   -- ^ Initial parameters
       -> [r]                             -- ^ Samples
       -> Integer                         -- ^ Maximum number of iterations
       -> LMA_I.Options r                 -- ^ Minimization options
       -> Maybe (SizedList n r)           -- ^ Optional lower bounds
       -> Maybe (SizedList n r)           -- ^ Optional upper bounds
       -> Maybe (LinearConstraints k n r) -- ^ Optional linear constraints
       -> Maybe (SizedList n r)           -- ^ Optional weights
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
      convertModel :: Model n r -> LMA_I.Model r
      (convertModel f) ps = fmap (\m -> powVal $ m undefined) (f $* pDs :: [r :~> r])
          where
            pDs :: SizedList n (r :~> r)
            pDs = unsafeFromList $ fmap constant ps

      jacobianOf :: Model n r -> LMA_I.Jacobian r
      (jacobianOf f) ps = fmap (\fs -> zipWith (firstDeriv .) fs ps)
                        . transpose
                        . fmap (\pD -> f $* (pD :: SizedList n (r :~> r))::[r :~> r])
                        $ pDs
          where
            pDs :: [SizedList n (r :~> r)]
            pDs = [unsafeFromList $ idDAt n ps | n <- [0 .. length ps - 1]]


-- The End ---------------------------------------------------------------------
