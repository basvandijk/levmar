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
      LMA.Model

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


import qualified LevMar              as LMA
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

import LevMar.Utils.AD  ( value, firstDeriv, constant, idDAt )

-- From vector-space:
import Data.Derivative  ( (:~>), (:>), powVal )
import Data.VectorSpace ( VectorSpace, Scalar )
import Data.Basis       ( HasBasis, Basis )

import Data.List        ( transpose )


--------------------------------------------------------------------------------
-- Levenberg-Marquardt algorithm.
--------------------------------------------------------------------------------

-- | The Levenberg-Marquardt algorithm that automatically computes the
-- 'Jacobian' using automatic differentiation of the model function.
--
-- /Warning/: Don't apply 'levmar' to 'LMA.Model's that apply methods of
-- the 'Eq' and 'Ord' classes to the parameters. These methods are
-- undefined for ':~>'!!!
levmar :: forall m n k r.
          ( Nat m
          , Nat n
          , Nat k
          , HasBasis r
          , Basis r ~ ()
          , VectorSpace (Scalar r)
          , LMA_I.LevMarable r
          )
       => (LMA.Model m n (r :~> r))       -- ^ Model. Note that ':~>'
                                          --   is overloaded for all the
                                          --   numeric classes.
       -> SizedList m r                   -- ^ Initial parameters
       -> SizedList n r                   -- ^ Samples
       -> Integer                         -- ^ Maximum number of iterations
       -> LMA_I.Options r                 -- ^ Minimization options
       -> Maybe (SizedList m r)           -- ^ Optional lower bounds
       -> Maybe (SizedList m r)           -- ^ Optional upper bounds
       -> Maybe (LinearConstraints k m r) -- ^ Optional linear constraints
       -> Maybe (SizedList m r)           -- ^ Optional weights
       -> Either LMA_I.LevMarError (SizedList m r, LMA_I.Info r, CovarMatrix m r)

levmar model params ys itMax opts mLowBs mUpBs mLinC mWghts =
    fmap convertResult $ LMA_I.levmar (convertModel model)
                                      (Just $ jacobianOf model)
                                      (toList params)
                                      (toList ys)
                                      itMax
                                      opts
                                      (fmap toList mLowBs)
                                      (fmap toList mUpBs)
                                      (fmap convertLinearConstraints mLinC)
                                      (fmap toList mWghts)
    where
      convertModel :: LMA.Model m n (r :~> r) -> LMA_I.Model r
      (convertModel mdl) ps = fmap value $ toList
                              (mdl $* pDs :: SizedList n (r :~> r))
          where
            pDs :: SizedList m (r :~> r)
            pDs = unsafeFromList $ fmap constant ps

      jacobianOf :: LMA.Model m n (r :~> r) -> LMA_I.Jacobian r
      (jacobianOf mdl) ps = fmap (\fs -> zipWith (firstDeriv .) fs ps)
                          . transpose
                          . fmap (\pD -> toList (mdl $* (pD :: SizedList m (r :~> r)) :: SizedList n (r :~> r)))
                          $ pDs
          where
            pDs :: [SizedList m (r :~> r)]
            pDs = [unsafeFromList $ idDAt n ps | n <- [0 .. length ps - 1]]


-- The End ---------------------------------------------------------------------
