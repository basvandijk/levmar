{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ScopedTypeVariables #-}

--------------------------------------------------------------------------------
-- |
-- Module      :  LevMar
-- Copyright   :  (c) 2009 Roel van Dijk & Bas van Dijk
-- License     :  BSD-style (see the file LICENSE)
--
-- Maintainer  :  vandijk.roel@gmail.com, v.dijk.bas@gmail.com
-- Stability   :  Experimental
--
--
--
-- For additional documentation see the documentation of the levmar C
-- library which this library is based on:
-- <http://www.ics.forth.gr/~lourakis/levmar/>
--
--------------------------------------------------------------------------------

module LevMar
    ( -- * Model & Jacobian.
      Model
    , Jacobian

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

import Data.Either


--------------------------------------------------------------------------------
-- Model & Jacobian.
--------------------------------------------------------------------------------

{- | A functional relation describing measurements represented as a function
from @m@ parameters to @n@ expected measurements.

An example from /Demo.hs/:

@
type N4 = 'S' ('S' ('S' ('S' 'Z')))

hatfldc :: Model N4 N4 Double
hatfldc (p0 ::: p1 ::: p2 ::: p3 ::: Nil) =     p0 - 1.0
                                            ::: p0 - sqrt p1
                                            ::: p1 - sqrt p2
                                            ::: p3 - 1.0
                                            ::: Nil
@
-}
type Model m n r = SizedList m r -> SizedList n r

{- | The jacobian of the 'Model' function. Expressed as a function
from @m@ parameters to a @n@/x/@m@ matrix which for each of the @n@
expected measurement describes the @m@ partial derivatives of the
parameters.

See: <http://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant>

For example the jacobian of the above @hatfldc@ model is:

@
type N4 = 'S' ('S' ('S' ('S' 'Z')))

hatfldc_jac :: Jacobian N4 N4 Double
hatfldc_jac (_ ::: p1 ::: p2 ::: _ ::: Nil) =     (1.0 ::: 0.0            ::: 0.0            ::: 0.0 ::: Nil)
                                              ::: (1.0 ::: -0.5 / sqrt p1 ::: 0.0            ::: 0.0 ::: Nil)
                                              ::: (0.0 ::: 1.0            ::: -0.5 / sqrt p2 ::: 0.0 ::: Nil)
                                              ::: (0.0 ::: 0.0            ::: 0.0            ::: 1.0 ::: Nil)
                                              ::: Nil
@
-}

type Jacobian m n r = SizedList m r -> Matrix n m r


--------------------------------------------------------------------------------
-- Levenberg-Marquardt algorithm.
--------------------------------------------------------------------------------

-- | The Levenberg-Marquardt algorithm.
levmar :: forall m n k r. (Nat m, Nat n, Nat k, LMA_I.LevMarable r)
       => (Model m n r)                   -- ^ Model
       -> Maybe (Jacobian m n r)          -- ^ Optional jacobian
       -> SizedList m r                   -- ^ Initial parameters
       -> SizedList n r                   -- ^ Samples
       -> Integer                         -- ^ Maximum number of iterations
       -> LMA_I.Options r                 -- ^ Minimization options
       -> Maybe (SizedList m r)           -- ^ Optional lower bounds
       -> Maybe (SizedList m r)           -- ^ Optional upper bounds
       -> Maybe (LinearConstraints k m r) -- ^ Optional linear constraints
       -> Maybe (SizedList m r)           -- ^ Optional weights
       -> Either LMA_I.LevMarError (SizedList m r, LMA_I.Info r, CovarMatrix m r)

levmar model mJac params ys itMax opts mLowBs mUpBs mLinC mWghts =
    fmap convertResult $ LMA_I.levmar (convertModel model)
                                      (fmap convertJacob mJac)
                                      (toList params)
                                      (toList ys)
                                      itMax
                                      opts
                                      (fmap toList mLowBs)
                                      (fmap toList mUpBs)
                                      (fmap convertLinearConstraints mLinC)
                                      (fmap toList mWghts)
    where
      convertModel f = toList .               f . unsafeFromList
      convertJacob f = toList . fmap toList . f . unsafeFromList


-- The End ---------------------------------------------------------------------
