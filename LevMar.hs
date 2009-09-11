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
    , NFunction
    )
    where


import qualified LevMar.Intermediate as LMA_I

import TypeLevelNat (Z, S, Nat)
import SizedList    (SizedList(..), toList, unsafeFromList)
import NFunction    (NFunction, ($*))

import Data.Either


--------------------------------------------------------------------------------
-- Model & Jacobian.
--------------------------------------------------------------------------------

{- | A functional relation describing measurements represented as a function
from @n@ parameters of type @r@ to a list of @r@.

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
type Model n r = NFunction n r [r]

{- | The jacobian of the 'Model' function. Expressed as a function from
@n@ parameters of type @r@ to a list of @n@-sized lists of @r@.

See: <http://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant>

 * Ensure that the length of the ouput list equals the length of the sample list
   in 'levmar'.

For example the jacobian of the above @hatfldc@ model is:

@
type N4 = 'S' ('S' ('S' ('S' 'Z')))

hatfldc_jac :: Jacobian N4 Double
hatfldc_jac _ p1 p2 _ = [ 1.0 ::: 0.0            ::: 0.0            ::: 0.0 ::: Nil
                        , 1.0 ::: -0.5 / sqrt p1 ::: 0.0            ::: 0.0 ::: Nil
                        , 0.0 ::: 1.0            ::: -0.5 / sqrt p2 ::: 0.0 ::: Nil
                        , 0.0 ::: 0.0            ::: 0.0            ::: 1.0 ::: Nil
                        ]
@
-}

type Jacobian n r = NFunction n r [SizedList n r]


--------------------------------------------------------------------------------
-- Levenberg-Marquardt algorithm.
--------------------------------------------------------------------------------

-- | The Levenberg-Marquardt algorithm.
levmar :: forall n k r. (Nat n, Nat k, LMA_I.LevMarable r)
       => (Model n r)                     -- ^ Model
       -> Maybe (Jacobian n r)            -- ^ Optional jacobian
       -> SizedList n r                   -- ^ Initial parameters
       -> [r]                             -- ^ Samples
       -> Integer                         -- ^ Maximum number of iterations
       -> LMA_I.Options r                 -- ^ Minimization options
       -> Maybe (SizedList n r)           -- ^ Optional lower bounds
       -> Maybe (SizedList n r)           -- ^ Optional upper bounds
       -> Maybe (LinearConstraints k n r) -- ^ Optional linear constraints
       -> Maybe (SizedList n r)           -- ^ Optional weights
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
                                      (fmap convertLinC mLinC)
                                      (fmap toList mWghts)
    where
      convertModel f = \ps ->              f $* (unsafeFromList ps :: SizedList n r)
      convertJacob f = \ps -> map toList ((f $* (unsafeFromList ps :: SizedList n r)) :: [SizedList n r])
      convertLinC (cMat, rhcVec) = ( map toList $ toList cMat
                                   , toList rhcVec
                                   )
      convertResult (psResult, info, covar) = ( unsafeFromList psResult
                                              , info
                                              , unsafeFromList $ map unsafeFromList covar
                                              )

-- | Linear constraints consisting of a constraints matrix, /kxn/ and
--   a right hand constraints vector, /kx1/ where /n/ is the number of
--   parameters and /k/ is the number of constraints.
type LinearConstraints k n r = (Matrix k n r, SizedList k r)

-- |Value to denote the absense of any linear constraints over the
-- parameters of the model function. Use this instead of 'Nothing'
-- because the type parameter which contains the number of constraints
-- can't be inferred.
noLinearConstraints :: Nat n => Maybe (LinearConstraints Z n r)
noLinearConstraints = Nothing

-- | A /nxm/ matrix is a sized list of /n/ sized lists of length /m/.
type Matrix n m r = SizedList n (SizedList m r)

-- | Covariance matrix corresponding to LS solution.
type CovarMatrix n r = Matrix n n r


-- The End ---------------------------------------------------------------------
