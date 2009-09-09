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

type Model    n r = NFunction n r [r]

-- | See: <http://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant>
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
