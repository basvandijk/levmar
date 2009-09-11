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
      Model

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

import Data.Derivative  ((:~>), (:>), powVal, idD, pureD, derivAtBasis)
import Data.VectorSpace (VectorSpace, Scalar, AdditiveGroup)
import Data.Basis       (HasBasis, Basis)
import Data.MemoTrie    (HasTrie)

import Data.List        (transpose)


--------------------------------------------------------------------------------
-- Model
--------------------------------------------------------------------------------

{- | A functional relation describing measurements.

 * Ensure that the length of the parameters list equals the length of the initial
   parameters list in 'levmar'.

 * Ensure that the length of the ouput list equals the length of the samples list
   in 'levmar'.

For example:

@
hatfldc :: Model Double
hatfldc [p0, p1, p2, p3] = [ p0 - 1.0
                           , p0 - sqrt p1
                           , p1 - sqrt p2
                           , p3 - 1.0
                           ]
@
-}
type Model r = [r :~> r] -> [r :~> r]


--------------------------------------------------------------------------------
-- Levenberg-Marquardt algorithm.
--------------------------------------------------------------------------------

-- | The Levenberg-Marquardt algorithm that automatically computes the
-- 'Jacobian' using automatic differentiation of the model function.
levmar :: ( HasBasis r
          , Basis r ~ ()
          , VectorSpace (Scalar r)
          , LMA_I.LevMarable r
          )
       => Model r                           -- ^ Model
       -> [r]                               -- ^ Initial parameters
       -> [r]                               -- ^ Samples
       -> Integer                           -- ^ Maximum iterations
       -> LMA_I.Options r                   -- ^ Minimization options
       -> Maybe [r]                         -- ^ Optional lower bounds
       -> Maybe [r]                         -- ^ Optional upper bounds
       -> Maybe (LMA_I.LinearConstraints r) -- ^ Optional linear constraints
       -> Maybe [r]                         -- ^ Optional weights
       -> Either LMA_I.LevMarError ([r], LMA_I.Info r, LMA_I.CovarMatrix r)

levmar model = LMA_I.levmar (convertModel model) $ Just $ jacobianOf model

convertModel :: (HasBasis r, HasTrie (Basis r)) => Model r -> LMA_I.Model r
convertModel model = \ps -> map (\m -> powVal $ m undefined) $ model $ map constant ps

jacobianOf :: ( HasBasis a
              , Basis a ~ ()
              , HasTrie (Basis a)
              , AdditiveGroup a
              , VectorSpace (Scalar a)
              ) => Model a -> LMA_I.Jacobian a
jacobianOf model =
    \ps -> let pDs = [idDAt n ps | n <- [0 .. length ps - 1]]
           in map (\fs -> zipWith (\f p -> firstDeriv $ f p) fs ps) $
                  transpose $ map model pDs

-- | @idDAt n ps@ maps each parameter in @ps@ to a /constant/
-- infinitely differentiable function (@const . pureD@), except the @n@th
-- parameter is replaced with the differentiable /identity/ function
-- (@idD@).
idDAt ::(HasBasis r, HasTrie (Basis r), VectorSpace (Scalar r)) => Int -> [r] -> [r :~> r]
idDAt n = replace n idD . map constant

-- | @firstDeriv f@ returns the first derivative of @f@.
firstDeriv :: (HasBasis a, Basis a ~ (), AdditiveGroup b) => (a :> b) -> b
firstDeriv f = powVal $ derivAtBasis f ()

-- | A constant infinitely differentiable function.
constant :: (AdditiveGroup b, HasBasis a, HasTrie (Basis a)) => b -> a:~>b
constant = const . pureD

-- | @replace i r xs@ replaces the @i@th element in @xs@ with @r@.
replace :: Int -> a -> [a] -> [a]
replace i r xs
    | i < 0     = xs
    | otherwise = rep i xs
  where rep _ [] = []
        rep i (x:xs)
          | i > 0     = x : rep (i - 1) xs
          | otherwise = r : xs


-- The End ---------------------------------------------------------------------
