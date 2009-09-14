{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleContexts #-}

module LevMar.Utils.AD where

import Data.Derivative  ((:~>), (:>), powVal, idD, pureD, derivAtBasis)
import Data.VectorSpace (VectorSpace, Scalar, AdditiveGroup)
import Data.Basis       (HasBasis, Basis)
import Data.MemoTrie    (HasTrie)

-- | @firstDeriv f@ returns the first derivative of @f@.
firstDeriv :: (HasBasis a, Basis a ~ (), AdditiveGroup b)
           => (a :> b) -> b
firstDeriv f = powVal $ derivAtBasis f ()

-- | A constant infinitely differentiable function.
constant :: (AdditiveGroup b, HasBasis a, HasTrie (Basis a))
         => b -> a:~>b
constant = const . pureD

-- | @idDAt n ps@ maps each parameter in @ps@ to a /constant/
-- infinitely differentiable function (@const . pureD@), except the @n@th
-- parameter is replaced with the differentiable /identity/ function
-- (@idD@).
idDAt :: (HasBasis r, HasTrie (Basis r), VectorSpace (Scalar r))
      => Int -> [r] -> [r :~> r]
idDAt n = replace n idD . map constant

-- | @replace i r xs@ replaces the @i@th element in @xs@ with @r@.
replace :: Int -> a -> [a] -> [a]
replace i r xs
    | i < 0     = xs
    | otherwise = rep i xs
  where rep _ [] = []
        rep i (x:xs)
          | i > 0     = x : rep (i - 1) xs
          | otherwise = r : xs
