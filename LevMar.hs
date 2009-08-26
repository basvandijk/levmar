-- Thanks to Ryan Ingram who wrote most of this module.
-- See: http://www.haskell.org/pipermail/haskell-cafe/2009-August/065674.html

{-# LANGUAGE GADTs #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE FlexibleContexts #-}

module LevMar
    ( levmar

    , LevMar
    , ParamFunc
    , Model
    , Jacobian
    , Z(..), S(..)
    , SizedList(..)
    , Nat(..)
    , CovarMatrix

    , LMA_I.Options(..)
    , LMA_I.defaultOpts
    , LMA_I.StopReason(..)
    , LMA_I.Info(..)

    -- TODO: The following is only needed for the test:
    , ($*)
    )
    where

import qualified LevMar.Intermediate as LMA_I

import Data.Maybe (fromJust)

data Z = Z
newtype S n = S n

class Nat n where
   caseNat :: forall r.
              n
           -> (n ~ Z => r)
           -> (forall p. (n ~ S p, Nat p) => p -> r)
           -> r

instance Nat Z where
   caseNat _ z _ = z

instance Nat n => Nat (S n) where
   caseNat (S n) _ s = s n

induction :: forall p n. Nat n
          => n
          -> p Z
          -> (forall x. Nat x => p x -> p (S x))
          -> p n
induction n z s = caseNat n isZ isS
    where
      isZ :: n ~ Z => p n
      isZ = z

      isS :: forall x. (n ~ S x, Nat x) => x -> p n
      isS x = s (induction x z s)

newtype Witness x = Witness { unWitness :: x }

witnessNat :: forall n. Nat n => n
witnessNat = theWitness
    where
      theWitness = unWitness $ induction (undefined `asTypeOf` theWitness)
                                         (Witness Z)
                                         (Witness . S . unWitness)

data SizedList a n where
   Nil   :: SizedList a Z
   (:::) :: a -> SizedList a n -> SizedList a (S n)

infixr 5 ::: -- Same precedence and associativity as (:)

consPrecedence :: Int
consPrecedence = 5

instance Show a => Show (SizedList a n) where
    showsPrec _ Nil        = showString "Nil"
    showsPrec p (x ::: xs) = showParen (p > consPrecedence)
                           $ showsPrec (consPrecedence + 1) x
                           . showString " ::: "
                           . showsPrec consPrecedence xs

newtype ToList a n = ToList { unToList :: SizedList a n -> [a] }

toList :: forall a n. Nat n => SizedList a n -> [a]
toList = unToList $ induction (witnessNat :: n)
                              (ToList tl0)
                              (ToList . tlS . unToList)
    where
      tl0 :: SizedList a Z -> [a]
      tl0 Nil = []

      tlS :: forall x. Nat x => (SizedList a x -> [a]) -> SizedList a (S x) -> [a]
      tlS f (x ::: xs) = x : f xs

newtype FromList a n = FromList { unFromList :: [a] -> Maybe (SizedList a n) }

fromList :: forall a n. Nat n => [a] -> Maybe (SizedList a n)
fromList = unFromList $ induction (witnessNat :: n)
                                  (FromList fl0)
                                  (FromList . flS . unFromList)
    where
      fl0 [] = Just Nil
      fl0 _  = Nothing

      flS _ []     = Nothing
      flS k (x:xs) = fmap (x :::) $ k xs

unsafeFromList :: forall a n. Nat n => [a] -> SizedList a n
unsafeFromList = fromJust . fromList

-- | @ParamFunc r n@ represents a function from @n@ @r@'s to a @r@.
-- For example: @ParamFunc Double (S (S (S Z))) ~ Double -> Double -> Double -> Double@
type family ParamFunc r n :: *

type instance ParamFunc r Z     = r
type instance ParamFunc r (S n) = r -> ParamFunc r n

($*) :: ParamFunc r n -> SizedList r n -> r
f $* Nil        = f
f $* (x ::: xs) = f x $* xs

type Model a r n = a -> ParamFunc r n

type Jacobian a r n = Model a r n

type CovarMatrix r n = SizedList (SizedList r n) n

type LevMar a r n =  (Model a r n)          -- model
                  -> Maybe (Jacobian a r n) -- jacobian
                  -> SizedList r n          -- init params
                  -> [(a, r)]               -- samples
                  -> Integer                -- max iterations
                  -> LMA_I.Options r
                  -> Maybe (SizedList r n)  -- lower bounds
                  -> Maybe (SizedList r n)  -- upper bounds
                  -> (SizedList r n, LMA_I.Info r, CovarMatrix r n)

levmar :: forall n r a. (Nat n, LMA_I.LevMarable r) => LevMar a r n
levmar model mJac params samples itMax opts mLowBs mUpBs =
    let mkModel f = \x ps -> f x $* (unsafeFromList ps :: SizedList r n)
        (psResult, info, covar) = LMA_I.levmar (mkModel model)
                                               (fmap mkModel mJac)
                                               (toList params)
                                               samples
                                               itMax
                                               opts
                                               (fmap toList mLowBs)
                                               (fmap toList mUpBs)
    in ( unsafeFromList psResult
       , info
       , unsafeFromList $ map unsafeFromList covar
       )
