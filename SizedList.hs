{-# LANGUAGE GADTs #-}
{-# LANGUAGE ScopedTypeVariables #-}

module SizedList
    ( SizedList(..)
    , toList
    , fromList
    , unsafeFromList
    , lengthSL
    , replicateSL
    ) where

import Data.Maybe (fromJust)

import TypeLevelNat (Z(..), S(..), Nat, induction, witnessNat, N(..))

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

replicateSL :: N n -> a -> SizedList a n
replicateSL Zero     _ = Nil
replicateSL (Succ n) x = x ::: replicateSL n x

lengthSL :: SizedList a n -> N n
lengthSL Nil        = Zero
lengthSL (_ ::: xs) = Succ (lengthSL xs)
