{-# LANGUAGE GADTs #-}
{-# LANGUAGE ScopedTypeVariables #-}

module SizedList
    ( SizedList(..)
    , toList
    , fromList
    , unsafeFromList
    , length
    , replicate
    ) where

import Prelude hiding (replicate, length)
import Data.Maybe     (fromJust)
import TypeLevelNat   (Z(..), S(..), Nat, induction, witnessNat, N(..))

-- | A list which is indexed with a type-level natural that denotes the size of
-- the list.
data SizedList n a where
   Nil   :: SizedList Z a
   (:::) :: a -> SizedList n a -> SizedList (S n) a

infixr 5 ::: -- Same precedence and associativity as (:)

consPrecedence :: Int
consPrecedence = 5

instance Show a => Show (SizedList n a) where
    showsPrec _ Nil        = showString "Nil"
    showsPrec p (x ::: xs) = showParen (p > consPrecedence)
                           $ showsPrec (consPrecedence + 1) x
                           . showString " ::: "
                           . showsPrec consPrecedence xs

newtype ToList a n = ToList { unToList :: SizedList n a -> [a] }

toList :: forall a n. Nat n => SizedList n a -> [a]
toList = unToList $ induction (witnessNat :: n)
                              (ToList tl0)
                              (ToList . tlS . unToList)
    where
      tl0 :: SizedList Z a -> [a]
      tl0 Nil = []

      tlS :: forall x. Nat x => (SizedList x a -> [a]) -> SizedList (S x) a -> [a]
      tlS f (x ::: xs) = x : f xs

newtype FromList a n = FromList { unFromList :: [a] -> Maybe (SizedList n a) }

fromList :: forall a n. Nat n => [a] -> Maybe (SizedList n a)
fromList = unFromList $ induction (witnessNat :: n)
                                  (FromList fl0)
                                  (FromList . flS . unFromList)
    where
      fl0 [] = Just Nil
      fl0 _  = Nothing

      flS _ []     = Nothing
      flS k (x:xs) = fmap (x :::) $ k xs

unsafeFromList :: forall a n. Nat n => [a] -> SizedList n a
unsafeFromList = fromJust . fromList

replicate :: N n -> a -> SizedList n a
replicate Zero     _ = Nil
replicate (Succ n) x = x ::: replicate n x

length :: SizedList n a -> N n
length Nil        = Zero
length (_ ::: xs) = Succ $ length xs
