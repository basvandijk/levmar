{-# LANGUAGE GADTs #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE Rank2Types #-}

module SizedList
    ( SizedList(..)
    , foldr
    , foldrN
    , toList
    , length
    , fromList
    , unsafeFromList
    , replicate
    ) where


import Prelude hiding ( foldr, replicate, length )
import Data.Maybe     ( fromMaybe )
import TypeLevelNat   ( Z(..), S(..), Nat, induction, witnessNat, N(..) )


--------------------------------------------------------------------------------

-- | A list which is indexed with a type-level natural that denotes the size of
-- the list.
data SizedList n a where
   Nil   :: SizedList Z a
   (:::) :: a -> SizedList n a -> SizedList (S n) a

instance Functor (SizedList n) where
    fmap _ Nil        = Nil
    fmap f (x ::: xs) = f x ::: fmap f xs

infixr 5 ::: -- Same precedence and associativity as (:)


--------------------------------------------------------------------------------

consPrecedence :: Int
consPrecedence = 5

instance Show a => Show (SizedList n a) where
    showsPrec _ Nil        = showString "Nil"
    showsPrec p (x ::: xs) = showParen (p > consPrecedence)
                           $ showsPrec (consPrecedence + 1) x
                           . showString " ::: "
                           . showsPrec consPrecedence xs


--------------------------------------------------------------------------------

-- | Fold a binary operator over a @SizedList@.
foldr :: forall a b n. (a -> b -> b) -> b -> SizedList n a -> b
foldr f z = foldr_f_z
    where
      foldr_f_z :: forall k. SizedList k a -> b
      foldr_f_z Nil        = z
      foldr_f_z (x ::: xs) = f x $ foldr_f_z xs

-- | Fold a binary operator yielding a value with a natural number
-- indexed type over a @SizedList@.
foldrN :: forall a b n. (forall m. a -> b m -> b (S m)) -> b Z -> SizedList n a -> b n
foldrN f z = foldrN_f_z
    where
      foldrN_f_z :: forall k. SizedList k a -> b k
      foldrN_f_z Nil        = z
      foldrN_f_z (x ::: xs) = f x $ foldrN_f_z xs

-- | Convert a @SizedList@ to a normal list.
toList :: SizedList n a -> [a]
toList = foldr (:) []

-- | Returns the length of the @SizedList@.
length :: SizedList n a -> N n
length = foldrN (const Succ) Zero


--------------------------------------------------------------------------------

newtype FromList a n = FL { unFL :: [a] -> Maybe (SizedList n a) }

-- | Convert a normal list to a @SizedList@. If the length of the given
-- list does not equal @n@, @Nothing@ is returned.
fromList :: forall a n. Nat n => [a] -> Maybe (SizedList n a)
fromList = unFL $ induction (witnessNat :: n) (FL flZ) (FL . flS . unFL)
    where
      flZ [] = Just Nil
      flZ _  = Nothing

      flS _ []     = Nothing
      flS k (x:xs) = fmap (x :::) $ k xs

-- | Convert a normal list to a @SizeList@. If the length of the given
-- list does not equal @n@, an error is thrown.
unsafeFromList :: forall a n. Nat n => [a] -> SizedList n a
unsafeFromList = fromMaybe (error "unsafeFromList xs: xs does not have the right length ") .
                 fromList


--------------------------------------------------------------------------------

newtype Replicate a n = R { unR :: SizedList n a}

-- | @replicate x :: SizedList n a@ returns a @SizedList@ of @n@ @x@s.
replicate :: forall a n. Nat n => a -> SizedList n a
replicate x = unR $ induction (witnessNat :: n) (R Nil) (R . (x :::) . unR)


-- The End ---------------------------------------------------------------------
