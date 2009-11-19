{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ScopedTypeVariables #-}

module NFunction
    ( NFunction
    , uncurry, ($*)
    , Curry, curry
    , ComposeN, compose
    ) where


import Prelude hiding ( curry, uncurry )
import TypeLevelNat   ( Z(..), S(..), Nat )
import SizedList      ( SizedList(..) )


--------------------------------------------------------------------------------

-- | A @NFunction n a b@ is a function which takes @n@ arguments of
-- type @a@ and returns a @b@.
-- For example: @NFunction (S (S (S Z))) a b ~ (a -> a -> a -> b)@
type family NFunction n a b :: *

type instance NFunction Z     a b = b
type instance NFunction (S n) a b = a -> NFunction n a b


--------------------------------------------------------------------------------

-- | @uncurry f xs@ applies the /n/-arity function @f@ to each of the arguments in
-- the /n/-sized list @xs@.
uncurry :: NFunction n a b -> SizedList n a -> b
uncurry f Nil        = f
uncurry f (x ::: xs) = uncurry (f x) xs

-- | Infix version of 'uncurry'.
($*) :: NFunction n a b -> SizedList n a -> b
($*) = uncurry

infixr 0 $* -- same as $


--------------------------------------------------------------------------------

class Nat n => Curry n where
    curry :: (SizedList n a -> b) -> NFunction n a b

instance Curry Z where
    curry = ($ Nil)

instance Curry n => Curry (S n) where
    curry f = \x -> curry $ f . (x :::)


--------------------------------------------------------------------------------

class Nat n => ComposeN n where
    -- | Composition of @NFunction@s.
    --
    -- Note that the @n@ and @a@ arguments are used by the type
    -- checker to select the right @ComposeN@ instance. They are
    -- usally given as @(witnessNat :: n)@ and @(undefined :: a)@.
    compose :: forall a b c. n -> a
            -> (b -> c) -> NFunction n a b -> NFunction n a c

instance ComposeN Z where
    compose Z _ = ($)

instance ComposeN n => ComposeN (S n) where
    compose (S n) (_ :: a) f g = compose n (undefined :: a) f . g


--------------------------------------------------------------------------------

{-
TODO: The following does not work as expected.
See: http://www.haskell.org/pipermail/haskell-cafe/2009-August/065850.html

-- | @f .* g@ composes @f@ with the /n/-arity function @g@.
(.*) :: forall n a b c. (ComposeN n) => (b -> c) -> NFunction n a b -> NFunction n a c
(.*) = compose (witnessNat :: n) (undefined :: a)

infixr 9 .* -- same as .
-}


-- The End ---------------------------------------------------------------------
