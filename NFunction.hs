{-# LANGUAGE TypeFamilies #-}

module NFunction
    ( NFunction
    , ($*)
    ) where

import TypeLevelNat (Z, S)
import SizedList    (SizedList(..))

-- | A @NFunction n a b@ is a function which takes @n@ arguments of
-- type @a@ and returns a @b@.
-- For example: NFunction (S (S (S Z))) a b ~ (a -> a -> a -> b)
type family NFunction n a b :: *

type instance NFunction Z     a b = b
type instance NFunction (S n) a b = a -> NFunction n a b

-- | @f $* ps@ applies the /n/-arity function @f@ to each of the arguments in
-- the /n/-sized list @xs@.
($*) :: NFunction n a b -> SizedList n a -> b
f $* Nil        = f
f $* (x ::: xs) = f x $* xs

infixr 0 $* -- same as $
