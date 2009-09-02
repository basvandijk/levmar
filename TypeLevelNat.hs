-- Thanks to Ryan Ingram who wrote most of this module.
-- See: http://www.haskell.org/pipermail/haskell-cafe/2009-August/065674.html

{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies #-}

module TypeLevelNat
    ( Z(..)
    , S(..)
    , Nat
    , caseNat
    , induction
    , witnessNat

    , N(..)
    ) where


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

-- | A value-level natural indexed with an equivalent type-level natural.
data N n where
    Zero :: N Z
    Succ :: N n -> N (S n)

{-
Template Haskell code to construct a type synonym for an arbitrary
type level natural number.

Instead of

> type N6 = S (S (S (S (S (S Z)))))

you can write

> $(mkNat "N6" 6)
-}

-- import Language.Haskell.TH.Syntax

-- mkNat :: String -> Int -> Q [Dec]
-- mkNat syn = runQ . return . (:[]) . TySynD (mkName syn) [] . go
--     where go 0 = ConT $ mkName "Z"
--           go n = AppT (ConT $ mkName "S") $ go (n - 1)

