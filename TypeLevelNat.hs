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


-- | Type-level natural denoting zero
data Z = Z deriving Show

-- | Type-level natural denoting the /S/uccessor of another type-level natural.
newtype S n = S n deriving Show

-- | Class of all type-level naturals.
class Nat n where
   -- | Case analysis on natural numbers.
   caseNat :: forall r.
              n                                      -- ^ The natural number to case analyse.
           -> (n ~ Z => r)                           -- ^ The result @r@ when @n@ equals zero.
           -> (forall p. (n ~ S p, Nat p) => p -> r) -- ^ Function to apply to the predecessor
                                                     --   of @n@ to yield the result @r@.
           -> r

instance Nat Z where
   caseNat _ z _ = z

instance Nat p => Nat (S p) where
   caseNat (S p) _ s = s p

-- | The axiom of induction on natural numbers.
-- See: <http://en.wikipedia.org/wiki/Mathematical_induction#Axiom_of_induction>
induction :: forall p n. Nat n
          => n
          -> p Z
          -> (forall m. Nat m => p m -> p (S m))
          -> p n
induction n z s = caseNat n isZ isS
    where
      isZ :: n ~ Z => p n
      isZ = z

      isS :: forall m. (n ~ S m, Nat m) => m -> p n
      isS m = s (induction m z s)

newtype Witness x = Witness { unWitness :: x }

-- | The value of @witnessNat :: n@ is the natural number of type @n@.
-- For example:
--
-- @
-- *TypeLevelNat> witnessNat :: S (S (S Z))
-- S (S (S Z))
-- @
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

