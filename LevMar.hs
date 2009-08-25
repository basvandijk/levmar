{-# LANGUAGE GADTs #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE EmptyDataDecls  #-}
{-# LANGUAGE RankNTypes  #-}
{-# LANGUAGE ScopedTypeVariables #-}

module LevMar
    ( dlevmar_dif
    , slevmar_dif
    , LevMarDif
    , ModelFunc
    , Vector(..)
    , Z, S
    , ($*)
    , CovarMatrix

    , LI.Options(..)
    , LI.defaultOpts
    , LI.StopReason(..)
    , LI.Info(..)
    ) where

import qualified LevMar.Intermediate as LI

import Unsafe.Coerce (unsafeCoerce)

data Z
data S n

data Vector n a where
   Nil   :: Vector Z a
   (:*:) :: a -> Vector n a -> Vector (S n) a

infixr 5 :*:

consPrecedence :: Int
consPrecedence = 5

instance Show a => Show (Vector n a) where
   showsPrec _ Nil        = showString "Nil"
   showsPrec p (x :*: xs) = showParen (p > consPrecedence) $
                            showsPrec (consPrecedence + 1) x .
                            showString " :*: "               .
                            showsPrec consPrecedence xs

toList :: Vector n a -> [a]
toList Nil        = []
toList (x :*: xs) = x : toList xs

vectorCPS :: [a] -> (forall n. Vector n a -> t) -> t
vectorCPS []       f = f Nil
vectorCPS (x : xs) f = vectorCPS xs (\ys -> f (x :*: ys))

-- | You have to make sure that the length of the input list equals
-- the length in the type of the output vector.
reallyUnsafeFromList :: [a] -> Vector n a
reallyUnsafeFromList xs = vectorCPS xs unsafeCoerce

type family ModelFunc n a :: *

type instance ModelFunc Z     a = a
type instance ModelFunc (S n) a = a -> ModelFunc n a

($*) :: ModelFunc n a -> Vector n a -> a
f $* Nil        = f
f $* (x :*: xs) = f x $* xs

type CovarMatrix n r = Vector n (Vector n r)

type LevMarDif n r a =  (a -> ModelFunc n r)
                     -> Vector n r
                     -> [(a, r)]
                     -> Integer
                     -> LI.Options r
                     -> (Vector n r, LI.Info r, CovarMatrix n r)

gen_levmar_dif :: forall n r a. LI.LevMarDif r a -> LevMarDif n r a
gen_levmar_dif levmar_dif model params samples itMax opts =
    let (psResult, info, covar) = levmar_dif (\x ps -> model x $* (reallyUnsafeFromList ps :: Vector n r))
                                             (toList params)
                                             samples
                                             itMax
                                             opts
    in ( reallyUnsafeFromList psResult
       , info
       , reallyUnsafeFromList $ map reallyUnsafeFromList covar
       )

dlevmar_dif :: LevMarDif n Double a
dlevmar_dif = gen_levmar_dif LI.dlevmar_dif

slevmar_dif :: LevMarDif n Float a
slevmar_dif = gen_levmar_dif LI.slevmar_dif
