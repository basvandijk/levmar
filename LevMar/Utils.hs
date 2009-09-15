module LevMar.Utils
    ( LinearConstraints
    , noLinearConstraints
    , Matrix
    , CovarMatrix
    , convertLinearConstraints
    , convertResult
    ) where

import qualified LevMar.Intermediate as LMA_I

import TypeLevelNat ( Nat, Z )
import SizedList    ( SizedList, toList, unsafeFromList )

-- | Linear constraints consisting of a constraints matrix, /kxn/ and
--   a right hand constraints vector, /kx1/ where /n/ is the number of
--   parameters and /k/ is the number of constraints.
type LinearConstraints k n r = (Matrix k n r, SizedList k r)

-- |Value to denote the absense of any linear constraints over the
-- parameters of the model function. Use this instead of 'Nothing'
-- because the type parameter which contains the number of constraints
-- can't be inferred.
noLinearConstraints :: Nat n => Maybe (LinearConstraints Z n r)
noLinearConstraints = Nothing

-- | A /nxm/ matrix is a sized list of /n/ sized lists of length /m/.
type Matrix n m r = SizedList n (SizedList m r)

-- | Covariance matrix corresponding to LS solution.
type CovarMatrix n r = Matrix n n r

convertLinearConstraints :: (Nat k, Nat n) => LinearConstraints k n r -> LMA_I.LinearConstraints r
convertLinearConstraints (cMat, rhcVec) = ( map toList $ toList cMat
                                          , toList rhcVec
                                          )

convertResult :: (Nat n)
              => ([r],           LMA_I.Info r, LMA_I.CovarMatrix r)
              -> (SizedList n r, LMA_I.Info r, CovarMatrix n r)
convertResult (psResult, info, covar) = ( unsafeFromList psResult
                                        , info
                                        , unsafeFromList $ map unsafeFromList covar
                                        )
