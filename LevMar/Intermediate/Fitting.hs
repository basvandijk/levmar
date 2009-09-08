module LevMar.Intermediate.Fitting
    ( Model
    , Jacobian
    , levmar

    , LMA_I.Options(..)
    , LMA_I.defaultOpts

    , LMA_I.LinearConstraints

    , LMA_I.Info(..)
    , LMA_I.StopReason(..)
    , LMA_I.CovarMatrix
    , LMA_I.LevMarError(..)

    , LMA_I.LevMarable
    ) where

import qualified LevMar.Intermediate as LMA_I

type Model r a = [r] -> a -> r

-- | See: <http://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant>
type Jacobian r a = [r] -> a -> [r]

levmar :: LMA_I.LevMarable r
       => Model r a                         -- ^ Model
       -> Maybe (Jacobian r a)              -- ^ Optional jacobian
       -> [r]                               -- ^ Initial parameters
       -> [(a, r)]                          -- ^ Samples
       -> Integer                           -- ^ Maximum iterations
       -> LMA_I.Options r                   -- ^ Options
       -> Maybe [r]                         -- ^ Optional lower bounds
       -> Maybe [r]                         -- ^ Optional upper bounds
       -> Maybe (LMA_I.LinearConstraints r) -- ^ Optional linear constraints
       -> Maybe [r]                         -- ^ Optional weights
       -> Either LMA_I.LevMarError ([r], LMA_I.Info r, LMA_I.CovarMatrix r)
levmar model mJac ps samples =
    LMA_I.levmar (\ps' -> map (model ps') xs)
                 (fmap (\jac -> \ps' -> map (jac ps') xs) mJac)
                 ps
                 ys
        where
          (xs, ys) = unzip samples
