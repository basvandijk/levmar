{-# LANGUAGE CPP
           , NoImplicitPrelude
           , UnicodeSyntax
           , ScopedTypeVariables
           , DeriveDataTypeable
  #-}

--------------------------------------------------------------------------------
-- |
-- Module:     Numeric.LevMar
-- Copyright:  (c) 2009 - 2011 Roel van Dijk & Bas van Dijk
-- License:    BSD-style (see the file LICENSE)
-- Maintainer: Roel van Dijk <vandijk.roel@gmail.com>
--             Bas van Dijk <v.dijk.bas@gmail.com>
-- Stability:  Experimental
--
-- For additional documentation see the documentation of the levmar C
-- library which this library is based on:
-- <http://www.ics.forth.gr/~lourakis/levmar/>
--
--------------------------------------------------------------------------------

module Numeric.LevMar
    ( -- * Model & Jacobian.
      Model
    , Jacobian

      -- * Levenberg-Marquardt algorithm.
    , LevMarable(levmar)

      -- * Minimization options.
    , Options(..)
    , defaultOpts

      -- * Constraints
    , Constraints(..)
    , LinearConstraints

      -- * Output
    , Info(..)
    , StopReason(..)
    , LevMarError(..)
    ) where


-------------------------------------------------------------------------------
-- Imports
-------------------------------------------------------------------------------

-- from base:
import Control.Monad         ( return, mplus )
import Control.Exception     ( Exception )
import Data.Typeable         ( Typeable )
import Data.Either           ( Either(Left, Right) )
import Data.Function         ( ($) )
import Data.Functor          ( (<$>) )
import Data.Int              ( Int )
import Data.List             ( lookup, (++) )
import Data.Maybe            ( Maybe(Nothing, Just), isJust, fromJust, fromMaybe )
import Data.Monoid           ( Monoid(mempty, mappend) )
import Data.Ord              ( (<) )
import Foreign.Marshal.Array ( allocaArray, withArray, peekArray, copyArray )
import Foreign.Ptr           ( Ptr, nullPtr )
import Foreign.ForeignPtr    ( ForeignPtr, newForeignPtr_, withForeignPtr )
import Foreign.Storable      ( Storable )
import Prelude               ( Enum, Fractional, RealFrac, Float, Double
                             , toEnum, (-), (*), error, floor
                             )
import System.IO             ( IO )
import System.IO.Unsafe      ( unsafePerformIO )
import Text.Read             ( Read )
import Text.Show             ( Show, show )

#if __GLASGOW_HASKELL__ >= 605
import GHC.ForeignPtr        ( mallocPlainForeignPtrBytes )
import Prelude               ( undefined )
import Foreign.Storable      ( sizeOf )
#else
import Foreign.ForeignPtr    ( mallocForeignPtrArray )
#endif

#if __GLASGOW_HASKELL__ < 700
import Prelude               ( fromInteger )
#endif

-- from base-unicode-symbols:
import Data.Bool.Unicode     ( (∧), (∨) )
import Data.Eq.Unicode       ( (≢) )
import Data.Function.Unicode ( (∘) )

-- from hmatrix:
import Data.Packed.Vector ( Vector )
import Data.Packed.Matrix ( Matrix, Element, flatten, rows, reshape )

-- from vector:
import qualified Data.Vector.Storable as VS ( unsafeWith, length
                                            , unsafeFromForeignPtr
                                            , length
                                            )

-- from bindings-levmar:
import Bindings.LevMar ( c'LM_INFO_SZ

                       , withModel
                       , withJacobian

                       , c'LM_ERROR
                       , c'LM_ERROR_LAPACK_ERROR
                       , c'LM_ERROR_FAILED_BOX_CHECK
                       , c'LM_ERROR_MEMORY_ALLOCATION_FAILURE
                       , c'LM_ERROR_CONSTRAINT_MATRIX_ROWS_GT_COLS
                       , c'LM_ERROR_CONSTRAINT_MATRIX_NOT_FULL_ROW_RANK
                       , c'LM_ERROR_TOO_FEW_MEASUREMENTS
                       , c'LM_ERROR_SINGULAR_MATRIX
                       , c'LM_ERROR_SUM_OF_SQUARES_NOT_FINITE

                       , c'LM_INIT_MU
                       , c'LM_STOP_THRESH
                       , c'LM_DIFF_DELTA
                       )
import qualified Bindings.LevMar ( Model, Jacobian )

-- from levmar (this package):
import Bindings.LevMar.CurryFriendly ( LevMarDer
                                     , LevMarDif
                                     , LevMarBCDer
                                     , LevMarBCDif
                                     , LevMarLecDer
                                     , LevMarLecDif
                                     , LevMarBLecDer
                                     , LevMarBLecDif
                                     , dlevmar_der,      slevmar_der
                                     , dlevmar_dif,      slevmar_dif
                                     , dlevmar_bc_der,   slevmar_bc_der
                                     , dlevmar_bc_dif,   slevmar_bc_dif
                                     , dlevmar_lec_der,  slevmar_lec_der
                                     , dlevmar_lec_dif,  slevmar_lec_dif
                                     , dlevmar_blec_der, slevmar_blec_der
                                     , dlevmar_blec_dif, slevmar_blec_dif
                                     )


--------------------------------------------------------------------------------
-- Model & Jacobian.
--------------------------------------------------------------------------------

{-| A functional relation describing measurements represented as a function
from a vector of parameters to a vector of expected measurements.

 * Ensure that the length of the parameters vector equals the length of the
   initial parameters vector in 'levmar'.

 * Ensure that the length of the ouput vector equals the length of the samples
   vector in 'levmar'.
-}
type Model r = Vector r → Vector r

{-| The jacobian of the 'Model' function. Expressed as a function from a vector
of parameters to a matrix which for each expected measurement describes
the partial derivatives of the parameters.

See: <http://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant>

 * Ensure that the length of the parameter vector equals the length of the initial
   parameter vector in 'levmar'.

 * Ensure that the output matrix has the dimension @n><m@ where @n@ is the
   number of samples and @m@ is the number of parameters.
-}
type Jacobian r = Vector r → Matrix r


--------------------------------------------------------------------------------
-- Levenberg-Marquardt algorithm.
--------------------------------------------------------------------------------

-- | The Levenberg-Marquardt algorithm is overloaded to work on 'Double' and 'Float'.
class LevMarable r where

    -- | The Levenberg-Marquardt algorithm.
    --
    -- Returns a tuple of the found parameters, a structure containing
    -- information about the minimization and the covariance matrix
    -- corresponding to LS solution.
    levmar ∷ Model r            -- ^ Model
           → Maybe (Jacobian r) -- ^ Optional jacobian
           → Vector r           -- ^ Initial parameters
           → Vector r           -- ^ Samples
           → Int                -- ^ Maximum iterations
           → Options r          -- ^ Minimization options
           → Constraints r      -- ^ Constraints
           → Either LevMarError (Vector r, Info r, Matrix r)

instance LevMarable Float where
    levmar = gen_levmar slevmar_der
                        slevmar_dif
                        slevmar_bc_der
                        slevmar_bc_dif
                        slevmar_lec_der
                        slevmar_lec_dif
                        slevmar_blec_der
                        slevmar_blec_dif

instance LevMarable Double where
    levmar = gen_levmar dlevmar_der
                        dlevmar_dif
                        dlevmar_bc_der
                        dlevmar_bc_dif
                        dlevmar_lec_der
                        dlevmar_lec_dif
                        dlevmar_blec_der
                        dlevmar_blec_dif

{-| @gen_levmar@ takes the low-level C functions as arguments and
executes one of them depending on the optional jacobian and constraints.

Preconditions:

@
  length ys >= length ps

     isJust mLowBs && length (fromJust mLowBs) == length ps
  && isJust mUpBs  && length (fromJust mUpBs)  == length ps

  boxConstrained && (all $ zipWith (<=) (fromJust mLowBs) (fromJust mUpBs))
@
-}
gen_levmar ∷ ∀ r. (Storable r, RealFrac r, Element r)
           ⇒ LevMarDer r
           → LevMarDif r
           → LevMarBCDer r
           → LevMarBCDif r
           → LevMarLecDer r
           → LevMarLecDif r
           → LevMarBLecDer r
           → LevMarBLecDif r

           → Model r            -- ^ Model
           → Maybe (Jacobian r) -- ^ Optional jacobian
           → Vector r           -- ^ Initial parameters
           → Vector r           -- ^ Samples
           → Int                -- ^ Maximum iterations
           → Options r          -- ^ Options
           → Constraints r      -- ^ Constraints
           → Either LevMarError (Vector r, Info r, Matrix r)
gen_levmar f_der
           f_dif
           f_bc_der
           f_bc_dif
           f_lec_der
           f_lec_dif
           f_blec_der
           f_blec_dif
           model mJac ps ys itMax opts (Constraints mLowBs mUpBs mWeights mLinC) =
  -- All effects are contained, so we can safely perform:
  unsafePerformIO $ do

    -- We need to pass the initial parameters 'ps' to the C function.
    -- However, we can't just pass a pointer to them because the C function
    -- will modify the parameters during execution which will violate
    -- referential transparanency. Instead we allocate new space
    -- and copy the parameters to it.
    --
    -- Note that, in the end, the array is returned from this function.
    -- This means that the only way to guarantee its finalisation
    -- is to allocate it using a ForeignPtr:
    psFP ← fastMallocForeignPtrArray lenPs
    withForeignPtr psFP $ \psPtr → do
      VS.unsafeWith ps $ \psPtrInp →
        copyArray psPtr psPtrInp lenPs

      -- Retrieve the (read-only) pointer 'ysPtr' to the samples vector 'ys'
      -- so we can pass it to the C function:
      VS.unsafeWith ys $ \ysPtr →

        -- Convert the Options 'opts' to a list and then to an array
        -- so we can pass the (read-only) pointer 'optsPtr' to the C function:
        withArray (optsToList opts) $ \optsPtr →

          -- Allocate space for the info array
          -- so we can pass it to the C function.
          -- Note that, in the end, this array is converted to an Info value
          -- and returned from this function.
          allocaArray c'LM_INFO_SZ $ \infoPtr → do

            -- Allocate space for the covariance matrix
            -- so we can pass it to the C function.
            -- Like the parameters array the matrix
            -- needs to be returned from this function.
            -- So we also allocate it using a ForeignPtr:
            covarFP ← fastMallocForeignPtrArray covarLen
            withForeignPtr covarFP $ \covarPtr →

              -- 'cmodel' is the low-level model function which is converted
              -- to the FunPtr 'modelFunPtr' and passed to the C function.
              -- 'cmodel' will first convert the parameters pointer 'parPtr'
              -- into a Vector after converting it into a ForeignPtr
              -- (without a finalizer).
              -- Then it will apply the high-level 'model' function
              -- to this parameter vector. The resulting vector is then copied
              -- to the output buffer 'hxPtr':
              let cmodel ∷ Bindings.LevMar.Model r
                  cmodel parPtr hxPtr _ _ _ = do
                    parFP ← newForeignPtr_ parPtr
                    let psV = VS.unsafeFromForeignPtr parFP 0 lenPs
                        vector = model psV
                    VS.unsafeWith vector $ \p → copyArray hxPtr p (VS.length vector)
              in withModel cmodel $ \modelFunPtr → do

                 -- All the low-level C functions share a common set of arguments.
                 -- 'runDif' applies these arguments to the given C function 'f':
                 let runDif ∷ LevMarDif r → IO Int
                     runDif f = f modelFunPtr
                                  psPtr
                                  ysPtr
                                  lenPs
                                  lenYs
                                  itMax
                                  optsPtr
                                  infoPtr
                                  nullPtr
                                  covarPtr
                                  nullPtr

                 err ← case mJac of
                   Nothing → if boxConstrained
                             then if linConstrained
                                  then withBoxConstraints
                                         (withLinConstraints $ withWeights runDif)
                                         f_blec_dif
                                  else withBoxConstraints runDif f_bc_dif
                             else if linConstrained
                                  then withLinConstraints runDif f_lec_dif
                                  else runDif f_dif

                   Just jac →
                     let cjacobian ∷ Bindings.LevMar.Jacobian r
                         cjacobian parPtr jPtr _ _ _ = do
                           parFP ← newForeignPtr_ parPtr
                           let psV    = VS.unsafeFromForeignPtr parFP 0 lenPs
                               matrix = jac psV
                               vector = flatten matrix
                           VS.unsafeWith vector $ \p →
                             copyArray jPtr p (VS.length vector)
                     in withJacobian cjacobian $ \jacobPtr →

                       let runDer ∷ LevMarDer r → IO Int
                           runDer f = runDif $ f jacobPtr
                       in if boxConstrained
                          then if linConstrained
                               then withBoxConstraints
                                      (withLinConstraints $ withWeights runDer)
                                      f_blec_der
                               else withBoxConstraints runDer f_bc_der
                          else if linConstrained
                               then withLinConstraints runDer f_lec_der
                               else runDer f_der

                 -- Handling errors:
                 if err < 0
                    -- we don't treat the following two as an error:
                    ∧ err ≢ c'LM_ERROR_SINGULAR_MATRIX
                    ∧ err ≢ c'LM_ERROR_SUM_OF_SQUARES_NOT_FINITE
                   then return $ Left $ convertLevMarError err

                   else do -- Converting results:
                           info ← listToInfo <$> peekArray c'LM_INFO_SZ infoPtr
                           let psV = VS.unsafeFromForeignPtr psFP 0 lenPs
                           let covarM = reshape lenPs $
                                          VS.unsafeFromForeignPtr covarFP 0 covarLen

                           return $ Right (psV, info, covarM)
  where
    lenPs          = VS.length ps
    lenYs          = VS.length ys
    covarLen       = lenPs*lenPs
    (cMat, rhcVec) = fromJust mLinC

    -- Whether the parameters are constrained by a linear equation.
    linConstrained = isJust mLinC

    -- Whether the parameters are constrained by a bounding box.
    boxConstrained = isJust mLowBs ∨ isJust mUpBs

    withBoxConstraints f g =
        maybeWithArray mLowBs $ \lBsPtr →
          maybeWithArray mUpBs $ \uBsPtr →
            f $ g lBsPtr uBsPtr

    withLinConstraints f g =
        VS.unsafeWith (flatten cMat) $ \cMatPtr →
          VS.unsafeWith rhcVec $ \rhcVecPtr →
            f ∘ g cMatPtr rhcVecPtr $ rows cMat

    withWeights f g = maybeWithArray mWeights $ f ∘ g

maybeWithArray ∷ (Storable α) ⇒ Maybe (Vector α) → (Ptr α → IO β) → IO β
maybeWithArray Nothing  f = f nullPtr
maybeWithArray (Just v) f = VS.unsafeWith v f

#if __GLASGOW_HASKELL__ >= 605
{-# INLINE fastMallocForeignPtrArray #-}
fastMallocForeignPtrArray ∷ ∀ α. Storable α ⇒ Int → IO (ForeignPtr α)
fastMallocForeignPtrArray n = mallocPlainForeignPtrBytes
                                (n * sizeOf (undefined ∷ α))
#else
fastMallocForeignPtrArray ∷ ∀ α. Storable α ⇒ Int → IO (ForeignPtr α)
fastMallocForeignPtrArray = mallocForeignPtrArray
#endif


--------------------------------------------------------------------------------
-- Minimization options.
--------------------------------------------------------------------------------

-- | Minimization options
data Options r =
    Opts { optScaleInitMu      ∷ !r -- ^ Scale factor for initial /mu/.
         , optStopNormInfJacTe ∷ !r -- ^ Stopping thresholds for @||J^T e||_inf@.
         , optStopNorm2Dp      ∷ !r -- ^ Stopping thresholds for @||Dp||_2@.
         , optStopNorm2E       ∷ !r -- ^ Stopping thresholds for @||e||_2@.
         , optDelta            ∷ !r -- ^ Step used in the difference
                                    -- approximation to the Jacobian. If
                                    -- @optDelta<0@, the Jacobian is approximated
                                    -- with central differences which are more
                                    -- accurate (but slower!)  compared to the
                                    -- forward differences employed by default.
         } deriving (Read, Show)

-- | Default minimization options
defaultOpts ∷ Fractional r ⇒ Options r
defaultOpts = Opts { optScaleInitMu      = c'LM_INIT_MU
                   , optStopNormInfJacTe = c'LM_STOP_THRESH
                   , optStopNorm2Dp      = c'LM_STOP_THRESH
                   , optStopNorm2E       = c'LM_STOP_THRESH
                   , optDelta            = c'LM_DIFF_DELTA
                   }

optsToList ∷ Options r → [r]
optsToList (Opts mu  eps1  eps2  eps3  delta) =
                [mu, eps1, eps2, eps3, delta]


--------------------------------------------------------------------------------
-- Constraints
--------------------------------------------------------------------------------

-- | Ensure that these vectors have the same length as the number of parameters.
data Constraints r = Constraints
    { lowerBounds       ∷ !(Maybe (Vector r))            -- ^ Optional lower bounds
    , upperBounds       ∷ !(Maybe (Vector r))            -- ^ Optional upper bounds
    , weights           ∷ !(Maybe (Vector r))            -- ^ Optional weights
    , linearConstraints ∷ !(Maybe (LinearConstraints r)) -- ^ Optional linear constraints
    }

-- | Linear constraints consisting of a constraints matrix, @k><m@ and
--   a right hand constraints vector, of length @k@ where @m@ is the number of
--   parameters and @k@ is the number of constraints.
type LinearConstraints r = (Matrix r, Vector r)

-- | * 'mempty' is defined as a 'Constraints' where all fields are 'Nothing'.
--
--   * 'mappend' merges two 'Constraints' by taking the first non-'Nothing' value
--     for each field.
instance Monoid (Constraints r) where
    mempty = Constraints Nothing Nothing Nothing Nothing
    mappend (Constraints lb1 ub1 w1 l1)
            (Constraints lb2 ub2 w2 l2) = Constraints (lb1 `mplus` lb2)
                                                      (ub1 `mplus` ub2)
                                                      (w1  `mplus` w2)
                                                      (l1  `mplus` l2)


--------------------------------------------------------------------------------
-- Output
--------------------------------------------------------------------------------

-- | Information regarding the minimization.
data Info r = Info
  { infNorm2initE      ∷ !r          -- ^ @||e||_2@             at initial parameters.
  , infNorm2E          ∷ !r          -- ^ @||e||_2@             at estimated parameters.
  , infNormInfJacTe    ∷ !r          -- ^ @||J^T e||_inf@       at estimated parameters.
  , infNorm2Dp         ∷ !r          -- ^ @||Dp||_2@            at estimated parameters.
  , infMuDivMax        ∷ !r          -- ^ @\mu/max[J^T J]_ii ]@ at estimated parameters.
  , infNumIter         ∷ !Int        -- ^ Number of iterations.
  , infStopReason      ∷ !StopReason -- ^ Reason for terminating.
  , infNumFuncEvals    ∷ !Int        -- ^ Number of function evaluations.
  , infNumJacobEvals   ∷ !Int        -- ^ Number of jacobian evaluations.
  , infNumLinSysSolved ∷ !Int        -- ^ Number of linear systems solved,
                                     --   i.e. attempts for reducing error.
  } deriving (Read, Show)

listToInfo ∷ (RealFrac r) ⇒ [r] → Info r
listToInfo [a,b,c,d,e,f,g,h,i,j] =
    Info { infNorm2initE      = a
         , infNorm2E          = b
         , infNormInfJacTe    = c
         , infNorm2Dp         = d
         , infMuDivMax        = e
         , infNumIter         = floor f
         , infStopReason      = toEnum $ floor g - 1
         , infNumFuncEvals    = floor h
         , infNumJacobEvals   = floor i
         , infNumLinSysSolved = floor j
         }
listToInfo _ = error "liftToInfo: wrong list length"

-- | Reason for terminating.
data StopReason
  = SmallGradient  -- ^ Stopped because of small gradient @J^T e@.
  | SmallDp        -- ^ Stopped because of small Dp.
  | MaxIterations  -- ^ Stopped because maximum iterations was reached.
  | SingularMatrix -- ^ Stopped because of singular matrix. Restart from current
                   --   estimated parameters with increased 'optScaleInitMu'.
  | SmallestError  -- ^ Stopped because no further error reduction is
                   --   possible. Restart with increased 'optScaleInitMu'.
  | SmallNorm2E    -- ^ Stopped because of small @||e||_2@.
  | InvalidValues  -- ^ Stopped because model function returned invalid values
                   --   (i.e. NaN or Inf). This is a user error.
    deriving (Read, Show, Enum)


--------------------------------------------------------------------------------
-- Error
--------------------------------------------------------------------------------

data LevMarError
    = LevMarError                    -- ^ Generic error (not one of the others)
    | LapackError                    -- ^ A call to a lapack subroutine failed
                                     --   in the underlying C levmar library.
    | FailedBoxCheck                 -- ^ At least one lower bound exceeds the
                                     --   upper one.
    | MemoryAllocationFailure        -- ^ A call to @malloc@ failed in the
                                     --   underlying C levmar library.
    | ConstraintMatrixRowsGtCols     -- ^ The matrix of constraints cannot have
                                     --   more rows than columns.
    | ConstraintMatrixNotFullRowRank -- ^ Constraints matrix is not of full row
                                     --   rank.
    | TooFewMeasurements             -- ^ Cannot solve a problem with fewer
                                     --   measurements than unknowns.  In case
                                     --   linear constraints are provided, this
                                     --   error is also returned when the number
                                     --   of measurements is smaller than the
                                     --   number of unknowns minus the number of
                                     --   equality constraints.
      deriving (Show, Typeable)

-- Handy in case you want to thow a LevMarError as an exception:
instance Exception LevMarError

levmarCErrorToLevMarError ∷ [(Int, LevMarError)]
levmarCErrorToLevMarError =
    [ (c'LM_ERROR,                                     LevMarError)
    , (c'LM_ERROR_LAPACK_ERROR,                        LapackError)
  --, (c'LM_ERROR_NO_JACOBIAN,                         can never happen)
  --, (c'LM_ERROR_NO_BOX_CONSTRAINTS,                  can never happen)
    , (c'LM_ERROR_FAILED_BOX_CHECK,                    FailedBoxCheck)
    , (c'LM_ERROR_MEMORY_ALLOCATION_FAILURE,           MemoryAllocationFailure)
    , (c'LM_ERROR_CONSTRAINT_MATRIX_ROWS_GT_COLS,      ConstraintMatrixRowsGtCols)
    , (c'LM_ERROR_CONSTRAINT_MATRIX_NOT_FULL_ROW_RANK, ConstraintMatrixNotFullRowRank)
    , (c'LM_ERROR_TOO_FEW_MEASUREMENTS,                TooFewMeasurements)
  --, (c'LM_ERROR_SINGULAR_MATRIX,                     we don't treat this as an error)
  --, (c'LM_ERROR_SUM_OF_SQUARES_NOT_FINITE,           we don't treat this as an error)
    ]

convertLevMarError ∷ Int → LevMarError
convertLevMarError err = fromMaybe (error $ "Unknown levmar error: " ++ show err)
                                   (lookup err levmarCErrorToLevMarError)


-- The End ---------------------------------------------------------------------
