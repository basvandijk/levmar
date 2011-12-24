{-# LANGUAGE NoImplicitPrelude, UnicodeSyntax #-}

module Bindings.LevMar.CurryFriendly
    ( -- * Handy type synonyms used in the curry friendly types.
      BoxConstraints
    , LinearConstraints

      -- * Curry friendly types of the Levenberg-Marquardt algorithms.
    , LevMarDer
    , LevMarDif
    , LevMarBCDer
    , LevMarBCDif
    , LevMarLecDer
    , LevMarLecDif
    , LevMarBLecDer
    , LevMarBLecDif

      -- * Curry friendly variants of the Levenberg-Marquardt
      -- algorithms in 'Bindings.Levmar'.
    , dlevmar_der,      slevmar_der
    , dlevmar_dif,      slevmar_dif
    , dlevmar_bc_der,   slevmar_bc_der
    , dlevmar_bc_dif,   slevmar_bc_dif
    , dlevmar_lec_der,  slevmar_lec_der
    , dlevmar_lec_dif,  slevmar_lec_dif
    , dlevmar_blec_der, slevmar_blec_der
    , dlevmar_blec_dif, slevmar_blec_dif
    ) where

import Prelude     ( Double, Float )
import Foreign.Ptr ( FunPtr )

import qualified Bindings.LevMar as BLM


--------------------------------------------------------------------------------
-- Handy type synonyms used in the curry friendly types.
--------------------------------------------------------------------------------

type BoxConstraints r α = BLM.LowerBounds r
                        → BLM.UpperBounds r
                        → α

type LinearConstraints r α = BLM.ConstraintsMatrix r
                           → BLM.ConstraintsVector r
                           → BLM.NrOfConstraints
                           → α


--------------------------------------------------------------------------------
-- Curry friendly types of the Levenberg-Marquardt algorithms.
--------------------------------------------------------------------------------

type LevMarDif     r = BLM.LevMarDif r
type LevMarDer     r = FunPtr (BLM.Jacobian r) → LevMarDif r
type LevMarBCDif   r = BoxConstraints r (LevMarDif r)
type LevMarBCDer   r = BoxConstraints r (LevMarDer r)
type LevMarLecDif  r = LinearConstraints r (LevMarDif r)
type LevMarLecDer  r = LinearConstraints r (LevMarDer r)
type LevMarBLecDif r = BoxConstraints r (LinearConstraints r (BLM.Weights r → LevMarDif r))
type LevMarBLecDer r = BoxConstraints r (LinearConstraints r (BLM.Weights r → LevMarDer r))


--------------------------------------------------------------------------------
-- Reordering arguments to create curry friendly variants.
--------------------------------------------------------------------------------

mk_levmar_der ∷ BLM.LevMarDer r → LevMarDer r
mk_levmar_der lma j f
            = lma f j

mk_levmar_bc_dif ∷ BLM.LevMarBCDif r → LevMarBCDif r
mk_levmar_bc_dif lma lb ub f p x m n
               = lma f p x m n lb ub

mk_levmar_bc_der ∷ BLM.LevMarBCDer r → LevMarBCDer r
mk_levmar_bc_der lma lb ub j f p x m n
               = lma f j p x m n lb ub

mk_levmar_lec_dif ∷ BLM.LevMarLecDif r → LevMarLecDif r
mk_levmar_lec_dif lma a b k f p x m n
                = lma f p x m n a b k

mk_levmar_lec_der ∷ BLM.LevMarLecDer r → LevMarLecDer r
mk_levmar_lec_der lma a b k j f p x m n
                = lma f j p x m n a b k

mk_levmar_blec_dif ∷ BLM.LevMarBLecDif r → LevMarBLecDif r
mk_levmar_blec_dif lma lb ub a b k wghts f p x m n
                 = lma f p x m n lb ub a b k wghts

mk_levmar_blec_der ∷ BLM.LevMarBLecDer r → LevMarBLecDer r
mk_levmar_blec_der lma lb ub a b k wghts j f p x m n
                 = lma f j p x m n lb ub a b k wghts


--------------------------------------------------------------------------------
-- Curry friendly variants of the Levenberg-Marquardt algorithms in
-- 'Bindings.Levmar'.
--------------------------------------------------------------------------------

slevmar_dif ∷ LevMarDif Float
slevmar_dif = BLM.c'slevmar_dif

dlevmar_dif ∷ LevMarDif Double
dlevmar_dif = BLM.c'dlevmar_dif

slevmar_der ∷ LevMarDer Float
slevmar_der = mk_levmar_der BLM.c'slevmar_der

dlevmar_der ∷ LevMarDer Double
dlevmar_der = mk_levmar_der BLM.c'dlevmar_der

slevmar_bc_dif ∷ LevMarBCDif Float
slevmar_bc_dif = mk_levmar_bc_dif BLM.c'slevmar_bc_dif

dlevmar_bc_dif ∷ LevMarBCDif Double
dlevmar_bc_dif = mk_levmar_bc_dif BLM.c'dlevmar_bc_dif

slevmar_bc_der ∷ LevMarBCDer Float
slevmar_bc_der = mk_levmar_bc_der BLM.c'slevmar_bc_der

dlevmar_bc_der ∷ LevMarBCDer Double
dlevmar_bc_der = mk_levmar_bc_der BLM.c'dlevmar_bc_der

slevmar_lec_dif ∷ LevMarLecDif Float
slevmar_lec_dif = mk_levmar_lec_dif BLM.c'slevmar_lec_dif

dlevmar_lec_dif ∷ LevMarLecDif Double
dlevmar_lec_dif = mk_levmar_lec_dif BLM.c'dlevmar_lec_dif

slevmar_lec_der ∷ LevMarLecDer Float
slevmar_lec_der = mk_levmar_lec_der BLM.c'slevmar_lec_der

dlevmar_lec_der ∷ LevMarLecDer Double
dlevmar_lec_der = mk_levmar_lec_der BLM.c'dlevmar_lec_der

slevmar_blec_dif ∷ LevMarBLecDif Float
slevmar_blec_dif = mk_levmar_blec_dif BLM.c'slevmar_blec_dif

dlevmar_blec_dif ∷ LevMarBLecDif Double
dlevmar_blec_dif = mk_levmar_blec_dif BLM.c'dlevmar_blec_dif

slevmar_blec_der ∷ LevMarBLecDer Float
slevmar_blec_der = mk_levmar_blec_der BLM.c'slevmar_blec_der

dlevmar_blec_der ∷ LevMarBLecDer Double
dlevmar_blec_der = mk_levmar_blec_der BLM.c'dlevmar_blec_der
