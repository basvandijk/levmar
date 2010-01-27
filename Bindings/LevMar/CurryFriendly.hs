{-# LANGUAGE UnicodeSyntax #-}

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
    , dlevmar_der
    , slevmar_der
    , dlevmar_dif
    , slevmar_dif
    , dlevmar_bc_der
    , slevmar_bc_der
    , dlevmar_bc_dif
    , slevmar_bc_dif
    , dlevmar_lec_der
    , slevmar_lec_der
    , dlevmar_lec_dif
    , slevmar_lec_dif
    , dlevmar_blec_der
    , slevmar_blec_der
    , dlevmar_blec_dif
    , slevmar_blec_dif
    ) where


import Foreign.C.Types ( CFloat, CDouble )
import Foreign.Ptr     ( FunPtr )

import qualified Bindings.LevMar as LMA_C


--------------------------------------------------------------------------------
-- Handy type synonyms used in the curry friendly types.
--------------------------------------------------------------------------------

type BoxConstraints cr a = LMA_C.LowerBounds cr
                         → LMA_C.UpperBounds cr
                         → a

type LinearConstraints cr a = LMA_C.ConstraintsMatrix cr
                            → LMA_C.ConstraintsVector cr
                            → LMA_C.NrOfConstraints
                            → a


--------------------------------------------------------------------------------
-- Curry friendly types of the Levenberg-Marquardt algorithms.
--------------------------------------------------------------------------------

type LevMarDif     cr = LMA_C.LevMarDif cr
type LevMarDer     cr = FunPtr (LMA_C.Jacobian cr) → LevMarDif cr
type LevMarBCDif   cr = BoxConstraints cr (LevMarDif cr)
type LevMarBCDer   cr = BoxConstraints cr (LevMarDer cr)
type LevMarLecDif  cr = LinearConstraints cr (LevMarDif cr)
type LevMarLecDer  cr = LinearConstraints cr (LevMarDer cr)
type LevMarBLecDif cr = BoxConstraints cr (LinearConstraints cr (LMA_C.Weights cr → (LevMarDif cr)))
type LevMarBLecDer cr = BoxConstraints cr (LinearConstraints cr (LMA_C.Weights cr → (LevMarDer cr)))


--------------------------------------------------------------------------------
-- Reordering arguments to create curry friendly variants.
--------------------------------------------------------------------------------

mk_levmar_der ∷ LMA_C.LevMarDer cr → LevMarDer cr
mk_levmar_der lma j f
            = lma f j

mk_levmar_bc_dif ∷ LMA_C.LevMarBCDif cr → LevMarBCDif cr
mk_levmar_bc_dif lma lb ub f p x m n
               = lma f p x m n lb ub

mk_levmar_bc_der ∷ LMA_C.LevMarBCDer cr → LevMarBCDer cr
mk_levmar_bc_der lma lb ub j f p x m n
               = lma f j p x m n lb ub

mk_levmar_lec_dif ∷ LMA_C.LevMarLecDif cr → LevMarLecDif cr
mk_levmar_lec_dif lma a b k f p x m n
                = lma f p x m n a b k

mk_levmar_lec_der ∷ LMA_C.LevMarLecDer cr → LevMarLecDer cr
mk_levmar_lec_der lma a b k j f p x m n
                = lma f j p x m n a b k

mk_levmar_blec_dif ∷ LMA_C.LevMarBLecDif cr → LevMarBLecDif cr
mk_levmar_blec_dif lma lb ub a b k wghts f p x m n
                 = lma f p x m n lb ub a b k wghts

mk_levmar_blec_der ∷ LMA_C.LevMarBLecDer cr → LevMarBLecDer cr
mk_levmar_blec_der lma lb ub a b k wghts j f p x m n
                 = lma f j p x m n lb ub a b k wghts


--------------------------------------------------------------------------------
-- Curry friendly variants of the Levenberg-Marquardt algorithms in
-- 'Bindings.Levmar'.
--------------------------------------------------------------------------------

slevmar_dif ∷ LevMarDif CFloat
slevmar_dif = LMA_C.c'slevmar_dif

dlevmar_dif ∷ LevMarDif CDouble
dlevmar_dif = LMA_C.c'dlevmar_dif

slevmar_der ∷ LevMarDer CFloat
slevmar_der = mk_levmar_der LMA_C.c'slevmar_der

dlevmar_der ∷ LevMarDer CDouble
dlevmar_der = mk_levmar_der LMA_C.c'dlevmar_der

slevmar_bc_dif ∷ LevMarBCDif CFloat
slevmar_bc_dif = mk_levmar_bc_dif LMA_C.c'slevmar_bc_dif

dlevmar_bc_dif ∷ LevMarBCDif CDouble
dlevmar_bc_dif = mk_levmar_bc_dif LMA_C.c'dlevmar_bc_dif

slevmar_bc_der ∷ LevMarBCDer CFloat
slevmar_bc_der = mk_levmar_bc_der LMA_C.c'slevmar_bc_der

dlevmar_bc_der ∷ LevMarBCDer CDouble
dlevmar_bc_der = mk_levmar_bc_der LMA_C.c'dlevmar_bc_der

slevmar_lec_dif ∷ LevMarLecDif CFloat
slevmar_lec_dif = mk_levmar_lec_dif LMA_C.c'slevmar_lec_dif

dlevmar_lec_dif ∷ LevMarLecDif CDouble
dlevmar_lec_dif = mk_levmar_lec_dif LMA_C.c'dlevmar_lec_dif

slevmar_lec_der ∷ LevMarLecDer CFloat
slevmar_lec_der = mk_levmar_lec_der LMA_C.c'slevmar_lec_der

dlevmar_lec_der ∷ LevMarLecDer CDouble
dlevmar_lec_der = mk_levmar_lec_der LMA_C.c'dlevmar_lec_der

slevmar_blec_dif ∷ LevMarBLecDif CFloat
slevmar_blec_dif = mk_levmar_blec_dif LMA_C.c'slevmar_blec_dif

dlevmar_blec_dif ∷ LevMarBLecDif CDouble
dlevmar_blec_dif = mk_levmar_blec_dif LMA_C.c'dlevmar_blec_dif

slevmar_blec_der ∷ LevMarBLecDer CFloat
slevmar_blec_der = mk_levmar_blec_der LMA_C.c'slevmar_blec_der

dlevmar_blec_der ∷ LevMarBLecDer CDouble
dlevmar_blec_der = mk_levmar_blec_der LMA_C.c'dlevmar_blec_der


-- The End ---------------------------------------------------------------------
