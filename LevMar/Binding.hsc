{-# LANGUAGE ForeignFunctionInterface #-}

module LevMar.Binding
    ( _LM_OPTS_SZ
    , _LM_INFO_SZ

    , _LM_ERROR

    , _LM_INIT_MU
    , _LM_STOP_THRESH
    , _LM_DIFF_DELTA

    , _LM_VERSION

    , Model
    , Jacobian

    , withModel
    , withJacobian

    , LevMarDer
    , LevMarDif
    , LevMarBCDer
    , LevMarBCDif
    , LevMarLecDer
    , LevMarLecDif
    , LevMarBLecDer
    , LevMarBLecDif

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

    , LevMarDer'
    , LevMarDif'
    , LevMarBCDer'
    , LevMarBCDif'
    , LevMarLecDer'
    , LevMarLecDif'
    , LevMarBLecDer'
    , LevMarBLecDif'

    , dlevmar_der'
    , slevmar_der'
    , dlevmar_dif'
    , slevmar_dif'
    , dlevmar_bc_der'
    , slevmar_bc_der'
    , dlevmar_bc_dif'
    , slevmar_bc_dif'
    , dlevmar_lec_der'
    , slevmar_lec_der'
    , dlevmar_lec_dif'
    , slevmar_lec_dif'
    , dlevmar_blec_der'
    , slevmar_blec_der'
    , dlevmar_blec_dif'
    , slevmar_blec_dif'
    ) where

import Foreign.C.Types   (CInt, CFloat, CDouble)
import Foreign.Ptr       (Ptr, FunPtr, freeHaskellFunPtr)
import Control.Exception (bracket)

#include <lm.h>

_LM_OPTS_SZ, _LM_INFO_SZ :: Int
_LM_OPTS_SZ = #const LM_OPTS_SZ
_LM_INFO_SZ = #const LM_INFO_SZ

_LM_ERROR :: CInt
_LM_ERROR = #const LM_ERROR

-- TODO: The HSC #const construct only works for integer values. We
-- need something similar for floating-point values. It is posible to
-- define your own constructs but we haven't done that yet:
_LM_INIT_MU, _LM_STOP_THRESH, _LM_DIFF_DELTA :: Fractional a => a
_LM_INIT_MU     = #const_real LM_INIT_MU
_LM_STOP_THRESH = #const_real LM_STOP_THRESH
_LM_DIFF_DELTA  = #const_real LM_DIFF_DELTA

_LM_VERSION :: String
_LM_VERSION = #const_str LM_VERSION

type Model r =  Ptr r  -- p
             -> Ptr r  -- hx
             -> CInt   -- m
             -> CInt   -- n
             -> Ptr () -- adata
             -> IO ()

type Jacobian a = Model a

foreign import ccall "wrapper"
  mkModel :: Model a -> IO (FunPtr (Model a))

mkJacobian :: Jacobian a -> IO (FunPtr (Jacobian a))
mkJacobian = mkModel

withModel :: Model a -> (FunPtr (Model a) -> IO b) -> IO b
withModel m = bracket (mkModel m) freeHaskellFunPtr

withJacobian :: Jacobian a -> (FunPtr (Jacobian a) -> IO b) -> IO b
withJacobian j = bracket (mkJacobian j) freeHaskellFunPtr

type LevMarDer cr =  FunPtr (Model cr)    -- func
                  -> FunPtr (Jacobian cr) -- jacf
                  -> Ptr cr               -- p
                  -> Ptr cr               -- x
                  -> CInt                 -- m
                  -> CInt                 -- n
                  -> CInt                 -- itmax
                  -> Ptr cr               -- opts
                  -> Ptr cr               -- info
                  -> Ptr cr               -- work
                  -> Ptr cr               -- covar
                  -> Ptr ()               -- adata
                  -> IO CInt

type LevMarDif cr =  FunPtr (Model cr) -- func
                  -> Ptr cr            -- p
                  -> Ptr cr            -- x
                  -> CInt              -- m
                  -> CInt              -- n
                  -> CInt              -- itmax
                  -> Ptr cr            -- opts
                  -> Ptr cr            -- info
                  -> Ptr cr            -- work
                  -> Ptr cr            -- covar
                  -> Ptr ()            -- adata
                  -> IO CInt

type LevMarBCDer cr =  FunPtr (Model cr)    -- func
                    -> FunPtr (Jacobian cr) -- jacf
                    -> Ptr cr               -- p
                    -> Ptr cr               -- x
                    -> CInt                 -- m
                    -> CInt                 -- n
                    -> Ptr cr               -- lb
                    -> Ptr cr               -- ub
                    -> CInt                 -- itmax
                    -> Ptr cr               -- opts
                    -> Ptr cr               -- info
                    -> Ptr cr               -- work
                    -> Ptr cr               -- covar
                    -> Ptr ()               -- adata
                    -> IO CInt

type LevMarBCDif cr =  FunPtr (Model cr) -- func
                    -> Ptr cr            -- p
                    -> Ptr cr            -- x
                    -> CInt              -- m
                    -> CInt              -- n
                    -> Ptr cr            -- lb
                    -> Ptr cr            -- ub
                    -> CInt              -- itmax
                    -> Ptr cr            -- opts
                    -> Ptr cr            -- info
                    -> Ptr cr            -- work
                    -> Ptr cr            -- covar
                    -> Ptr ()            -- adata
                    -> IO CInt

type LevMarLecDer cr =  FunPtr (Model cr)    -- func
                     -> FunPtr (Jacobian cr) -- jacf
                     -> Ptr cr               -- p
                     -> Ptr cr               -- x
                     -> CInt                 -- m
                     -> CInt                 -- n
                     -> Ptr cr               -- A
                     -> Ptr cr               -- B
                     -> CInt                 -- k
                     -> CInt                 -- itmax
                     -> Ptr cr               -- opts
                     -> Ptr cr               -- info
                     -> Ptr cr               -- work
                     -> Ptr cr               -- covar
                     -> Ptr ()               -- adata
                     -> IO CInt

type LevMarLecDif cr =  FunPtr (Model cr) -- func
                     -> Ptr cr            -- p
                     -> Ptr cr            -- x
                     -> CInt              -- m
                     -> CInt              -- n
                     -> Ptr cr            -- A
                     -> Ptr cr            -- B
                     -> CInt              -- k
                     -> CInt              -- itmax
                     -> Ptr cr            -- opts
                     -> Ptr cr            -- info
                     -> Ptr cr            -- work
                     -> Ptr cr            -- covar
                     -> Ptr ()            -- adata
                     -> IO CInt

type LevMarBLecDer cr =  FunPtr (Model cr)    -- func
                      -> FunPtr (Jacobian cr) -- jacf
                      -> Ptr cr               -- p
                      -> Ptr cr               -- x
                      -> CInt                 -- m
                      -> CInt                 -- n
                      -> Ptr cr               -- lb
                      -> Ptr cr               -- ub
                      -> Ptr cr               -- A
                      -> Ptr cr               -- B
                      -> CInt                 -- k
                      -> Ptr cr               -- wghts
                      -> CInt                 -- itmax
                      -> Ptr cr               -- opts
                      -> Ptr cr               -- info
                      -> Ptr cr               -- work
                      -> Ptr cr               -- covar
                      -> Ptr ()               -- adata
                      -> IO CInt

type LevMarBLecDif cr =  FunPtr (Model cr) -- func
                      -> Ptr cr            -- p
                      -> Ptr cr            -- x
                      -> CInt              -- m
                      -> CInt              -- n
                      -> Ptr cr            -- lb
                      -> Ptr cr            -- ub
                      -> Ptr cr            -- A
                      -> Ptr cr            -- B
                      -> CInt              -- k
                      -> Ptr cr            -- wghts
                      -> CInt              -- itmax
                      -> Ptr cr            -- opts
                      -> Ptr cr            -- info
                      -> Ptr cr            -- work
                      -> Ptr cr            -- covar
                      -> Ptr ()            -- adata
                      -> IO CInt

foreign import ccall "slevmar_der"      slevmar_der      :: LevMarDer     CFloat
foreign import ccall "dlevmar_der"      dlevmar_der      :: LevMarDer     CDouble
foreign import ccall "slevmar_dif"      slevmar_dif      :: LevMarDif     CFloat
foreign import ccall "dlevmar_dif"      dlevmar_dif      :: LevMarDif     CDouble
foreign import ccall "slevmar_bc_der"   slevmar_bc_der   :: LevMarBCDer   CFloat
foreign import ccall "dlevmar_bc_der"   dlevmar_bc_der   :: LevMarBCDer   CDouble
foreign import ccall "slevmar_bc_dif"   slevmar_bc_dif   :: LevMarBCDif   CFloat
foreign import ccall "dlevmar_bc_dif"   dlevmar_bc_dif   :: LevMarBCDif   CDouble
foreign import ccall "slevmar_lec_der"  slevmar_lec_der  :: LevMarLecDer  CFloat
foreign import ccall "dlevmar_lec_der"  dlevmar_lec_der  :: LevMarLecDer  CDouble
foreign import ccall "slevmar_lec_dif"  slevmar_lec_dif  :: LevMarLecDif  CFloat
foreign import ccall "dlevmar_lec_dif"  dlevmar_lec_dif  :: LevMarLecDif  CDouble
foreign import ccall "slevmar_blec_der" slevmar_blec_der :: LevMarBLecDer CFloat
foreign import ccall "dlevmar_blec_der" dlevmar_blec_der :: LevMarBLecDer CDouble
foreign import ccall "slevmar_blec_dif" slevmar_blec_dif :: LevMarBLecDif CFloat
foreign import ccall "dlevmar_blec_dif" dlevmar_blec_dif :: LevMarBLecDif CDouble

-------------------------------------------------------------------------------
-- Curry friendly variants
-------------------------------------------------------------------------------

type BoxConstraints    cr a = Ptr cr -> Ptr cr -> a
type LinearConstraints cr a = Ptr cr -> Ptr cr -> CInt -> a

type LevMarDif'     cr = LevMarDif cr
type LevMarDer'     cr = FunPtr (Jacobian cr) -> LevMarDif' cr
type LevMarBCDif'   cr = BoxConstraints cr (LevMarDif' cr)
type LevMarBCDer'   cr = BoxConstraints cr (LevMarDer' cr)
type LevMarLecDif'  cr = LinearConstraints cr (LevMarDif' cr)
type LevMarLecDer'  cr = LinearConstraints cr (LevMarDer' cr)
type LevMarBLecDif' cr = BoxConstraints cr (LinearConstraints cr (Ptr cr -> LevMarDif' cr))
type LevMarBLecDer' cr = BoxConstraints cr (LinearConstraints cr (Ptr cr -> LevMarDer' cr))


mk_levmar_der' :: LevMarDer cr -> LevMarDer' cr
mk_levmar_der' lma j f -- p x m n itMax opts info work covar aData
             = lma f j -- p x m n itMax opts info work covar aData

mk_levmar_bc_dif' :: LevMarBCDif cr -> LevMarBCDif' cr
mk_levmar_bc_dif' lma lb ub f p x m n -- itMax opts info work covar aData
                = lma f p x m n lb ub -- itMax opts info work covar aData

mk_levmar_bc_der' :: LevMarBCDer cr -> LevMarBCDer' cr
mk_levmar_bc_der' lma lb ub j f p x m n -- itMax opts info work covar aData
                = lma f j p x m n lb ub -- itMax opts info work covar aData

mk_levmar_lec_dif' :: LevMarLecDif cr -> LevMarLecDif' cr
mk_levmar_lec_dif' lma a b k f p x m n -- itMax opts info work covar aData
                 = lma f p x m n a b k -- itMax opts info work covar aData

mk_levmar_lec_der' :: LevMarLecDer cr -> LevMarLecDer' cr
mk_levmar_lec_der' lma a b k j f p x m n -- itMax opts info work covar aData
                 = lma f j p x m n a b k -- itMax opts info work covar aData

mk_levmar_blec_dif' :: LevMarBLecDif cr -> LevMarBLecDif' cr
mk_levmar_blec_dif' lma lb ub a b k wghts f p x m n -- itMax opts info work covar aData
                  = lma f p x m n lb ub a b k wghts -- itMax opts info work covar aData

mk_levmar_blec_der' :: LevMarBLecDer cr -> LevMarBLecDer' cr
mk_levmar_blec_der' lma lb ub a b k wghts j f p x m n -- itMax opts info work covar aData
                  = lma f j p x m n lb ub a b k wghts -- itMax opts info work covar aData


slevmar_dif' :: LevMarDif' CFloat
slevmar_dif' = slevmar_dif

dlevmar_dif' :: LevMarDif' CDouble
dlevmar_dif' = dlevmar_dif

slevmar_der' :: LevMarDer' CFloat
slevmar_der' = mk_levmar_der' slevmar_der

dlevmar_der' :: LevMarDer' CDouble
dlevmar_der' = mk_levmar_der' dlevmar_der

slevmar_bc_dif' :: LevMarBCDif' CFloat
slevmar_bc_dif' = mk_levmar_bc_dif' slevmar_bc_dif

dlevmar_bc_dif' :: LevMarBCDif' CDouble
dlevmar_bc_dif' = mk_levmar_bc_dif' dlevmar_bc_dif

slevmar_bc_der' :: LevMarBCDer' CFloat
slevmar_bc_der' = mk_levmar_bc_der' slevmar_bc_der

dlevmar_bc_der' :: LevMarBCDer' CDouble
dlevmar_bc_der' = mk_levmar_bc_der' dlevmar_bc_der

slevmar_lec_dif' :: LevMarLecDif' CFloat
slevmar_lec_dif' = mk_levmar_lec_dif' slevmar_lec_dif

dlevmar_lec_dif' :: LevMarLecDif' CDouble
dlevmar_lec_dif' = mk_levmar_lec_dif' dlevmar_lec_dif

slevmar_lec_der' :: LevMarLecDer' CFloat
slevmar_lec_der' = mk_levmar_lec_der' slevmar_lec_der

dlevmar_lec_der' :: LevMarLecDer' CDouble
dlevmar_lec_der' = mk_levmar_lec_der' dlevmar_lec_der

slevmar_blec_dif' :: LevMarBLecDif' CFloat
slevmar_blec_dif' = mk_levmar_blec_dif' slevmar_blec_dif

dlevmar_blec_dif' :: LevMarBLecDif' CDouble
dlevmar_blec_dif' = mk_levmar_blec_dif' dlevmar_blec_dif

slevmar_blec_der' :: LevMarBLecDer' CFloat
slevmar_blec_der' = mk_levmar_blec_der' slevmar_blec_der

dlevmar_blec_der' :: LevMarBLecDer' CDouble
dlevmar_blec_der' = mk_levmar_blec_der' dlevmar_blec_der
