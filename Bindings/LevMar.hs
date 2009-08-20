{-# LANGUAGE ForeignFunctionInterface #-}

module Bindings.LevMar
    ( Model
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
    ) where

import Foreign.C.Types   (CInt, CFloat, CDouble)
import Foreign.Ptr       (Ptr, FunPtr, freeHaskellFunPtr)
import Control.Exception (bracket)

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

foreign import ccall "dlevmar_der"      dlevmar_der      :: LevMarDer     CDouble
foreign import ccall "slevmar_der"      slevmar_der      :: LevMarDer     CFloat
foreign import ccall "dlevmar_dif"      dlevmar_dif      :: LevMarDif     CDouble
foreign import ccall "slevmar_dif"      slevmar_dif      :: LevMarDif     CFloat
foreign import ccall "dlevmar_bc_der"   dlevmar_bc_der   :: LevMarBCDer   CDouble
foreign import ccall "slevmar_bc_der"   slevmar_bc_der   :: LevMarBCDer   CFloat
foreign import ccall "dlevmar_bc_dif"   dlevmar_bc_dif   :: LevMarBCDif   CDouble
foreign import ccall "slevmar_bc_dif"   slevmar_bc_dif   :: LevMarBCDif   CFloat
foreign import ccall "dlevmar_lec_der"  dlevmar_lec_der  :: LevMarLecDer  CDouble
foreign import ccall "slevmar_lec_der"  slevmar_lec_der  :: LevMarLecDer  CFloat
foreign import ccall "dlevmar_lec_dif"  dlevmar_lec_dif  :: LevMarLecDif  CDouble
foreign import ccall "slevmar_lec_dif"  slevmar_lec_dif  :: LevMarLecDif  CFloat
foreign import ccall "dlevmar_blec_der" dlevmar_blec_der :: LevMarBLecDer CDouble
foreign import ccall "slevmar_blec_der" slevmar_blec_der :: LevMarBLecDer CFloat
foreign import ccall "dlevmar_blec_dif" dlevmar_blec_dif :: LevMarBLecDif CDouble
foreign import ccall "slevmar_blec_dif" slevmar_blec_dif :: LevMarBLecDif CFloat
