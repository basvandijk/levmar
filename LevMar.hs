{-# LANGUAGE ForeignFunctionInterface #-}

module LevMar ( Model
              , Jacobian

              , mkModel
              , mkJacobian

              , c_dlevmar_der
              , c_dlevmar_dif
              , c_dlevmar_bc_der
              , c_dlevmar_bc_dif
              , c_dlevmar_lec_der
              , c_dlevmar_lec_dif
              , c_dlevmar_blec_der
              , c_dlevmar_blec_dif

              , c_slevmar_der
              , c_slevmar_dif
              , c_slevmar_bc_der
              , c_slevmar_bc_dif
              , c_slevmar_lec_der
              , c_slevmar_lec_dif
              , c_slevmar_blec_der
              , c_slevmar_blec_dif
              ) where

import Foreign.C.Types
import Foreign.Ptr

type Model a =  Ptr a  -- p
             -> Ptr a  -- hx
             -> CInt   -- m
             -> CInt   -- n
             -> Ptr () -- adata
             -> IO ()

type Jacobian a = Model a

foreign import ccall "wrapper"
  mkModel :: Model a -> IO (FunPtr (Model a))

mkJacobian :: Jacobian a -> IO (FunPtr (Jacobian a))
mkJacobian = mkModel

foreign import ccall "c-levmar/lm.h dlevmar_der"
  c_dlevmar_der :: FunPtr (Model CDouble)    -- func
                -> FunPtr (Jacobian CDouble) -- jacf
                -> Ptr CDouble               -- p
                -> Ptr CDouble               -- x
                -> CInt                      -- m
                -> CInt                      -- n
                -> CInt                      -- itmax
                -> Ptr CDouble               -- opts
                -> Ptr CDouble               -- info
                -> Ptr CDouble               -- work
                -> Ptr CDouble               -- covar
                -> Ptr ()                    -- adata
                -> IO CInt

foreign import ccall "c-levmar/lm.h dlevmar_dif"
  c_dlevmar_dif :: FunPtr (Model CDouble) -- func
                -> Ptr CDouble            -- p
                -> Ptr CDouble            -- x
                -> CInt                   -- m
                -> CInt                   -- n
                -> CInt                   -- itmax
                -> Ptr CDouble            -- opts
                -> Ptr CDouble            -- info
                -> Ptr CDouble            -- work
                -> Ptr CDouble            -- covar
                -> Ptr ()                 -- adata
                -> IO CInt

foreign import ccall "c-levmar/lm.h dlevmar_bc_der"
  c_dlevmar_bc_der :: FunPtr (Model CDouble)    -- func
                   -> FunPtr (Jacobian CDouble) -- jacf
                   -> Ptr CDouble               -- p
                   -> Ptr CDouble               -- x
                   -> CInt                      -- m
                   -> CInt                      -- n
                   -> Ptr CDouble               -- lb
                   -> Ptr CDouble               -- ub
                   -> CInt                      -- itmax
                   -> Ptr CDouble               -- opts
                   -> Ptr CDouble               -- info
                   -> Ptr CDouble               -- work
                   -> Ptr CDouble               -- covar
                   -> Ptr ()                    -- adata
                   -> IO CInt

foreign import ccall "c-levmar/lm.h dlevmar_bc_dif"
  c_dlevmar_bc_dif :: FunPtr (Model CDouble) -- func
                   -> Ptr CDouble            -- p
                   -> Ptr CDouble            -- x
                   -> CInt                   -- m
                   -> CInt                   -- n
                   -> Ptr CDouble            -- lb
                   -> Ptr CDouble            -- ub
                   -> CInt                   -- itmax
                   -> Ptr CDouble            -- opts
                   -> Ptr CDouble            -- info
                   -> Ptr CDouble            -- work
                   -> Ptr CDouble            -- covar
                   -> Ptr ()                 -- adata
                   -> IO CInt

foreign import ccall "c-levmar/lm.h dlevmar_lec_der"
  c_dlevmar_lec_der :: FunPtr (Model CDouble)    -- func
                    -> FunPtr (Jacobian CDouble) -- jacf
                    -> Ptr CDouble               -- p
                    -> Ptr CDouble               -- x
                    -> CInt                      -- m
                    -> CInt                      -- n
                    -> Ptr CDouble               -- A
                    -> Ptr CDouble               -- B
                    -> CInt                      -- k
                    -> CInt                      -- itmax
                    -> Ptr CDouble               -- opts
                    -> Ptr CDouble               -- info
                    -> Ptr CDouble               -- work
                    -> Ptr CDouble               -- covar
                    -> Ptr ()                    -- adata
                    -> IO CInt

foreign import ccall "c-levmar/lm.h dlevmar_lec_dif"
  c_dlevmar_lec_dif :: FunPtr (Model CDouble)    -- func
                    -> Ptr CDouble               -- p
                    -> Ptr CDouble               -- x
                    -> CInt                      -- m
                    -> CInt                      -- n
                    -> Ptr CDouble               -- A
                    -> Ptr CDouble               -- B
                    -> CInt                      -- k
                    -> CInt                      -- itmax
                    -> Ptr CDouble               -- opts
                    -> Ptr CDouble               -- info
                    -> Ptr CDouble               -- work
                    -> Ptr CDouble               -- covar
                    -> Ptr ()                    -- adata
                    -> IO CInt

foreign import ccall "c-levmar/lm.h dlevmar_blec_der"
  c_dlevmar_blec_der :: FunPtr (Model CDouble)    -- func
                     -> FunPtr (Jacobian CDouble) -- jacf
                     -> Ptr CDouble               -- p
                     -> Ptr CDouble               -- x
                     -> CInt                      -- m
                     -> CInt                      -- n
                     -> Ptr CDouble               -- lb
                     -> Ptr CDouble               -- ub
                     -> Ptr CDouble               -- A
                     -> Ptr CDouble               -- B
                     -> CInt                      -- k
                     -> Ptr CDouble               -- wghts
                     -> CInt                      -- itmax
                     -> Ptr CDouble               -- opts
                     -> Ptr CDouble               -- info
                     -> Ptr CDouble               -- work
                     -> Ptr CDouble               -- covar
                     -> Ptr ()                    -- adata
                     -> IO CInt

foreign import ccall "c-levmar/lm.h dlevmar_blec_dif"
  c_dlevmar_blec_dif :: FunPtr (Model CDouble)    -- func
                     -> Ptr CDouble               -- p
                     -> Ptr CDouble               -- x
                     -> CInt                      -- m
                     -> CInt                      -- n
                     -> Ptr CDouble               -- lb
                     -> Ptr CDouble               -- ub
                     -> Ptr CDouble               -- A
                     -> Ptr CDouble               -- B
                     -> CInt                      -- k
                     -> Ptr CDouble               -- wghts
                     -> CInt                      -- itmax
                     -> Ptr CDouble               -- opts
                     -> Ptr CDouble               -- info
                     -> Ptr CDouble               -- work
                     -> Ptr CDouble               -- covar
                     -> Ptr ()                    -- adata
                     -> IO CInt

foreign import ccall "c-levmar/lm.h slevmar_der"
  c_slevmar_der :: FunPtr (Model CFloat)    -- func
                -> FunPtr (Jacobian CFloat) -- jacf
                -> Ptr CFloat               -- p
                -> Ptr CFloat               -- x
                -> CInt                     -- m
                -> CInt                     -- n
                -> CInt                     -- itmax
                -> Ptr CFloat               -- opts
                -> Ptr CFloat               -- info
                -> Ptr CFloat               -- work
                -> Ptr CFloat               -- covar
                -> Ptr ()                   -- adata
                -> IO CInt

foreign import ccall "c-levmar/lm.h slevmar_dif"
  c_slevmar_dif :: FunPtr (Model CFloat) -- func
                -> Ptr CFloat            -- p
                -> Ptr CFloat            -- x
                -> CInt                  -- m
                -> CInt                  -- n
                -> CInt                  -- itmax
                -> Ptr CFloat            -- opts
                -> Ptr CFloat            -- info
                -> Ptr CFloat            -- work
                -> Ptr CFloat            -- covar
                -> Ptr ()                -- adata
                -> IO CInt

foreign import ccall "c-levmar/lm.h slevmar_bc_der"
  c_slevmar_bc_der :: FunPtr (Model CFloat)    -- func
                   -> FunPtr (Jacobian CFloat) -- jacf
                   -> Ptr CFloat               -- p
                   -> Ptr CFloat               -- x
                   -> CInt                     -- m
                   -> CInt                     -- n
                   -> Ptr CFloat               -- lb
                   -> Ptr CFloat               -- ub
                   -> CInt                     -- itmax
                   -> Ptr CFloat               -- opts
                   -> Ptr CFloat               -- info
                   -> Ptr CFloat               -- work
                   -> Ptr CFloat               -- covar
                   -> Ptr ()                   -- adata
                   -> IO CInt

foreign import ccall "c-levmar/lm.h slevmar_bc_dif"
  c_slevmar_bc_dif :: FunPtr (Model CFloat) -- func
                   -> Ptr CFloat            -- p
                   -> Ptr CFloat            -- x
                   -> CInt                  -- m
                   -> CInt                  -- n
                   -> Ptr CFloat            -- lb
                   -> Ptr CFloat            -- ub
                   -> CInt                  -- itmax
                   -> Ptr CFloat            -- opts
                   -> Ptr CFloat            -- info
                   -> Ptr CFloat            -- work
                   -> Ptr CFloat            -- covar
                   -> Ptr ()                -- adata
                   -> IO CInt

foreign import ccall "c-levmar/lm.h slevmar_lec_der"
  c_slevmar_lec_der :: FunPtr (Model CFloat)    -- func
                    -> FunPtr (Jacobian CFloat) -- jacf
                    -> Ptr CFloat               -- p
                    -> Ptr CFloat               -- x
                    -> CInt                     -- m
                    -> CInt                     -- n
                    -> Ptr CFloat               -- A
                    -> Ptr CFloat               -- B
                    -> CInt                     -- k
                    -> CInt                     -- itmax
                    -> Ptr CFloat               -- opts
                    -> Ptr CFloat               -- info
                    -> Ptr CFloat               -- work
                    -> Ptr CFloat               -- covar
                    -> Ptr ()                   -- adata
                    -> IO CInt

foreign import ccall "c-levmar/lm.h slevmar_lec_dif"
  c_slevmar_lec_dif :: FunPtr (Model CFloat) -- func
                    -> Ptr CFloat            -- p
                    -> Ptr CFloat            -- x
                    -> CInt                  -- m
                    -> CInt                  -- n
                    -> Ptr CFloat            -- A
                    -> Ptr CFloat            -- B
                    -> CInt                  -- k
                    -> CInt                  -- itmax
                    -> Ptr CFloat            -- opts
                    -> Ptr CFloat            -- info
                    -> Ptr CFloat            -- work
                    -> Ptr CFloat            -- covar
                    -> Ptr ()                -- adata
                    -> IO CInt

foreign import ccall "c-levmar/lm.h slevmar_blec_der"
  c_slevmar_blec_der :: FunPtr (Model CFloat)    -- func
                     -> FunPtr (Jacobian CFloat) -- jacf
                     -> Ptr CFloat               -- p
                     -> Ptr CFloat               -- x
                     -> CInt                     -- m
                     -> CInt                     -- n
                     -> Ptr CFloat               -- lb
                     -> Ptr CFloat               -- ub
                     -> Ptr CFloat               -- A
                     -> Ptr CFloat               -- B
                     -> CInt                     -- k
                     -> Ptr CFloat               -- wghts
                     -> CInt                     -- itmax
                     -> Ptr CFloat               -- opts
                     -> Ptr CFloat               -- info
                     -> Ptr CFloat               -- work
                     -> Ptr CFloat               -- covar
                     -> Ptr ()                   -- adata
                     -> IO CInt

foreign import ccall "c-levmar/lm.h slevmar_blec_dif"
  c_slevmar_blec_dif :: FunPtr (Model CFloat) -- func
                     -> Ptr CFloat            -- p
                     -> Ptr CFloat            -- x
                     -> CInt                  -- m
                     -> CInt                  -- n
                     -> Ptr CFloat            -- lb
                     -> Ptr CFloat            -- ub
                     -> Ptr CFloat            -- A
                     -> Ptr CFloat            -- B
                     -> CInt                  -- k
                     -> Ptr CFloat            -- wghts
                     -> CInt                  -- itmax
                     -> Ptr CFloat            -- opts
                     -> Ptr CFloat            -- info
                     -> Ptr CFloat            -- work
                     -> Ptr CFloat            -- covar
                     -> Ptr ()                -- adata
                     -> IO CInt
