{-# LANGUAGE ForeignFunctionInterface #-}

module LevMar ( Model
              , Jacobian

              , withModel
              , withJacobian

              , C_LevMarDer
              , C_LevMarDif
              , C_LevMarBCDer
              , C_LevMarBCDif
              , C_LevMarLecDer
              , C_LevMarLecDif
              , C_LevMarBLecDer
              , C_LevMarBLecDif

              , c_dlevmar_der
              , c_slevmar_der
              , c_dlevmar_dif
              , c_slevmar_dif
              , c_dlevmar_bc_der
              , c_slevmar_bc_der
              , c_dlevmar_bc_dif
              , c_slevmar_bc_dif
              , c_dlevmar_lec_der
              , c_slevmar_lec_der
              , c_dlevmar_lec_dif
              , c_slevmar_lec_dif
              , c_dlevmar_blec_der
              , c_slevmar_blec_der
              , c_dlevmar_blec_dif
              , c_slevmar_blec_dif

              , PureModel
              , Options(..)
              , StopReason(..)
              , Info(..)
              , CoVarMatrix
              , LevMarDif
              , defaultOpts
              , dlevmar_dif
              , slevmar_dif
              ) where

import Foreign.Marshal.Array
import Foreign.C.Types
import Foreign.Ptr
import Foreign.Storable
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

type C_LevMarDer cr =  FunPtr (Model cr)    -- func
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

type C_LevMarDif cr =  FunPtr (Model cr) -- func
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

type C_LevMarBCDer cr =  FunPtr (Model cr)    -- func
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

type C_LevMarBCDif cr =  FunPtr (Model cr) -- func
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

type C_LevMarLecDer cr =  FunPtr (Model cr)    -- func
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

type C_LevMarLecDif cr =  FunPtr (Model cr) -- func
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

type C_LevMarBLecDer cr =  FunPtr (Model cr)    -- func
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

type C_LevMarBLecDif cr =  FunPtr (Model cr) -- func
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

foreign import ccall "dlevmar_der"      c_dlevmar_der      :: C_LevMarDer     CDouble
foreign import ccall "slevmar_der"      c_slevmar_der      :: C_LevMarDer     CFloat
foreign import ccall "dlevmar_dif"      c_dlevmar_dif      :: C_LevMarDif     CDouble
foreign import ccall "slevmar_dif"      c_slevmar_dif      :: C_LevMarDif     CFloat
foreign import ccall "dlevmar_bc_der"   c_dlevmar_bc_der   :: C_LevMarBCDer   CDouble
foreign import ccall "slevmar_bc_der"   c_slevmar_bc_der   :: C_LevMarBCDer   CFloat
foreign import ccall "dlevmar_bc_dif"   c_dlevmar_bc_dif   :: C_LevMarBCDif   CDouble
foreign import ccall "slevmar_bc_dif"   c_slevmar_bc_dif   :: C_LevMarBCDif   CFloat
foreign import ccall "dlevmar_lec_der"  c_dlevmar_lec_der  :: C_LevMarLecDer  CDouble
foreign import ccall "slevmar_lec_der"  c_slevmar_lec_der  :: C_LevMarLecDer  CFloat
foreign import ccall "dlevmar_lec_dif"  c_dlevmar_lec_dif  :: C_LevMarLecDif  CDouble
foreign import ccall "slevmar_lec_dif"  c_slevmar_lec_dif  :: C_LevMarLecDif  CFloat
foreign import ccall "dlevmar_blec_der" c_dlevmar_blec_der :: C_LevMarBLecDer CDouble
foreign import ccall "slevmar_blec_der" c_slevmar_blec_der :: C_LevMarBLecDer CFloat
foreign import ccall "dlevmar_blec_dif" c_dlevmar_blec_dif :: C_LevMarBLecDif CDouble
foreign import ccall "slevmar_blec_dif" c_slevmar_blec_dif :: C_LevMarBLecDif CFloat

-------------------------------------------------------------------------------

type PureModel r = [r] -> Int -> r

data Options r = Opts { opt_mu       :: r
                      , opt_epsilon1 :: r
                      , opt_epsilon2 :: r
                      , opt_epsilon3 :: r
                      , opt_delta    :: r
                      } deriving Show

data StopReason = SmallGradient
                | SmallDp
                | MaxIterations
                | SingularMatrix
                | SmallestError
                | SmallE_2
                | InvalidValues
                  deriving (Show, Enum)

data Info r = Info { inf_values          :: [r]
                   , inf_numIter         :: Integer
                   , inf_stopReason      :: StopReason
                   , inf_numFuncEvals    :: Integer
                   , inf_numJacobEvals   :: Integer
                   , inf_numLinSysSolved :: Integer
                   } deriving Show

type CoVarMatrix r = [r]

type LevMarDif r =  PureModel r
                 -> [r]       -- initial parameters
                 -> [r]       -- samples
                 -> Integer   -- itmax
                 -> Options r -- opts
                 -> IO ([r], Info r, CoVarMatrix r)

defaultOpts :: Fractional r => Options r
defaultOpts = Opts { opt_mu       = 1e-3
                   , opt_epsilon1 = 1e-17
                   , opt_epsilon2 = 1e-17
                   , opt_epsilon3 = 1e-17
                   , opt_delta    = 1e-6
                   }

optsToList :: Options r -> [r]
optsToList (Opts mu eps1 eps2 eps3 delta) = [mu, eps1, eps2, eps3, delta]

listToInfo :: RealFrac r => [r] -> Info r
listToInfo [a,b,c,d,e,f,g,h,i,j] =
    Info { inf_values          = [a,b,c,d,e]
         , inf_numIter         = floor f
         , inf_stopReason      = toEnum $ floor g - 1
         , inf_numFuncEvals    = floor h
         , inf_numJacobEvals   = floor i
         , inf_numLinSysSolved = floor j
         }
listToInfo _ = error "liftToInfo: wrong list length"

convertModel :: (Real h, Fractional h, Storable c, Real c, Fractional c)
             => PureModel h -> Model c
convertModel f = \p hx m n _ -> do
                   ps <- peekArray (fromIntegral m) p
                   let ps' = map realToFrac ps
                       ns  = [0 .. fromIntegral n - 1]
                       hx' = map (realToFrac . f ps') ns
                   pokeArray hx hx'


gen_levmar_dif :: (Storable cr, Real cr, Fractional cr, RealFrac r)
               => C_LevMarDif cr -> LevMarDif r
gen_levmar_dif minimise f ps xs itMax opts = let lenPs = length ps in
    withArray (map realToFrac ps) $ \psPtr ->
      withArray (map realToFrac xs) $ \xsPtr ->
        withArray (map realToFrac $ optsToList opts) $ \optsPtr ->
          allocaArray 10 $ \infoPtr ->
            allocaArray (lenPs * lenPs) $ \coVarPtr ->
              withModel (convertModel f) $ \modelPtr -> do

                _ <- minimise modelPtr
                              psPtr
                              xsPtr
                              (fromIntegral $ lenPs)
                              (fromIntegral $ length xs)
                              (fromIntegral itMax)
                              optsPtr
                              infoPtr
                              nullPtr -- work
                              coVarPtr
                              nullPtr -- adata

                result <- peekArray lenPs psPtr
                info   <- peekArray 10 infoPtr
                coVar  <- peekArray (lenPs * lenPs) coVarPtr

                return $ ( map realToFrac result
                         , listToInfo $ map realToFrac info
                         , map realToFrac coVar
                         )


dlevmar_dif :: LevMarDif Double
dlevmar_dif = gen_levmar_dif c_dlevmar_dif

slevmar_dif :: LevMarDif Float
slevmar_dif = gen_levmar_dif c_slevmar_dif
