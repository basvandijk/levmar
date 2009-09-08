-- This module is a Haskell translation of lmdemo.c from the C levmar library.

module Demo where

import LevMar ( levmar

              , Model
              , Jacobian

              , Options(..), defaultOpts

              , LinearConstraints, noLinearConstraints

              , LevMarError

              , Info, CovarMatrix

              , S, Z
              , SizedList(..)
              )

import qualified LevMar.Fitting as Fitting

type Result n = Either LevMarError
                       ( SizedList n Double
                       , Info Double
                       , CovarMatrix n Double
                       )

--------------------------------------------------------------------------------
-- Handy type synonyms for type-level naturals:

type N0 = Z
type N1 = S N0
type N2 = S N1
type N3 = S N2
type N4 = S N3
type N5 = S N4

--------------------------------------------------------------------------------
-- Default options:

opts :: Options Double
opts = defaultOpts { optEpsilon1 = 1e-15
                   , optEpsilon2 = 1e-15
                   , optEpsilon3 = 1e-20
                   }

--------------------------------------------------------------------------------
-- Rosenbrock function,
-- global minimum at (1, 1)

ros :: Model N2 Double
ros p0 p1 = replicate ros_n ((1.0 - p0)**2 + ros_d*m**2)
    where
      m = p1 - p0**2

ros_jac :: Jacobian N2 Double
ros_jac p0 p1 = replicate ros_n (p0d ::: p1d ::: Nil)
    where
      p0d = -2 + 2*p0 - 4*ros_d*m*p0
      p1d = 2*ros_d*m
      m   = p1 - p0**2

ros_d :: Double
ros_d = 105.0

ros_n :: Int
ros_n = 2

ros_params :: SizedList N2 Double
ros_params = -1.2 ::: 1.0 ::: Nil

ros_samples :: [Double]
ros_samples = replicate ros_n 0.0

run_ros :: Result N2
run_ros = levmar ros
                 (Just ros_jac)
                 ros_params
                 ros_samples
                 1000
                 opts
                 Nothing
                 Nothing
                 noLinearConstraints
                 Nothing

{- run_ros =>
Right ( 0.9439745509911915 ::: 0.8908311893490423 ::: Nil
      , Info { infValues = [ 1266.8564479999995
                           , 1.9791779654508303e-5
                           , 2.8077694910189473e-4
                           , 6.626244087850897e-8
                           , 108.82463493455907
                           ]
             , infNumIter         = 1000
             , infStopReason      = MaxIterations
             , infNumFuncEvals    = 1283
             , infNumJacobEvals   = 1000
             , infNumLinSysSolved = 1282
             }
      ,     (1.033011505608052e-3    ::: 1.6140246295044776e-3   ::: Nil)
        ::: (1.2624987557696709e-247 ::: 1.2625052477848744e-247 ::: Nil)
        ::: Nil
      )
-}

--------------------------------------------------------------------------------
-- Modified Rosenbrock problem,
-- global minimum at (1, 1)

modros :: Model N2 Double
modros p0 p1 = [ 10*(p1 - p0**2)
               , 1.0 - p0
               , modros_lam
               ]

modros_jac :: Jacobian N2 Double
modros_jac p0 _ = [ -20*p0 ::: 10.0 ::: Nil
                  , -1.0   ::: 0.0  ::: Nil
                  , 0.0    ::: 0.0  ::: Nil
                  ]

modros_lam :: Double
modros_lam = 1e02

modros_n :: Int
modros_n = 3

modros_params :: SizedList N2 Double
modros_params = -1.2 ::: 1.0 ::: Nil

modros_samples :: [Double]
modros_samples = replicate modros_n 0.0

run_modros :: Result N2
run_modros = levmar modros
                    (Just modros_jac)
                    modros_params
                    modros_samples
                    1000
                    opts
                    Nothing
                    Nothing
                    noLinearConstraints
                    Nothing

{- run_modros =>
Right ( 0.999999200793784 ::: 0.99999839572234 ::: Nil
      , Info { infValues = [ 10024.2
                           , 10000.0
                           , 5.865866792298391e-7
                           , 3.7449494489959514e-34
                           , 8.964512911192581e7
                           ]
             , infNumIter         = 14
             , infStopReason      = SmallDp
             , infNumFuncEvals    = 25
             , infNumJacobEvals   = 14
             , infNumLinSysSolved = 25
             }
      ,     ( 9999.999999999893         ::: 19999.984015875463        ::: Nil)
        ::: (   -2.7989760846248645e-34 :::    -2.799005471954965e-34 ::: Nil)
        ::: Nil
      )
-}

--------------------------------------------------------------------------------
-- Powell's function,
-- minimum at (0, 0)

powell :: Model N2 Double
powell p0 p1 = [ p0
               , 10.0*p0 / m + 2*p1**2
               ]
    where
      m = p0 + 0.1

powell_jac :: Jacobian N2 Double
powell_jac p0 p1 = [ 1.0        ::: 0.0    ::: Nil
                   , 1.0 / m**2 ::: 4.0*p1 ::: Nil
                   ]
    where
      m = p0 + 0.1

powell_n :: Int
powell_n = 2

powell_params :: SizedList N2 Double
powell_params = -1.2 ::: 1.0 ::: Nil

powell_samples :: [Double]
powell_samples = replicate powell_n 0.0

run_powell :: Result N2
run_powell = levmar powell
                    (Just powell_jac)
                    powell_params
                    powell_samples
                    1000
                    opts
                    Nothing
                    Nothing
                    noLinearConstraints
                    Nothing

{- run_powell =>
Right ( -8.768525913529536e-11 ::: -6.631546508520969e-5 ::: Nil
      , Info { infValues = [ 168.08462809917359
                           , 8.415325091358066e-21
                           , 2.8996497517975197e-9
                           , 1.303952261193396e-11
                           , 6.108016441269007e-15
                           ]
             , infNumIter         = 182
             , infStopReason      = SmallE_2
             , infNumFuncEvals    = 185
             , infNumJacobEvals   = 182
             , infNumLinSysSolved = 184
             }
      ,     (Infinity              ::: Infinity              ::: Nil)
        ::: (6.94996061456946e-310 ::: 6.94996061456946e-310 ::: Nil)
        ::: Nil
      )
-}

--------------------------------------------------------------------------------
-- Wood's function,
-- minimum at (1, 1, 1, 1)

wood :: Model N4 Double
wood p0 p1 p2 p3 = [ 10.0*(p1 - p0**2)
                   , 1.0 - p0
                   , sqrt 90.0*(p3 - p2**2)
                   , 1.0 - p2
                   , sqrt 10.0*(p1 + p3 - 2.0)
                   , (p1 - p3) / sqrt 10.0
                   ]

wood_n :: Int
wood_n = 6

wood_params :: SizedList N4 Double
wood_params =  -3.0 ::: -1.0 ::: -3.0 ::: -1.0 ::: Nil

wood_samples :: [Double]
wood_samples = replicate wood_n 0.0

run_wood :: Result N4
run_wood = levmar wood
                  Nothing
                  wood_params
                  wood_samples
                  1000
                  opts
                  Nothing
                  Nothing
                  noLinearConstraints
                  Nothing

{- run_wood =>
Right ( 0.9999999999991299 ::: 0.999999999998257 ::: 1.0000000000008702 ::: 1.0000000000017433 ::: Nil
      , Info { infValues = [ 19192.0
                           , 2.731305324788736e-24
                           , 5.26965526926442e-10
                           , 1.4466573640990308e-18
                           , 6.850442855162863e-7
                           ]
             , infNumIter         = 113
             , infStopReason      = SmallE_2
             , infNumFuncEvals    = 158
             , infNumJacobEvals   = 11
             , infNumLinSysSolved = 113
             }
      ,     (3.8980269258203757e-25  ::: 7.747272602655087e-25   ::: -3.695213667010158e-25  ::: -7.410945761969915e-25  ::: Nil)
        ::: (-3.997387753537078e250  ::: -1.7465702694268607e221 ::: -9.137793733109449e269  ::: -3.9973836573137235e250 ::: Nil)
        ::: (7.747272602655087e-25   ::: -3.695213667010158e-25  ::: -7.410945761969915e-25  ::: 7.747272602655085e-25   ::: Nil)
        ::: (-1.7465702694268607e221 ::: -9.137793733109449e269  ::: -3.9973836573137235e250 ::: -1.74909068114927e10    ::: Nil)
        ::: Nil
      )
-}

--------------------------------------------------------------------------------
-- Meyer's (reformulated) data fitting problem,
-- minimum at (2.48, 6.18, 3.45)

meyer :: Fitting.Model N3 Double Double
meyer p0 p1 p2 x = p0*exp (10.0*p1 / (ui + p2) - 13.0)
    where
      ui = 0.45 + 0.05*x

meyer_jac :: Fitting.Jacobian N3 Double Double
meyer_jac p0 p1 p2 x =     tmp
                       ::: 10.0*p0*tmp / (ui + p2)
                       ::: -10.0*p0*p1*tmp / ((ui + p2)*(ui + p2))
                       ::: Nil
    where
      tmp = exp (10.0*p1 / (ui + p2) - 13.0)
      ui = 0.45 + 0.05*x

meyer_n :: Int
meyer_n = 16

meyer_params :: SizedList N3 Double
meyer_params = 8.85 ::: 4.0 ::: 2.5 ::: Nil

meyer_samples :: [(Double, Double)]
meyer_samples =  zip [0..] [ 34.780
                           , 28.610
                           , 23.650
                           , 19.630
                           , 16.370
                           , 13.720
                           , 11.540
                           , 9.744
                           , 8.261
                           , 7.030
                           , 6.005
                           , 5.147
                           , 4.427
                           , 3.820
                           , 3.307
                           , 2.872
                           ]

run_meyer_jac :: Result N3
run_meyer_jac = Fitting.levmar meyer
                               (Just meyer_jac)
                               meyer_params
                               meyer_samples
                               1000
                               opts
                               Nothing
                               Nothing
                               noLinearConstraints
                               Nothing

{- run_meyer_jac =>
Right ( 2.4817783023532543 ::: 6.18134634518888 ::: 3.5022363458721086 ::: Nil
      , Info { infValues = [ 1308.2529235783877
                           , 8.794585516991546e-5
                           , 3.931110528787064e-11
                           , 1.2607019026963552e-30
                           , 8.313370657962425e-2
                           ]
             , infNumIter         = 83
             , infStopReason      = SmallDp
             , infNumFuncEvals    = 93
             , infNumJacobEvals   = 83
             , infNumLinSysSolved = 93
             }
      ,     (4.817100848339499e-3   ::: -1.6173017206961584e-3 ::: -5.441368285215037e-4  ::: Nil)
        ::: (-8.783861814211665e162 ::: 2.375201093181698e-271 ::: -6.265745407357833e162 ::: Nil)
        ::: (8.997040982923626e305  ::: -1.9781990768516175e99 ::: 8.997040782279326e305  ::: Nil)
        ::: Nil
      )
-}

run_meyer :: Result N3
run_meyer = Fitting.levmar meyer
                           Nothing
                           meyer_params
                           meyer_samples
                           1000
                           opts
                           Nothing
                           Nothing
                           noLinearConstraints
                           Nothing

{- run_meyer =>
Right ( 2.481778352780817 ::: 6.181346326650312 ::: 3.5022363391220463 ::: Nil
      , Info { infValues = [ 1308.2529235783877
                           , 8.794585517057738e-5
                           , 1.1000124576130386e-7
                           , 2.897020822609385e-33
                           , 38949.63759231421
                           ]
             , infNumIter         = 209
             , infStopReason      = SmallDp
             , infNumFuncEvals    = 273
             , infNumJacobEvals   = 21
             , infNumLinSysSolved = 210
             }
      ,     ( 4.8351423319113965e-3  ::: -1.6244524410218796e-3  ::: -5.481136912988729e-4   ::: Nil)
        ::: (-4.6422588659052177e207 :::  1.4304490480687128e263 ::: -9.579860655012518e207  ::: Nil)
        ::: (-7.854672410559325e-167 ::: -9.682389150573423e259  ::: -7.854672838122991e-167 ::: Nil)
        ::: Nil
      )
-}

--------------------------------------------------------------------------------
-- helical valley function,
-- minimum at (1.0, 0.0, 0.0)

helval :: Model  N3 Double
helval p0 p1 p2 = [ 10.0*(p2 - 10.0*theta)
                  , 10.0*sqrt tmp - 1.0
                  , p2
                  ]
    where
      m = atan (p1 / p0) / (2.0*pi)

      tmp = p0**2 + p1**2

      theta | p0 < 0.0  = m + 0.5
            | 0.0 < p0  = m
            | p1 >= 0   = 0.25
            | otherwise = -0.25

heval_jac :: Jacobian N3 Double
heval_jac p0 p1 _ = [ 50.0*p1 / (pi*tmp) ::: -50.0*p0 / (pi*tmp) ::: 10.0 ::: Nil
                    , 10.0*p0 / sqrt tmp :::  10.0*p1 / sqrt tmp ::: 0.0  ::: Nil
                    , 0.0                ::: 0.0                 ::: 1.0  ::: Nil
                    ]
    where
      tmp = p0**2 + p1**2

helval_n :: Int
helval_n = 3

helval_params :: SizedList N3 Double
helval_params = -1.0 ::: 0.0 ::: 0.0 ::: Nil

helval_samples :: [Double]
helval_samples = replicate helval_n 0.0

run_helval :: Result N3
run_helval = levmar helval
                    (Just heval_jac)
                    helval_params
                    helval_samples
                    1000
                    opts
                    Nothing
                    Nothing
                    noLinearConstraints
                    Nothing

{- run_helval =>
Right ( 0.1 ::: 2.880510934753633e-13 ::: 4.584295271122304e-12 ::: Nil
      , Info { infValues = [ 2581.0
                           , 2.1015766382483756e-23
                           , 3.3595896608561804e-9
                           , 1.1366676576712879e-17
                           , 1.7858243426552675e-8
                           ]
             , infNumIter         = 9
             , infStopReason      = SmallE_2
             , infNumFuncEvals    = 10
             , infNumJacobEvals   = 9
             , infNumLinSysSolved = 9
             }
      , (   Infinity             ::: Infinity            ::: -Infinity         ::: Nil)
        ::: (2.714828365066e-312 ::: 2.714828365066e-312 ::: 5.4309829775e-312 ::: Nil)
        ::: (1.61816e-319        ::: 1.61816e-319        ::: 3.2371e-319       ::: Nil)
        ::: Nil
      )
-}

--------------------------------------------------------------------------------
-- Boggs - Tolle problem 3 (linearly constrained),
-- minimum at (-0.76744, 0.25581, 0.62791, -0.11628, 0.25581)
--
-- constr1: p0 + 3*p1      = 0
-- constr2: p2 + p3 - 2*p4 = 0
-- constr3: p1 - p4          = 0

bt3 :: Model N5 Double
bt3 p0 p1 p2 p3 p4 = replicate bt3_n (t1**2 + t2**2 + t3**2 + t4**2)
    where
      t1 = p0 - p1
      t2 = p1 + p2 - 2.0
      t3 = p3 - 1.0
      t4 = p4 - 1.0

bt3_jac :: Jacobian N5 Double
bt3_jac p0 p1 p2 p3 p4 = replicate bt3_n (   2.0*t1
                                         ::: 2.0*(t2 - t1)
                                         ::: 2.0*t2
                                         ::: 2.0*t3
                                         ::: 2.0*t4
                                         ::: Nil
                                         )
    where
      t1 = p0 - p1
      t2 = p1 + p2 - 2.0
      t3 = p3 - 1.0
      t4 = p4 - 1.0

bt3_n :: Int
bt3_n = 5

bt3_params :: SizedList N5 Double
bt3_params = 2.0 ::: 2.0 ::: 2.0 :::2.0 ::: 2.0 ::: Nil

bt3_samples :: [Double]
bt3_samples = replicate bt3_n 0.0

bt3_linear_constraints :: LinearConstraints N3 N5 Double
bt3_linear_constraints = (     (1.0 ::: 3.0 ::: 0.0 ::: 0.0 :::  0.0 ::: Nil)
                           ::: (0.0 ::: 0.0 ::: 1.0 ::: 1.0 ::: -2.0 ::: Nil)
                           ::: (0.0 ::: 1.0 ::: 0.0 ::: 0.0 ::: -1.0 ::: Nil)
                           ::: Nil
                         , 0.0 ::: 0.0 ::: 0.0 ::: Nil
                         )

run_bt3 :: Result N5
run_bt3 = levmar bt3
                (Just bt3_jac)
                bt3_params
                bt3_samples
                1000
                opts
                Nothing
                Nothing
                (Just bt3_linear_constraints)
                Nothing

{- run_bt3 =>
Right ( -0.7674418626339431 ::: 0.2558139542113145 ::: 0.6279069703819862 ::: -0.11627906195935662 ::: 0.2558139542113146 ::: Nil
      , Info { infValues = [ 109.25737894331408
                           , 83.76419686316918
                           , 3.840097938822601e-7
                           , 1.2133299385969255e-33
                           , 6.553371761836997e24
                           ]
             , infNumIter         = 19
             , infStopReason      = SmallDp
             , infNumFuncEvals    = 31
             , infNumJacobEvals   = 19
             , infNumLinSysSolved = 31
             }
      ,     (6.516494898922502e-2   ::: 5.924086429976833e-3    ::: 7.108903541920185e-2     ::: 7.108903456282224e-2    ::: 4.7392689994008044e-2  ::: Nil)
        ::: (-4.090874204742538e129 ::: -1.6218479364970034e174 ::: -1.3944200639589396e174  ::: -2.355682872237576e129  ::: -4.090874308031206e129 ::: Nil)
        ::: (-3.47923313434499e52   ::: 1.1432257313549235e-70  ::: -3.7295837473563526e296  ::: 1.1675986418530346e127  ::: -2.991871737350886e150 ::: Nil)
        ::: (-7.007070698633969e-67 ::: -7.00704988119892e-67   ::: -1.0077995925212267e-114 ::: 2.0978419878214583e271  ::: 1.8213966943589208e8   ::: Nil)
        ::: (4.767292522510342e-200 ::: -1.5320295036637725e54  ::: -1.0524953945184693e114  ::: 1.6526571269395338e261  ::: -7.90012780402775e-124 ::: Nil)
        ::: Nil
      )
-}

--------------------------------------------------------------------------------
-- Hock - Schittkowski problem 28 (linearly constrained),
-- minimum at (0.5, -0.5, 0.5)
--
-- constr1: p0 + 2*p1 + 3*p2 = 1

hs28 :: Model N3 Double
hs28 p0 p1 p2 = replicate hs28_n (t1**2 + t2**2)
    where
      t1 = p0 + p1
      t2 = p1 + p2

hs28_jac :: Jacobian N3 Double
hs28_jac p0 p1 p2 = replicate hs28_n (     2.0*t1
                                       ::: 2.0*(t1 + t2)
                                       ::: 2.0*t2
                                       ::: Nil
                                     )
    where
      t1 = p0 + p1
      t2 = p1 + p2

hs28_n :: Int
hs28_n = 3

hs28_params :: SizedList N3 Double
hs28_params = -4.0 ::: 1.0 ::: 1.0 ::: Nil

hs28_samples :: [Double]
hs28_samples = replicate hs28_n 0.0

hs28_linear_constraints :: LinearConstraints N1 N3 Double
hs28_linear_constraints = ( ((1.0 ::: 2.0 ::: 3.0 ::: Nil) ::: Nil)
                          , 1.0 ::: Nil
                          )

run_hs28 :: Result N3
run_hs28 = levmar hs28
                  (Just hs28_jac)
                  hs28_params
                  hs28_samples
                  1000
                  opts
                  Nothing
                  Nothing
                  (Just hs28_linear_constraints)
                  Nothing

{- run_hs28 =>
Right ( 0.49999041471671685 ::: -0.4999905696535188 ::: 0.4999969081967735 ::: Nil
      , Info { infValues = [ 506.9999999999996
                           , 4.848394008150168e-21
                           , 1.816909882595112e-15
                           , 2.6576610168803227e-11
                           , 0.43002370636291093
                           ]
             , infNumIter         = 41
             , infStopReason      = SmallE_2
             , infNumFuncEvals    = 42
             , infNumJacobEvals   = 41
             , infNumLinSysSolved = 41
             }
      ,     (1.1225364286455017e-14 ::: 8.628338628002167e-14  ::: 7.505802199356666e-14 ::: Nil)
        ::: (8.05790282502707e181   ::: 3.5992518324782235e102 ::: 8.048178482448445e181 ::: Nil)
        ::: (2.8138652636265918e45  ::: 3.847881895971628e-152 ::: 2.8138652635019906e45 ::: Nil)
        ::: Nil
      )
-}

--------------------------------------------------------------------------------
-- Hock - Schittkowski problem 48 (linearly constrained),
-- minimum at (1.0, 1.0, 1.0, 1.0, 1.0)
--
-- constr1: sum [p0, p1, p2, p3, p4] = 5
-- constr2: p2 - 2*(p3 + p4)       = -3

hs48 :: Model N5 Double
hs48 p0 p1 p2 p3 p4 = replicate hs48_n (t1**2 + t2**2 + t3**2)
    where
      t1 = p0 - 1.0
      t2 = p1 - p2
      t3 = p3 - p4

hs48_jac :: Jacobian N5 Double
hs48_jac p0 p1 p2 p3 p4 = replicate hs48_n (     2.0*t1
                                             ::: 2.0*t2
                                             ::: 2.0*t2
                                             ::: 2.0*t3
                                             ::: 2.0*t3
                                             ::: Nil
                                           )
    where
      t1 = p0 - 1.0
      t2 = p1 - p2
      t3 = p3 - p4

hs48_n :: Int
hs48_n = 3

hs48_params :: SizedList N5 Double
hs48_params = 3.0 ::: 5.0 ::: -3.0 ::: 2.0 ::: -2.0 ::: Nil

hs48_samples :: [Double]
hs48_samples = replicate hs48_n 0.0

hs48_linear_constraints :: LinearConstraints N2 N5 Double
hs48_linear_constraints = (     (1.0 ::: 1.0 ::: 1.0 :::  1.0 :::  1.0 ::: Nil)
                            ::: (0.0 ::: 0.0 ::: 1.0 ::: -2.0 ::: -2.0 ::: Nil)
                            ::: Nil
                          , 5.0 ::: -3.0 ::: Nil
                          )

run_hs48 :: Result N5
run_hs48 = levmar hs48
                  (Just hs48_jac)
                  hs48_params
                  hs48_samples
                  1000
                  opts
                  Nothing
                  Nothing
                  (Just hs48_linear_constraints)
                  Nothing

{- run_hs48 =>
Right ( 3.0191758865932465 ::: 4.989156301590405 ::: -3.0055547921224344 ::: 1.9986113019693914 ::: -2.0013886980306084 ::: Nil
      , Info { infValues = [ 21168.0
                           , 21164.208443197764
                           , 1505.7000055524745
                           , 3.3227604671999094e-29
                           , 3.7039401275645075e15
                           ]
             , infNumIter         = 48
             , infStopReason      = SmallDp
             , infNumFuncEvals    = 55
             , infNumJacobEvals   = 48
             , infNumLinSysSolved = 55
             }
      ,     (0.13382566584315234     ::: 0.5298684193108543      ::: 0.5298684193108543       ::: 0.26510947705292354      ::: 0.26510947705292354    ::: Nil)
        ::: (-3.7983076882441762e-84 ::: -3.798307688247901e-84  ::: -1.7252659936150142e-242 ::: -1.7252659936141623e-242 ::: -3.798307734268781e-84 ::: Nil)
        ::: (-7.001288809770295e242  ::: 1.1602709897084776e94   ::: 1.1602582671592755e94    ::: -8.668035929883373e242   ::: 9.979072408262554e99   ::: Nil)
        ::: (-1.1358269589262246e159 ::: -1.117113915045771e82   ::: -1.117113915045771e82    ::: -1.1358269589265243e159  ::: 7.011066374900991e-305 ::: Nil)
        ::: (-1.0985572256787802e-12 ::: -1.0985572254671704e-12 ::: -1.6913178375115177e-166 ::: -8.291774969382475e151   ::: -8.291774971808222e151 ::: Nil)
        ::: Nil
      )
-}

--------------------------------------------------------------------------------
-- Hock - Schittkowski problem 51 (linearly constrained),
-- minimum at (1.0, 1.0, 1.0, 1.0, 1.0)
--
-- constr1: p0 + 3*p1      = 4
-- constr2: p2 + p3 - 2*p4 = 0
-- constr3: p1 - p4          = 0

hs51 :: Model N5 Double
hs51 p0 p1 p2 p3 p4 = replicate hs51_n (t1**2 + t2**2 + t3**2 + t4**2)
    where
      t1 = p0 - p1
      t2 = p1 + p2 - 2.0
      t3 = p3 - 1.0
      t4 = p4 - 1.0

hs51_jac :: Jacobian N5 Double
hs51_jac p0 p1 p2 p3 p4 = replicate hs51_n (     2.0*t1
                                             ::: 2.0*(t2 - t1)
                                             ::: 2.0*t2
                                             ::: 2.0*t3
                                             ::: 2.0*t4
                                             ::: Nil
                                           )
    where
      t1 = p0 - p1
      t2 = p1 + p2 - 2.0
      t3 = p3 - 1.0
      t4 = p4 - 1.0

hs51_n :: Int
hs51_n = 5

hs51_params :: SizedList N5 Double
hs51_params = 2.5 ::: 0.5 ::: 2.0 ::: -1.0 ::: 0.5 ::: Nil

hs51_samples :: [Double]
hs51_samples = replicate hs51_n 0.0

hs51_linear_constraints :: LinearConstraints N3 N5 Double
hs51_linear_constraints = (     (1.0 ::: 3.0 ::: 0.0 ::: 0.0 :::  0.0 ::: Nil)
                            ::: (0.0 ::: 0.0 ::: 1.0 ::: 1.0 ::: -2.0 ::: Nil)
                            ::: (0.0 ::: 1.0 ::: 0.0 ::: 0.0 ::: -1.0 ::: Nil)
                            ::: Nil
                          , 4.0 ::: 0.0 ::: 0.0 ::: Nil
                          )

run_hs51 :: Result N5
run_hs51 = levmar hs51
                  (Just hs51_jac)
                  hs51_params
                  hs51_samples
                  1000
                  opts
                  Nothing
                  Nothing
                  (Just hs51_linear_constraints)
                  Nothing

{- run_hs51 =>
Right ( 1.000000975349995 ::: 0.999999674883335 ::: 1.0000028194247943 ::: 0.999996530341875 ::: 0.9999996748833347 ::: Nil
      , Info { infValues = [ 361.25
                           , 2.0114192392578446e-21
                           , 3.3349798931125986e-15
                           , 7.288731283255759e-12
                           , 0.18336864567974345
                           ]
             , infNumIter         = 21
             , infStopReason      = SmallE_2
             , infNumFuncEvals    = 22
             , infNumJacobEvals   = 21
             , infNumLinSysSolved = 21
             }
      ,     (3.6582684719817334e-14  ::: 3.358327237732954e-14   ::: 7.016595709714686e-14   ::: -9.760297064930916e-14 ::: -9.145671183398134e-15  ::: Nil)
        ::: (-2.88187611099691e135   ::: -2.0472185018081177e296 ::: 3.2909149569972555e-265 ::: 4.352803377565251e-160 ::: -2.88187610287121e135   ::: Nil)
        ::: (-8.424787603029281e-46  ::: 5.650025491469458e-261  ::: -3.276673354691396e-279 ::: 4.132487030008717e-51  ::: -1.928724057589641e89   ::: Nil)
        ::: (1.9907458102175382e-59  ::: 3.313385888533072e-22   ::: 1.5840550679903494e-286 ::: 1.565292349747021e-140 ::: 5.286708966900749e-148  ::: Nil)
        ::: (1.6838187431942556e-283 ::: -1.672151245959948e118  ::: 6.091131102164253e-233  ::: 4.601067906299745e-80  ::: 1.6763426770031812e-156 ::: Nil)
        ::: Nil
      )
-}

--------------------------------------------------------------------------------
-- Hock - Schittkowski problem 01 (box constrained),
-- minimum at (1.0, 1.0)
--
-- constr1: p1 >= -1.5

hs01 :: Model N2 Double
hs01 p0 p1 = [ 10.0*(p1 - p0**2)
             , 1.0 - p0
             ]

hs01_jac :: Jacobian N2 Double
hs01_jac p0 _ = [ -20.0*p0 ::: 10.0 ::: Nil
                , -1.0     ::: 0.0  ::: Nil
                ]

hs01_n :: Int
hs01_n = 2

hs01_params :: SizedList N2 Double
hs01_params = -2.0 ::: 1.0 ::: Nil

hs01_samples :: [Double]
hs01_samples = replicate hs01_n 0.0

hs01_lb, hs01_ub :: SizedList N2 Double
hs01_lb = -_DBL_MAX ::: -1.5     ::: Nil
hs01_ub =  _DBL_MAX ::: _DBL_MAX ::: Nil

_DBL_MAX :: Double
_DBL_MAX = 1e+37 -- TODO: Get this directly from <float.h>.

run_hs01 :: Result N2
run_hs01 = levmar hs01
                  (Just hs01_jac)
                  hs01_params
                  hs01_samples
                  1000
                  opts
                  (Just hs01_lb)
                  (Just hs01_ub)
                  noLinearConstraints
                  Nothing

{-run_hs01 =>
Right ( 0.9999999999551645 ::: 0.9999999999101493 ::: Nil
      , Info { infValues = [ 909.0
                           , 2.013450081327987e-21
                           , 8.457767819436413e-9
                           , 1.8499798003293137e-15
                           , 3.875138636756722e-7
                           ]
             , infNumIter         = 14
             , infStopReason      = SmallE_2
             , infNumFuncEvals    = 23
             , infNumJacobEvals   = 14
             , infNumLinSysSolved = 14
             }
      ,     (Infinity              ::: Infinity              ::: Nil)
        ::: (6.94996061456946e-310 ::: 6.94996061456946e-310 ::: Nil)
        ::: Nil
      )
-}

--------------------------------------------------------------------------------
-- Hock - Schittkowski MODIFIED problem 21 (box constrained),
-- minimum at (2.0, 0.0)
--
-- constr1: 2 <= p0 <=50
-- constr2: -50 <= p1 <=50
--
-- Original HS21 has the additional constraint 10*p0 - p1 >= 10
-- which is inactive at the solution, so it is dropped here.

hs21 :: Model N2 Double
hs21 p0 p1 = [p0 / 10.0, p1]

hs21_jac :: Jacobian N2 Double
hs21_jac _ _ = [ 0.1 ::: 0.0 ::: Nil
               , 0.0 ::: 1.0 ::: Nil
               ]

hs21_n :: Int
hs21_n = 2

hs21_params :: SizedList N2 Double
hs21_params = -1.0 ::: -1.0 ::: Nil

hs21_samples :: [Double]
hs21_samples = replicate hs21_n 0.0

hs21_lb, hs21_ub :: SizedList N2 Double
hs21_lb = 2.0  ::: -50.0 ::: Nil
hs21_ub = 50.0 :::  50.0 ::: Nil

run_hs21 :: Result N2
run_hs21 = levmar hs21
                  (Just hs21_jac)
                  hs21_params
                  hs21_samples
                  1000
                  opts
                  (Just hs21_lb)
                  (Just hs21_ub)
                  noLinearConstraints
                  Nothing

{- run_hs21 =>
Right ( 2.0 ::: -4.688186073682305e-19 ::: Nil
      , Info { infValues = [ 1.04
                           , 4.000000000000001e-2
                           , 4.688186073682305e-19
                           , 0.0
                           , 1.7333333336800002e-4
                           ]
             , infNumIter         = 5
             , infStopReason      = SmallGradient
             , infNumFuncEvals    = 10
             , infNumJacobEvals   = 6
             , infNumLinSysSolved = 5
             }
      ,     (Infinity              ::: -Infinity              ::: Nil)
        ::: (6.94996061456946e-310 ::: 1.390501401903724e-309 ::: Nil)
        ::: Nil
      )
-}

--------------------------------------------------------------------------------
-- Problem hatfldb (box constrained),
-- minimum at (0.947214, 0.8, 0.64, 0.4096)
--
-- constri: pi >= 0.0 (i=1..4)
-- constr5: p1 <= 0.8

hatfldb :: Model N4 Double
hatfldb p0 p1 p2 p3 = [ p0 - 1.0
                      , p0 - sqrt p1
                      , p1 - sqrt p2
                      , p2 - sqrt p3
                      ]

hatfldb_jac :: Jacobian N4 Double
hatfldb_jac _ p1 p2 p3 = [ 1.0 ::: 0.0            ::: 0.0            ::: 0.0            ::: Nil
                         , 1.0 ::: -0.5 / sqrt p1 ::: 0.0            ::: 0.0            ::: Nil
                         , 0.0 ::: 1.0            ::: -0.5 / sqrt p2 ::: 0.0            ::: Nil
                         , 0.0 ::: 0.0            ::: 1.0            ::: -0.5 / sqrt p3 ::: Nil
                         ]

hatfldb_n :: Int
hatfldb_n = 4

hatfldb_params :: SizedList N4 Double
hatfldb_params = 0.1 ::: 0.1 ::: 0.1 ::: 0.1 ::: Nil

hatfldb_samples :: [Double]
hatfldb_samples = replicate hatfldb_n 0.0

hatfldb_lb, hatfldb_ub :: SizedList N4 Double
hatfldb_lb = 0.0      ::: 0.0 ::: 0.0      ::: 0.0      ::: Nil
hatfldb_ub = _DBL_MAX ::: 0.8 ::: _DBL_MAX ::: _DBL_MAX ::: Nil

run_hatfldb :: Result N4
run_hatfldb = levmar hatfldb
                     (Just hatfldb_jac)
                     hatfldb_params
                     hatfldb_samples
                     1000
                     opts
                     (Just hatfldb_lb)
                     (Just hatfldb_ub)
                     noLinearConstraints
                     Nothing

{- run_hatfldb =>
Right ( 0.9472135954999581 ::: 0.8 ::: 0.6400000000000067 ::: 0.40960000000001073 ::: Nil
      , Info { infValues = [ 0.9502633403898972
                           , 5.572809000084124e-3
                           , 1.3010426069825884e-15
                           , 2.0738413641161756e-30
                           , 9.460108956138705e-5
                           ]
             , infNumIter         = 938
             , infStopReason      = SmallDp
             , infNumFuncEvals    = 3151
             , infNumJacobEvals   = 938
             , infNumLinSysSolved = 938
             }
      ,     (Infinity        ::: Infinity        ::: Infinity        ::: Infinity        ::: Nil)
        ::: (1.06047983e-314 ::: 1.06047983e-314 ::: 1.06047983e-314 ::: 1.06047983e-314 ::: Nil)
        ::: (Infinity        ::: Infinity        ::: Infinity        ::: Infinity        ::: Nil)
        ::: (1.06047983e-314 ::: 1.06047983e-314 ::: 1.06047983e-314 ::: 1.06047983e-314 ::: Nil)
        ::: Nil
      )
-}

--------------------------------------------------------------------------------
-- Problem hatfldc (box constrained),
-- minimum at (1.0, 1.0, 1.0, 1.0)
--
-- constri:   pi >= 0.0  (i=1..4)
-- constri+4: pi <= 10.0 (i=1..4)

hatfldc :: Model N4 Double
hatfldc p0 p1 p2 p3 = [ p0 - 1.0
                      , p0 - sqrt p1
                      , p1 - sqrt p2
                      , p3 - 1.0
                      ]

hatfldc_jac :: Jacobian N4 Double
hatfldc_jac _ p1 p2 _ = [ 1.0 ::: 0.0            ::: 0.0            ::: 0.0 ::: Nil
                        , 1.0 ::: -0.5 / sqrt p1 ::: 0.0            ::: 0.0 ::: Nil
                        , 0.0 ::: 1.0            ::: -0.5 / sqrt p2 ::: 0.0 ::: Nil
                        , 0.0 ::: 0.0            ::: 0.0            ::: 1.0 ::: Nil
                        ]

hatfldc_n :: Int
hatfldc_n = 4

hatfldc_params :: SizedList N4 Double
hatfldc_params = 0.9 ::: 0.9 ::: 0.9 ::: 0.9 ::: Nil

hatfldc_samples :: [Double]
hatfldc_samples = replicate hatfldc_n 0.0

hatfldc_lb, hatfldc_ub :: SizedList N4 Double
hatfldc_lb =  0.0 :::  0.0 :::  0.0 :::  0.0 ::: Nil
hatfldc_ub = 10.0 ::: 10.0 ::: 10.0 ::: 10.0 ::: Nil

run_hatfldc :: Result N4
run_hatfldc = levmar hatfldc
                     (Just hatfldc_jac)
                     hatfldc_params
                     hatfldc_samples
                     1000
                     opts
                     (Just hatfldc_lb)
                     (Just hatfldc_ub)
                     noLinearConstraints
                     Nothing

{- run_hatfldc =>
Right ( 0.9999999999999974 ::: 0.99999999999999 ::: 0.9999999999999779 ::: 1.0 ::: Nil
      , Info { infValues = [ 2.47401270181503e-2
                           , 1.3718784179859158e-29
                           , 5.11353737418956e-11
                           , 1.4884674673286283e-18
                           , 7.635841678377282e-8
                           ]
             , infNumIter = 4
             , infStopReason = SmallE_2
             , infNumFuncEvals = 5
             , infNumJacobEvals = 4
             , infNumLinSysSolved = 4
             }
      ,     (Infinity        ::: Infinity :      ::: Infinity          ::: -Infinity         ::: Nil)
        ::: (1.06047983e-314 ::: 1.06047983e-314 ::: 1.06047983e-314   ::: 2.1217367583e-314 ::: Nil)
        ::: (Infinity        ::: Infinity        ::: -Infinity         ::: Infinity          ::: Nil)
        ::: (1.06047983e-314 ::: 1.06047983e-314 ::: 2.1217367583e-314 ::: 1.06047983e-314   ::: Nil)
        ::: Nil
      )
-}

--------------------------------------------------------------------------------
-- Hock - Schittkowski (modified) problem 52 (box/linearly constrained),
-- minimum at (-0.09, 0.03, 0.25, -0.19, 0.03)
--
-- constr1: p0 + 3*p1 = 0
-- constr2: p2 +   p3 - 2*p4 = 0
-- constr3: p1 -   p4 = 0
--
-- To the above 3 constraints, we add the following 5:
-- constr4: -0.09 <= p0
-- constr5:   0.0 <= p1 <= 0.3
-- constr6:          p2 <= 0.25
-- constr7:  -0.2 <= p3 <= 0.3
-- constr8:   0.0 <= p4 <= 0.3

modhs52 :: Model N5 Double
modhs52 p0 p1 p2 p3 p4 = [ 4.0*p0 - p1
                         , p1 + p2 - 2.0
                         , p3 - 1.0
                         , p4 - 1.0
                         ]

modhs52_jac :: Jacobian N5 Double
modhs52_jac _ _ _ _ _ = [ 4.0 ::: -1.0 ::: 0.0 ::: 0.0 ::: 0.0 ::: Nil
                        , 0.0 :::  1.0 ::: 1.0 ::: 0.0 ::: 0.0 ::: Nil
                        , 0.0 :::  0.0 ::: 0.0 ::: 1.0 ::: 0.0 ::: Nil
                        , 0.0 :::  0.0 ::: 0.0 ::: 0.0 ::: 1.0 ::: Nil
                        ]

modhs52_n :: Int
modhs52_n = 4

modhs52_params :: SizedList N5 Double
modhs52_params = 2.0 ::: 2.0 ::: 2.0 ::: 2.0 ::: 2.0 ::: Nil

modhs52_samples :: [Double]
modhs52_samples = replicate modhs52_n 0.0

modhs52_linear_constraints :: LinearConstraints N3 N5 Double
modhs52_linear_constraints = (     (1.0 ::: 3.0 ::: 0.0 ::: 0.0 :::  0.0 ::: Nil)
                               ::: (0.0 ::: 0.0 ::: 1.0 ::: 1.0 ::: -2.0 ::: Nil)
                               ::: (0.0 ::: 1.0 ::: 0.0 ::: 0.0 ::: -1.0 ::: Nil)
                               ::: Nil
                             , 0.0 ::: 0.0 ::: 0.0 ::: Nil
                             )

modhs52_weights :: SizedList N5 Double
modhs52_weights = 2000.0 ::: 2000.0 ::: 2000.0 ::: 2000.0 ::: 2000.0 ::: Nil

modhs52_lb, modhs52_ub :: SizedList N5 Double
modhs52_lb = -0.09    ::: 0.0 ::: -_DBL_MAX ::: -0.2 ::: 0.0 ::: Nil
modhs52_ub = _DBL_MAX ::: 0.3 ::: 0.25      :::  0.3 ::: 0.3 ::: Nil

run_modhs52 :: Result N5
run_modhs52 = levmar modhs52
                     (Just modhs52_jac)
                     modhs52_params
                     modhs52_samples
                     1000
                     opts
                     (Just modhs52_lb)
                     (Just modhs52_ub)
                     (Just modhs52_linear_constraints)
                     (Just modhs52_weights)

{- run_modhs52 =>
Right ( -0.10726256881484843 ::: 3.575418960494949e-2 ::: 0.2715083811346546 ::: -0.2000000019247555 ::: 3.5754189604949506e-2 ::: Nil
      , Info { infValues = [ 8.295857988165704
                           , 5.451173183409142
                           , 1.1266725658680453e-8
                           , 7.001683639593823e-35
                           , 14.796555841095126
                           ]
             , infNumIter         = 37
             , infStopReason      = SmallDp
             , infNumFuncEvals    = 76
             , infNumJacobEvals   = 37
             , infNumLinSysSolved = 76
             }
      ,     (8.00906987461398e-2     ::: 5.606348912229788e-2   ::: 7.608616380883285e-2     ::: 0.0                     ::: 0.0                     ::: Nil)
        ::: (1.3300064540394536e128  ::: 1.2407670163481725e164 ::: 2.062583e-317            ::: 0.0                     ::: 1.330006443552782e128   ::: Nil)
        ::: (-2.9355119850619e-85    ::: 3.46043998001964e-310  ::: 0.0                      ::: 1.2433080910244667e86   ::: -4.2290595329703815e223 ::: Nil)
        ::: (-1.0261297670251651e-60 ::: 3.1e-322               ::: 0.0                      ::: -2.4394031186789343e-93 ::: -4.346705625329091e79   ::: Nil)
        ::: (5.28021237e-315         ::: 0.0                    ::: -1.8970981212800587e-267 ::: 6.005574375491326e-36   ::: 1.459531513804533e-246  ::: Nil)
        ::: Nil
      )
-}

--------------------------------------------------------------------------------
-- Schittkowski (modified) problem 235 (box/linearly constrained),
-- minimum at (-1.725, 2.9, 0.725)
--
-- constr1: p0 + p2 = -1.0;
--
-- To the above constraint, we add the following 2:
-- constr2: p1 - 4*p2 = 0
-- constr3: 0.1 <= p1 <= 2.9
-- constr4: 0.7 <= p2

mods235 :: Model N3 Double
mods235 p0 p1 _ = [ 0.1*(p0 - 1.0)
                  , p1 - p0**2
                  ]

mods235_jac :: Jacobian N3 Double
mods235_jac p0 _ _ = [ 0.1     ::: 0.0 ::: 0.0 ::: Nil
                     , -2.0*p0 ::: 1.0 ::: 0.0 ::: Nil
                     ]

mods235_n :: Int
mods235_n = 2

mods235_params :: SizedList N3 Double
mods235_params = -2.0 ::: 3.0 ::: 1.0 ::: Nil

mods235_samples :: [Double]
mods235_samples = replicate mods235_n 0.0

mods235_linear_constraints :: LinearConstraints N2 N3 Double
mods235_linear_constraints = (     (1.0 ::: 0.0 :::  1.0 ::: Nil)
                               ::: (0.0 ::: 1.0 ::: -4.0 ::: Nil)
                               ::: Nil
                             , -1.0 ::: 0.0 ::: Nil
                             )

mods235_lb, mods235_ub :: SizedList N3 Double
mods235_lb = -_DBL_MAX ::: 0.1 ::: 0.7      ::: Nil
mods235_ub =  _DBL_MAX ::: 2.9 ::: _DBL_MAX ::: Nil

run_mods235 :: Result N3
run_mods235 = levmar mods235
                     (Just mods235_jac)
                     mods235_params
                     mods235_samples
                     1000
                     opts
                     (Just mods235_lb)
                     (Just mods235_ub)
                     (Just mods235_linear_constraints)
                     Nothing

{- run_mods235 =>
Right ( -1.7250000000043926 ::: 2.9000000000175707 ::: 0.7250000000043927 ::: Nil
      , Info { infValues = [ 1.051873723955268e7
                           , 7.997539062493696e-2
                           , 5.005265614795585e-8
                           , 7.612840476704719e-32
                           , 5.449901349526238e-6
                           ]
             , infNumIter         = 6
             , infStopReason      = SmallDp
             , infNumFuncEvals    = 6
             , infNumJacobEvals   = 6
             , infNumLinSysSolved = 6
             }
      ,     (2.2378563418280644e-3  ::: -3.7830961156569275e-11 ::: 0.0                   ::: Nil)
        ::: (1.662705184422498e11   ::: 4.02688470445e-312      ::: 1.8108435015223048e95 ::: Nil)
        ::: (-1.7588779497923699e59 ::: 2.40017e-319            ::: -1.759139859665786e59 ::: Nil)
        ::: Nil
      )
-}

--------------------------------------------------------------------------------
-- Boggs and Tolle modified problem 7 (box/linearly constrained),
-- minimum at (0.7, 0.49, 0.19, 1.19, -0.2)
--
-- We keep the original objective function & starting point and use the
-- following constraints:
--
-- subject to cons1:
--  x[1]+x[2] - x[3] = 1.0;
-- subject to cons2:
--   x[2] - x[4] + x[1] = 0.0;
-- subject to cons3:
--   x[5] + x[1] = 0.5;
-- subject to cons4:
--   x[5]>=-0.3;
-- subject to cons5:
--    x[1]<=0.7;

modbt7 :: Model N5 Double
modbt7 p0 p1 _ _ _ = replicate modbt7_n (100.0*m**2 + n**2)
    where
      m = p1 - p0**2
      n = p0 - 1.0

modbt7_jac :: Jacobian N5 Double
modbt7_jac p0 p1 _ _ _ = replicate modbt7_n
                         (    -400.0*m*p0 + 2.0*p0 - 2.0
                           ::: 200.0*m
                           ::: 0.0
                           ::: 0.0
                           ::: 0.0
                           ::: Nil
                         )
    where
      m = p1 - p0**2

modbt7_n :: Int
modbt7_n = 5

modbt7_params :: SizedList N5 Double
modbt7_params = -2.0 ::: 1.0 ::: 1.0 ::: 1.0 ::: 1.0 ::: Nil

modbt7_samples :: [Double]
modbt7_samples = replicate modbt7_n 0.0

modbt7_linear_constraints :: LinearConstraints N3 N5 Double
modbt7_linear_constraints = (     (1.0 ::: 1.0 ::: -1.0 :::  0.0 ::: 0.0 ::: Nil)
                              ::: (1.0 ::: 1.0 :::  0.0 ::: -1.0 ::: 0.0 ::: Nil)
                              ::: (1.0 ::: 0.0 :::  0.0 :::  0.0 ::: 1.0 ::: Nil)
                              ::: Nil
                            , 1.0 ::: 0.0 ::: 0.5 ::: Nil
                            )

modbt7_lb, modbt7_ub :: SizedList N5 Double
modbt7_lb = -_DBL_MAX ::: -_DBL_MAX ::: -_DBL_MAX ::: -_DBL_MAX ::: -0.3     ::: Nil
modbt7_ub = 0.7       ::: _DBL_MAX  ::: _DBL_MAX  ::: _DBL_MAX  ::: _DBL_MAX ::: Nil

run_modbt7 :: Result N5
run_modbt7 = levmar modbt7
                     (Just modbt7_jac)
                     modbt7_params
                     modbt7_samples
                     1000
                     opts
                     (Just modbt7_lb)
                     (Just modbt7_ub)
                     (Just modbt7_linear_constraints)
                     Nothing

-- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
-- !! TODO: This returns with: infStopReason = MaxIterations !!
-- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

{- run_modbt7 =>
Right ( 0.8027784157269001 ::: 0.6437482395532486 ::: 0.44652665528014895 ::: 1.4465266552801486 ::: -0.3027784157269002 ::: Nil
      , Info { infValues = [ 113808.31671508029
                           , 7.583973304165455e-3
                           , 1.83706477847966e-2
                           , 8.73511123476937e-8
                           , 477.47907973252035
                           ]
             , infNumIter         = 1000
             , infStopReason      = MaxIterations
             , infNumFuncEvals    = 1290
             , infNumJacobEvals   = 1000
             , infNumLinSysSolved = 1289
             }
      ,     (1.416096391982882e-3   ::: 1.8524500981821651e-3 ::: 0.0                    ::: 0.0                     ::: 0.0                    ::: Nil)
        ::: (-3.855852680731862e-14 ::: 2.0518156e-317        ::: 0.0                    ::: 0.0                     ::: -3.855852678112155e-14 ::: Nil)
        ::: (3.4423759570554e-310   ::: 0.0                   ::: 0.0                    ::: -2.8969826546932408e147 ::: 8.276140723879028e296  ::: Nil)
        ::: (3.1e-322               ::: 0.0                   ::: 0.0                    ::: 3.2140095318438895e146  ::: 1.0027299938276747e173 ::: Nil)
        ::: (0.0                    ::: 0.0                   ::: 4.817154016378489e-127 ::: 1.0252290669946354e212  ::: 5.25442997e-315        ::: Nil)
        ::: Nil
      )
-}

--------------------------------------------------------------------------------
-- Equilibrium combustion problem, constrained nonlinear equation from the book
-- by Floudas et al.
--
-- Minimum at (0.0034, 31.3265, 0.0684, 0.8595, 0.0370)
--
-- constri:   pi>=0.0001 (i=1..5)
-- constri+5: pi<=100.0  (i=1..5)

combust :: Model N5 Double
combust p0 p1 p2 p3 p4 =
    [ p0*p1 + p0 - 3*p4
    , 2*p0*p1 + p0 + 3*r10*p1*p1 + p1*p2*p2 + r7*p1*p2 + r9*p1*p3 + r8*p1 - r*p4
    , 2*p1*p2*p2 + r7*p1*p2 + 2*r5*p2*p2 + r6*p2-8*p4
    , r9*p1*p3 + 2*p3*p3 - 4*r*p4
    , p0*p1 + p0 + r10*p1*p1 + p1*p2*p2 + r7*p1*p2 + r9*p1*p3 + r8*p1 + r5*p2*p2 + r6*p2 + p3*p3 - 1.0
    ]

r, r5, r6, r7, r8, r9, r10 :: Double
r   = 10
r5  = 0.193
r6  = 4.10622*1e-4
r7  = 5.45177*1e-4
r8  = 4.4975 *1e-7
r9  = 3.40735*1e-5
r10 = 9.615  *1e-7

combust_jac :: Jacobian N5 Double
combust_jac p0 p1 p2 p3 _ =
    [     p1 + 1
      ::: p0
      ::: 0.0
      ::: 0.0
      ::: -3
      ::: Nil
    ,     2*p1 + 1
      ::: 2*p0 + 6*r10*p1 + p2*p2 + r7*p2 + r9*p3 + r8
      ::: 2*p1*p2 + r7*p1
      ::: r9*p1
      ::: -r
      ::: Nil
    ,     0.0
      ::: 2*p2*p2 + r7*p2
      ::: 4*p1*p2 + r7*p1 + 4*r5*p2 + r6
      ::: 0.0
      ::: -8
      ::: Nil
    ,     0.0
      ::: r9*p3
      ::: 0.0
      ::: r9*p1 + 4*p3
      ::: -4*r
      ::: Nil
    ,     p1 + 1
      ::: p0 + 2*r10*p1 + p2*p2 + r7*p2 + r9*p3 + r8
      ::: 2*p1*p2 + r7*p1 + 2*r5*p2 + r6
      ::: r9*p1 + 2*p3
      ::: 0.0
      ::: Nil
    ]

combust_n :: Int
combust_n = 5

combust_params :: SizedList N5 Double
combust_params = 0.0001 ::: 0.0001 ::: 0.0001 ::: 0.0001 ::: 0.0001 ::: Nil

combust_samples :: [Double]
combust_samples = replicate combust_n 0.0

combust_lb, combust_ub :: SizedList N5 Double
combust_lb =   0.0001 :::   0.0001 :::   0.0001 :::   0.0001 :::   0.0001 ::: Nil
combust_ub = 100.0    ::: 100.0    ::: 100.0    ::: 100.0    ::: 100.0    ::: Nil

run_combust :: Result N5
run_combust = levmar combust
                     (Just combust_jac)
                     combust_params
                     combust_samples
                     1000
                     opts
                     (Just combust_lb)
                     (Just combust_ub)
                     noLinearConstraints
                     Nothing

{- run_combust =>
Right ( 3.430230156696717e-3 ::: 31.326496799034423 ::: 6.835040137781885e-2 ::: 0.8595289964750606 ::: 3.696244139315891e-2 ::: Nil
      , Info { infValues = [ 0.999817373652678
                           , 1.1459395059856032e-24
                           , 5.507913614570493e-10
                           , 2.6959712108107423e-12
                           , 5.541757943902946e-15
                           ]
             , infNumIter         = 68
             , infStopReason      = SmallE_2
             , infNumFuncEvals    = 87
             , infNumJacobEvals   = 68
             , infNumLinSysSolved = 68
             }
      ,     (Infinity               ::: -Infinity               ::: Infinity               ::: -Infinity              ::: -Infinity              ::: Nil)
        ::: (4.1424993e-317         ::: 8.2870224e-317          ::: 4.1424993e-317         ::: 8.2870224e-317         ::: 8.2870224e-317         ::: Nil)
        ::: (1.390331642240446e-309 ::: 6.94996061456946e-310   ::: 1.390331642240446e-309 ::: 1.390331642240446e-309 ::: 1.390331642240446e-309 ::: Nil)
        ::: (-3.105036184601594e231 ::: -3.1050361846015054e231 ::: -3.105036184601594e231 ::: -3.105036184601594e231 ::: -3.105036184601594e231 ::: Nil)
        ::: (1.06047983e-314        ::: 2.1214777256e-314       ::: 2.1214777256e-314      ::: 2.1214777256e-314      ::: 1.06047983e-314        ::: Nil)
        ::: Nil
      )
-}

-- The End ---------------------------------------------------------------------
