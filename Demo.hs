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

sqr :: Num a => a -> a
sqr x = x*x

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
opts = defaultOpts { optStopNormInfJacTe = 1e-15
                   , optStopNorm2Dp      = 1e-15
                   , optStopNorm2E       = 1e-20
                   }

--------------------------------------------------------------------------------
-- Rosenbrock function,
-- global minimum at (1, 1)

ros :: Model N2 Double
ros p0 p1 = replicate ros_n (sqr (1.0 - p0) + ros_d*sqr m)
    where
      m = p1 - sqr p0

ros_jac :: Jacobian N2 Double
ros_jac p0 p1 = replicate ros_n (p0d ::: p1d ::: Nil)
    where
      p0d = -2 + 2*p0 - 4*ros_d*m*p0
      p1d = 2*ros_d*m
      m   = p1 - sqr p0

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

--------------------------------------------------------------------------------
-- Modified Rosenbrock problem,
-- global minimum at (1, 1)

modros :: Model N2 Double
modros p0 p1 = [ 10*(p1 - sqr p0)
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

--------------------------------------------------------------------------------
-- Powell's function,
-- minimum at (0, 0)

powell :: Model N2 Double
powell p0 p1 = [ p0
               , 10.0*p0 / m + 2*sqr p1
               ]
    where
      m = p0 + 0.1

powell_jac :: Jacobian N2 Double
powell_jac p0 p1 = [ 1.0         ::: 0.0    ::: Nil
                   , 1.0 / sqr m ::: 4.0*p1 ::: Nil
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

--------------------------------------------------------------------------------
-- Wood's function,
-- minimum at (1, 1, 1, 1)

wood :: Model N4 Double
wood p0 p1 p2 p3 = [ 10.0*(p1 - sqr p0)
                   , 1.0 - p0
                   , sqrt 90.0*(p3 - sqr p2)
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

      tmp = sqr p0 + sqr p1

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
      tmp = sqr p0 + sqr p1

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

--------------------------------------------------------------------------------
-- Boggs - Tolle problem 3 (linearly constrained),
-- minimum at (-0.76744, 0.25581, 0.62791, -0.11628, 0.25581)
--
-- constr1: p0 + 3*p1      = 0
-- constr2: p2 + p3 - 2*p4 = 0
-- constr3: p1 - p4          = 0

bt3 :: Model N5 Double
bt3 p0 p1 p2 p3 p4 = replicate bt3_n ( sqr t1
                                     + sqr t2
                                     + sqr t3
                                     + sqr t4
                                     )
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

--------------------------------------------------------------------------------
-- Hock - Schittkowski problem 28 (linearly constrained),
-- minimum at (0.5, -0.5, 0.5)
--
-- constr1: p0 + 2*p1 + 3*p2 = 1

hs28 :: Model N3 Double
hs28 p0 p1 p2 = replicate hs28_n ( sqr t1
                                 + sqr t2
                                 )
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

--------------------------------------------------------------------------------
-- Hock - Schittkowski problem 48 (linearly constrained),
-- minimum at (1.0, 1.0, 1.0, 1.0, 1.0)
--
-- constr1: sum [p0, p1, p2, p3, p4] = 5
-- constr2: p2 - 2*(p3 + p4)       = -3

hs48 :: Model N5 Double
hs48 p0 p1 p2 p3 p4 = replicate hs48_n ( sqr t1
                                       + sqr t2
                                       + sqr t3
                                       )
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

--------------------------------------------------------------------------------
-- Hock - Schittkowski problem 51 (linearly constrained),
-- minimum at (1.0, 1.0, 1.0, 1.0, 1.0)
--
-- constr1: p0 + 3*p1      = 4
-- constr2: p2 + p3 - 2*p4 = 0
-- constr3: p1 - p4          = 0

hs51 :: Model N5 Double
hs51 p0 p1 p2 p3 p4 = replicate hs51_n ( sqr t1
                                       + sqr t2
                                       + sqr t3
                                       + sqr t4
                                       )
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

--------------------------------------------------------------------------------
-- Hock - Schittkowski problem 01 (box constrained),
-- minimum at (1.0, 1.0)
--
-- constr1: p1 >= -1.5

hs01 :: Model N2 Double
hs01 p0 p1 = [ 10.0*(p1 - sqr p0)
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
                  , p1 - sqr p0
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
modbt7 p0 p1 _ _ _ = replicate modbt7_n (100.0*sqr m + sqr n)
    where
      m = p1 - sqr p0
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
      m = p1 - sqr p0

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

-- The End ---------------------------------------------------------------------
