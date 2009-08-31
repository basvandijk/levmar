module Demo where

import LevMar

--------------------------------------------------------------------------------

type N0 = Z
type N1 = S N0
type N2 = S N1
type N3 = S N2
type N4 = S N3
type N5 = S N4

--------------------------------------------------------------------------------

opts :: Options Double
opts = defaultOpts { optEpsilon1 = 1e-15
                   , optEpsilon2 = 1e-15
                   , optEpsilon3 = 1e-20
                   }

--------------------------------------------------------------------------------
-- Rosenbrock function,
-- global minimum at (1, 1)

ros :: Model' N2 Double
ros p0 p1 = replicate ros_n ((1.0 - p0)^2 + ros_d * m^2)
    where
      m = p1 - p0^2

ros_jac :: Jacobian' N2 Double
ros_jac p0 p1 = replicate ros_n (p0d ::: p1d ::: Nil)
    where
      p0d = -2 + 2 * p0 - 4 * ros_d * m * p0
      p1d = 2 * ros_d * m
      m   = p1 - p0^2

ros_d :: Double
ros_d = 105.0

ros_n = 2

ros_params = -1.2 ::: 1.0 ::: Nil

ros_samples = replicate ros_n 0.0

run_ros = levmar' ros
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

modros :: Model' N2 Double
modros p0 p1 = [ 10 * (p1 - p0^2)
               , 1.0 - p0
               , modros_lam
               ]

modros_jac :: Jacobian' N2 Double
modros_jac p0 _ = [ -20 * p0 ::: 10.0 ::: Nil
                  , -1.0     ::: 0.0  ::: Nil
                  , 0.0      ::: 0.0  ::: Nil
                  ]

modros_lam :: Double
modros_lam = 1e02

modros_n = 3

modros_params = -1.2 ::: 1.0 ::: Nil

modros_samples = replicate modros_n 0.0

run_modros = levmar' modros
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

powell :: Model' N2 Double
powell p0 p1 = [ p0
               , 10.0 * p0 / m + 2 * p1^2
               ]
    where
      m = p0 + 0.1

powell_jac :: Jacobian' N2 Double
powell_jac p0 p1 = [ 1.0       ::: 0.0      ::: Nil
                   , 1.0 / m^2 ::: 4.0 * p1 ::: Nil
                   ]
    where
      m = p0 + 0.1

powell_n = 2

powell_params = -1.2 ::: 1.0 ::: Nil

powell_samples = replicate powell_n 0.0

run_powell = levmar' powell
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

wood :: Model' N4 Double
wood p0 p1 p2 p3 = [ 10.0 * (p1 - p0^2)
                   , 1.0 - p0
                   , sqrt 90.0 * (p3 - p2^2)
                   , 1.0 - p2
                   , sqrt 10.0 * (p1 + p3 - 2.0)
                   , (p1 - p3) / sqrt 10.0
                   ]

wood_n = 6

wood_params =  -3.0 ::: -1.0 ::: -3.0 ::: -1.0 ::: Nil

wood_samples = replicate wood_n 0.0

run_wood = levmar' wood
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
-- Meyer's (reformulated) problem,
-- minimum at (2.48, 6.18, 3.45)

meyer :: Model N3 Double Double
meyer p0 p1 p2 x = p0 * exp (10.0 * p1 / (ui + p2) - 13.0)
    where
      ui = 0.45 + 0.05 * x

meyer_jac :: Jacobian N3 Double Double
meyer_jac p0 p1 p2 x =     tmp
                       ::: 10.0 * p0 * tmp / (ui + p2)
                       ::: -10.0 * p0 * p1 * tmp / ((ui + p2) * (ui + p2))
                       ::: Nil
    where
      tmp = exp (10.0 * p1 / (ui + p2) - 13.0)
      ui = 0.45 + 0.05 * x

meyer_n = 16

meyer_params = 8.85 ::: 4.0 ::: 2.5 ::: Nil

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

run_meyer_jac = levmar meyer
                       (Just meyer_jac)
                       meyer_params
                       meyer_samples
                       1000
                       opts
                       Nothing
                       Nothing
                       noLinearConstraints
                       Nothing

run_meyer = levmar meyer
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

helval :: Model'  N3 Double
helval p0 p1 p2 = [ 10.0 * (p2 - 10.0 * theta)
                  , 10.0 * sqrt tmp - 1.0
                  , p2
                  ]
    where
      m = atan (p1 / p0) / (2.0 * pi)

      tmp = p0^2 + p1^2

      theta | p0 < 0.0  = m + 0.5
            | 0.0 < p0  = m
            | p1 >= 0   = 0.25
            | otherwise = -0.25

heval_jac :: Jacobian' N3 Double
heval_jac p0 p1 p2 = [ 50.0 * p1 / (pi * tmp) ::: -50.0 * p0 / (pi * tmp) ::: 10.0 ::: Nil
                     , 10.0 * p0 / sqrt tmp   :::  10.0 * p1 / sqrt tmp   ::: 0.0  ::: Nil
                     , 0.0                    ::: 0.0                     ::: 1.0  ::: Nil
                     ]
    where
      tmp = p0^2 + p1^2

helval_n = 3

helval_params = -1.0 ::: 0.0 ::: 0.0 ::: Nil

helval_samples = replicate helval_n 0.0

run_helval = levmar' helval
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
-- constr1: p0 + 3 * p1      = 0
-- constr2: p2 + p3 - 2 * p4 = 0
-- constr3: p1 - p4          = 0

bt3 :: Model' N5 Double
bt3 p0 p1 p2 p3 p4 = replicate bt3_n (t1^2 + t2^2 + t3^2 + t4^2)
    where
      t1 = p0 - p1
      t2 = p1 + p2 - 2.0
      t3 = p3 - 1.0
      t4 = p4 - 1.0

bt3_jac :: Jacobian' N5 Double
bt3_jac p0 p1 p2 p3 p4 = replicate bt3_n (   2.0 * t1
                                         ::: 2.0 * (t2 - t1)
                                         ::: 2.0 * t2
                                         ::: 2.0 * t3
                                         ::: 2.0 * t4
                                         ::: Nil
                                         )
    where
      t1 = p0 - p1
      t2 = p1 + p2 - 2.0
      t3 = p3 - 1.0
      t4 = p4 - 1.0

bt3_n = 5

bt3_params = 2.0 ::: 2.0 ::: 2.0 :::2.0 ::: 2.0 ::: Nil

bt3_samples = replicate bt3_n 0.0

bt3_linear_constraints :: LinearConstraints N3 N5 Double
bt3_linear_constraints = (     (1.0 ::: 3.0 ::: 0.0 ::: 0.0 :::  0.0 ::: Nil)
                           ::: (0.0 ::: 0.0 ::: 1.0 ::: 1.0 ::: -2.0 ::: Nil)
                           ::: (0.0 ::: 1.0 ::: 0.0 ::: 0.0 ::: -1.0 ::: Nil)
                           ::: Nil
                         , 0.0 ::: 0.0 ::: 0.0 ::: Nil
                         )

run_bt3 = levmar' bt3
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
-- constr1: p0 + 2 * p1 + 3 * p2 = 1

hs28 :: Model' N3 Double
hs28 p0 p1 p2 = replicate hs28_n (t1^2 + t2^2)
    where
      t1 = p0 + p1
      t2 = p1 + p2

hs28_jac :: Jacobian' N3 Double
hs28_jac p0 p1 p2 = replicate hs28_n (     2.0 * t1
                                       ::: 2.0 * (t1 + t2)
                                       ::: 2.0 * t2
                                       ::: Nil
                                     )
    where
      t1 = p0 + p1
      t2 = p1 + p2

hs28_n = 3

hs28_params = -4.0 ::: 1.0 ::: 1.0 ::: Nil

hs28_samples = replicate hs28_n 0.0

hs28_linear_constraints :: LinearConstraints N1 N3 Double
hs28_linear_constraints = ( ((1.0 ::: 2.0 ::: 3.0 ::: Nil) ::: Nil)
                          , 1.0 ::: Nil
                          )

run_hs28 = levmar' hs28
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
-- constr1: sum [p0, p1, p2, p3, p4] = 5
-- constr2: p2 - 2 * (p3 + p4)       = -3

hs48 :: Model' N5 Double
hs48 p0 p1 p2 p3 p4 = replicate hs48_n (t1^2 + t2^2 + t3^2)
    where
      t1 = p0 - 1.0
      t2 = p1 - p2
      t3 = p3 - p4

hs48_jac :: Jacobian' N5 Double
hs48_jac p0 p1 p2 p3 p4 = replicate hs48_n (     2.0 * t1
                                             ::: 2.0 * t2
                                             ::: 2.0 * t2
                                             ::: 2.0 * t3
                                             ::: 2.0 * t3
                                             ::: Nil
                                           )
    where
      t1 = p0 - 1.0
      t2 = p1 - p2
      t3 = p3 - p4

hs48_n = 3

hs48_params = 3.0 ::: 5.0 ::: -3.0 ::: 2.0 ::: -2.0 ::: Nil

hs48_samples = replicate hs48_n 0.0

hs48_linear_constraints :: LinearConstraints N2 N5 Double
hs48_linear_constraints = (     (1.0 ::: 1.0 ::: 1.0 :::  1.0 :::  1.0 ::: Nil)
                            ::: (0.0 ::: 0.0 ::: 1.0 ::: -2.0 ::: -2.0 ::: Nil)
                            ::: Nil
                          , 5.0 ::: -3.0 ::: Nil
                          )

run_hs48 = levmar' hs48
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
-- constr1: p0 + 3 * p1      = 4
-- constr2: p2 + p3 - 2 * p4 = 0
-- constr3: p1 - p4          = 0

hs51 :: Model' N5 Double
hs51 p0 p1 p2 p3 p4 = replicate hs51_n (t1^2 + t2^2 + t3^2 + t4^2)
    where
      t1 = p0 - p1
      t2 = p1 + p2 - 2.0
      t3 = p3 - 1.0
      t4 = p4 - 1.0

hs51_jac :: Jacobian' N5 Double
hs51_jac p0 p1 p2 p3 p4 = replicate hs51_n (     2.0 * t1
                                             ::: 2.0 * (t2 - t1)
                                             ::: 2.0 * t2
                                             ::: 2.0 * t3
                                             ::: 2.0 * t4
                                             ::: Nil
                                           )
    where
      t1 = p0 - p1
      t2 = p1 + p2 - 2.0
      t3 = p3 - 1.0
      t4 = p4 - 1.0

hs51_n = 5

hs51_params = 2.5 ::: 0.5 ::: 2.0 ::: -1.0 ::: 0.5 ::: Nil

hs51_samples = replicate hs51_n 0.0

hs51_linear_constraints :: LinearConstraints N3 N5 Double
hs51_linear_constraints = (     (1.0 ::: 3.0 ::: 0.0 ::: 0.0 :::  0.0 ::: Nil)
                            ::: (0.0 ::: 0.0 ::: 1.0 ::: 1.0 ::: -2.0 ::: Nil)
                            ::: (0.0 ::: 1.0 ::: 0.0 ::: 0.0 ::: -1.0 ::: Nil)
                            ::: Nil
                          , 4.0 ::: 0.0 ::: 0.0 ::: Nil
                          )

run_hs51 = levmar' hs51
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
-- constr1: p1 >= -1.5

hs01 :: Model' N2 Double
hs01 p0 p1 = [ 10.0 * (p1 - p0^2)
             , 1.0 - p0
             ]

hs01_jac :: Jacobian' N2 Double
hs01_jac p0 p1 = [ -20.0 * p0 ::: 10.0 ::: Nil
                 , -1.0       ::: 0.0  ::: Nil
                 ]

hs01_n = 2

hs01_params = -2.0 ::: 1.0 ::: Nil

hs01_samples = replicate hs01_n 0.0

hs01_lb = -_DBL_MAX ::: -1.5     ::: Nil
hs01_ub =  _DBL_MAX ::: _DBL_MAX ::: Nil

_DBL_MAX = 1e+37 -- TODO: Get this directly from <float.h>.

run_hs01 = levmar' hs01
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

-- hs21

--------------------------------------------------------------------------------

-- hatfldb

--------------------------------------------------------------------------------

-- hatfldc

--------------------------------------------------------------------------------

-- modhs52

--------------------------------------------------------------------------------

-- mods235

--------------------------------------------------------------------------------

-- modbt7

--------------------------------------------------------------------------------

-- combust

--------------------------------------------------------------------------------
