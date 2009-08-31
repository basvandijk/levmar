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
--
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
--
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
--
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
--
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
--
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
-- Hock - Schittkowski MODIFIED problem 21 (box constrained),
-- minimum at (2.0, 0.0)
--
-- constr1: 2 <= p0 <=50
-- constr2: -50 <= p1 <=50
--
-- Original HS21 has the additional constraint 10*p0 - p1 >= 10
-- which is inactive at the solution, so it is dropped here.

hs21 :: Model' N2 Double
hs21 p0 p1 = [p0 / 10.0, p1]

hs21_jac :: Jacobian' N2 Double
hs21_jac p0 p1 = [ 0.1 ::: 0.0 ::: Nil
                 , 0.0 ::: 1.0 ::: Nil
                 ]

hs21_n = 2

hs21_params = -1.0 ::: -1.0 ::: Nil

hs21_samples = replicate hs21_n 0.0

hs21_lb = 2.0  ::: -50.0 ::: Nil
hs21_ub = 50.0 :::  50.0 ::: Nil

run_hs21 = levmar' hs21
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

hatfldb :: Model' N4 Double
hatfldb p0 p1 p2 p3 = [ p0 - 1.0
                      , p0 - sqrt p1
                      , p1 - sqrt p2
                      , p2 - sqrt p3
                      ]

hatfldb_jac :: Jacobian' N4 Double
hatfldb_jac p0 p1 p2 p3 = [ 1.0 ::: 0.0            ::: 0.0            ::: 0.0            ::: Nil
                          , 1.0 ::: -0.5 / sqrt p1 ::: 0.0            ::: 0.0            ::: Nil
                          , 0.0 ::: 1.0            ::: -0.5 / sqrt p2 ::: 0.0            ::: Nil
                          , 0.0 ::: 0.0            ::: 1.0            ::: -0.5 / sqrt p3 ::: Nil
                          ]

hatfldb_n = 4

hatfldb_params = 0.1 ::: 0.1 ::: 0.1 ::: 0.1 ::: Nil

hatfldb_samples = replicate hatfldb_n 0.0

hatfldb_lb = 0.0      ::: 0.0 ::: 0.0      ::: 0.0      ::: Nil
hatfldb_ub = _DBL_MAX ::: 0.8 ::: _DBL_MAX ::: _DBL_MAX ::: Nil

run_hatfldb = levmar' hatfldb
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

hatfldc :: Model' N4 Double
hatfldc p0 p1 p2 p3 = [ p0 - 1.0
                      , p0 - sqrt p1
                      , p1 - sqrt p2
                      , p3 - 1.0
                      ]

hatfldc_jac :: Jacobian' N4 Double
hatfldc_jac p0 p1 p2 p3 = [ 1.0 ::: 0.0            ::: 0.0            ::: 0.0 ::: Nil
                          , 1.0 ::: -0.5 / sqrt p1 ::: 0.0            ::: 0.0 ::: Nil
                          , 0.0 ::: 1.0            ::: -0.5 / sqrt p2 ::: 0.0 ::: Nil
                          , 0.0 ::: 0.0            ::: 0.0            ::: 1.0 ::: Nil
                          ]

hatfldc_n = 4

hatfldc_params = 0.9 ::: 0.9 ::: 0.9 ::: 0.9 ::: Nil

hatfldc_samples = replicate hatfldc_n 0.0

hatfldc_lb =  0.0 :::  0.0 :::  0.0 :::  0.0 ::: Nil
hatfldc_ub = 10.0 ::: 10.0 ::: 10.0 ::: 10.0 ::: Nil

run_hatfldc = levmar' hatfldc
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
-- constr1: p[0] + 3*p[1] = 0;
-- constr2: p[2] +   p[3] - 2*p[4] = 0;
-- constr3: p[1] -   p[4] = 0;
--
-- To the above 3 constraints, we add the following 5:
-- constr4: -0.09 <= p[0];
-- constr5:   0.0 <= p[1] <= 0.3;
-- constr6:          p[2] <= 0.25;
-- constr7:  -0.2 <= p[3] <= 0.3;
-- constr8:   0.0 <= p[4] <= 0.3;

modhs52 :: Model' N5 Double
modhs52 p0 p1 p2 p3 p4 = [ 4.0 * p0 - p1
                         , p1 + p2 - 2.0
                         , p3 - 1.0
                         , p4 - 1.0
                         ]

modhs52_jac :: Jacobian' N5 Double
modhs52_jac p0 p1 p2 p3 p4 = [ 4.0 ::: -1.0 ::: 0.0 ::: 0.0 ::: 0.0 ::: Nil
                             , 0.0 :::  1.0 ::: 1.0 ::: 0.0 ::: 0.0 ::: Nil
                             , 0.0 :::  0.0 ::: 0.0 ::: 1.0 ::: 0.0 ::: Nil
                             , 0.0 :::  0.0 ::: 0.0 ::: 0.0 ::: 1.0 ::: Nil
                             ]

modhs52_n = 4

modhs52_params = 2.0 ::: 2.0 ::: 2.0 ::: 2.0 ::: 2.0 ::: Nil

modhs52_samples = replicate modhs52_n 0.0

modhs52_lb = -0.09    ::: 0.0 ::: -_DBL_MAX ::: -0.2 ::: 0.0 ::: Nil
modhs52_ub = _DBL_MAX ::: 0.3 ::: 0.25      :::  0.3 ::: 0.3 ::: Nil

modhs52_linear_constraints :: LinearConstraints N3 N5 Double
modhs52_linear_constraints = (     (1.0 ::: 3.0 ::: 0.0 ::: 0.0 :::  0.0 ::: Nil)
                               ::: (0.0 ::: 0.0 ::: 1.0 ::: 1.0 ::: -2.0 ::: Nil)
                               ::: (0.0 ::: 1.0 ::: 0.0 ::: 0.0 ::: -1.0 ::: Nil)
                               ::: Nil
                             , 0.0 ::: 0.0 ::: 0.0 ::: Nil
                             )

modhs52_weights = 2000.0 ::: 2000.0 ::: 2000.0 ::: 2000.0 ::: 2000.0 ::: Nil

run_modhs52 = levmar' modhs52
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

-- mods235

--------------------------------------------------------------------------------

-- modbt7

--------------------------------------------------------------------------------

-- combust

--------------------------------------------------------------------------------
