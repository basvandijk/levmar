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

ros_d :: Double
ros_d = 105.0

ros_n = 2

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

run_ros = levmar' ros
                  (Just ros_jac)
                  params
                  samples
                  1000
                  opts
                  Nothing
                  Nothing
                  noLinearConstraints
                  Nothing
    where
      params = -1.2 ::: 1.0 ::: Nil
      samples = replicate ros_n 0.0

--------------------------------------------------------------------------------

modros_lam :: Double
modros_lam = 1e02

modros_n = 3

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

run_modros = levmar' modros
                     (Just modros_jac)
                     params
                     samples
                     1000
                     opts
                     Nothing
                     Nothing
                     noLinearConstraints
                     Nothing
    where
      params = -1.2 ::: 1.0 ::: Nil
      samples = replicate modros_n 0.0

--------------------------------------------------------------------------------

powell_n = 2

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

run_powell = levmar' powell
                     (Just powell_jac)
                     params
                     samples
                     1000
                     opts
                     Nothing
                     Nothing
                     noLinearConstraints
                     Nothing
    where
      params = -1.2 ::: 1.0 ::: Nil
      samples = replicate powell_n 0.0

--------------------------------------------------------------------------------

wood_n = 6

wood :: Model' N4 Double
wood p0 p1 p2 p3 = [ 10.0 * (p1 - p0^2)
                   , 1.0 - p0
                   , sqrt 90.0 * (p3 - p2^2)
                   , 1.0 - p2
                   , sqrt 10.0 * (p1 + p3 - 2.0)
                   , (p1 - p3) / sqrt 10.0
                   ]

run_wood = levmar' wood
                   Nothing
                   params
                   samples
                   1000
                   opts
                   Nothing
                   Nothing
                   noLinearConstraints
                   Nothing
    where
      params =  -3.0 ::: -1.0 ::: -3.0 ::: -1.0 ::: Nil
      samples = replicate wood_n 0.0

--------------------------------------------------------------------------------

meyer_n = 16

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

--------------------------------------------------------------------------------
