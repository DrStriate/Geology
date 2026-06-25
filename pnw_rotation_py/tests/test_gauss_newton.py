import numpy as np
import pytest
import src.gauss_newton_2d.gauss_newton as gn
import test_euler_kinematics as tek

longs = [512, 384, 512, 640]
lats = [384, 512, 640, 512]

def test_translation_east():
  #test 1 = translate east
  t_e = 1.0
  v_e = [t_e, t_e, t_e, t_e]
  v_n = [0.0, 0.0, 0.0, 0.0]
  x = gn.solve_gauss_newton_2D_transform(lats, longs, v_n, v_e)
  #print(f"east translate X: {x}\n")
  assert x['t_x'] == t_e

def test_translation_north():
  #test 2 = translate north
  t_n = 1.0
  v_e = [0.0, 0.0, 0.0, 0.0]
  v_n = [t_n, t_n, t_n, t_n]
  x = gn.solve_gauss_newton_2D_transform(lats, longs, v_n, v_e)
  #print(f"North translate X: {x}\n")
  assert x['t_y'] == t_n

def test_rotation():
  #test 3 = rotate
  dtheta = 0.01
  dp = 128 * np.tan(dtheta)
  v_e = [-dp, 0.0, dp, 0.0]
  v_n = [0.0, dp, 0.0, -dp]
  x = gn.solve_gauss_newton_2D_transform(lats, longs, v_n, v_e)
  #print(f"Rotate X:  {x}")
  assert x['r'] == pytest.approx(dtheta, abs=1e-6)

def test_euler_test_quad():
  #test setup
  euler_pole = {"lat" : 45.0,  "long" : -90, "omega" : 1.23 }
  bearings  = [45.0, 135.0, 225.0, 315.0]
  sample_dist = 50000 # m

  lats, lons, v_east, v_north = tek.create_simple_sample_quad(euler_pole, bearings, sample_dist)
  #pole_result = epr.fit_euler_pole_linear(sample_lats, sample_lons, sample_v_east, sample_v_north, True)
  x = gn.solve_gauss_newton_2D_transform(lats, lons, v_north, v_east)
  print(f"x: {x}")

def test_rot_disk(): # Need to reconcile gauss_newton rotations w. Euler rotations
  #test 4 = run euler kinematics random disk
  euler_pole = {"lat" : 45.0,  "long" : -90, "omega" : 1.23 }
  sample_count = 400
  sample_dist = 50000 # m
  test_omega = 1.23

  lats, lons, v_east, v_north = tek.create_random_sample_ring(euler_pole, sample_count, sample_dist, test_omega)
  x = gn.solve_gauss_newton_2D_transform(lats, lons, v_north, v_east)
  print(f"x: {x}")
