import numpy as np
import pytest
import src.gauss_newton as gn
import test_euler_kinematics as tek
import test_base as tb

# center at lat 512, long 512. Distance to center is 128
OC_NA_Pole = {"lat" : 45.54,  "long" : -119.60, "omega" : 1.32 }

sample_e = [512, 384, 512, 640]
sample_n = [384, 512, 640, 512]
dist_to_center = 128.0

def test_translation_east():
  #test 1 = translate east
  t_e = 1.0
  v_e = [t_e, t_e, t_e, t_e]
  v_n = [0.0, 0.0, 0.0, 0.0]
  x = gn.solve_gauss_newton_2D_transform(sample_e, sample_n, v_e, v_n)
  #print(f"east translate X: {x}\n")
  assert x['t_x'] == t_e

def test_translation_north():
  #test 2 = translate north
  t_n = 1.0
  v_e = [0.0, 0.0, 0.0, 0.0]
  v_n = [t_n, t_n, t_n, t_n]
  x = gn.solve_gauss_newton_2D_transform(sample_e, sample_n, v_e, v_n)
  #print(f"North translate X: {x}\n")
  assert x['t_y'] == t_n

def test_rotation_simple():
  #test 3 = rotate
  dtheta = 0.01 # radians 
  dp = dist_to_center * np.tan(dtheta)
  v_e = [-dp, 0.0, dp, 0.0]
  v_n = [0.0, dp, 0.0, -dp]
  x = gn.solve_gauss_newton_2D_transform(sample_e, sample_n, v_e, v_n)
  #print(f"Rotate X:  {x}")
  assert x['r'] == pytest.approx(dtheta, abs=1e-6)

def test_euler_test_quad():
  #test setup
  euler_pole = {"lat" : 45.0,  "long" : -90, "omega" : 1.23 }
  bearings  = [45.0, 135.0, 225.0, 315.0]
  sample_dist = 50000 # m

  sample_e, sample_n, v_east, v_north = tek.create_simple_sample_quad(euler_pole, bearings, sample_dist)
  #pole_result = epr.fit_euler_pole_linear(sample_lats, sample_lons, sample_v_east, sample_v_north, True)
  x = gn.solve_gauss_newton_2D_transform_geo(sample_e, sample_n, v_east, v_north, euler_pole)
  #print(f"\ntest_euler_test_quad x: {x}\n")
  assert x['r'] == pytest.approx(np.radians(euler_pole["omega"]), abs=1e-4)

def test_random_rot_disk(): 
  #test 4 = run euler kinematics random disk
  euler_pole = {"lat" : 45.0,  "long" : -90, "omega" : 1.23 }
  sample_count = 400
  sample_dist = 50000 # m
  test_omega = 1.23 
  
  sample_n, sample_e, v_east, v_north = \
      tek.create_random_sample_ring(
        euler_pole, 
        sample_count, 
        sample_dist, 
        test_omega, 
        1.0, 
        None)
  x = gn.solve_gauss_newton_2D_transform_geo(sample_e, sample_n, v_east, v_north, euler_pole)
  # print(f"test_rot_disk x: {x}\n")
  assert x['r'] == pytest.approx(np.radians(test_omega), abs=1e-4)

def test_random_cropped_rot_disk(): 
  #test 4 = run euler kinematics random disk
  euler_pole = {"lat" : 45.0,  "long" : -90, "omega" : 1.23 }
  sample_count = 400
  sample_dist = 50000 # m
  test_omega = 1.23 
  crop = 0.5 # 50% cropped out
  
  sample_n, sample_e, v_east, v_north = \
      tek.create_random_sample_ring(
        euler_pole, 
        sample_count, 
        sample_dist, 
        test_omega,
        crop,
        None)
  x = gn.solve_gauss_newton_2D_transform_geo(sample_e, sample_n, v_east, v_north, euler_pole)
  # print(f"test_rot_disk x: {x}\n")
  assert x['r'] == pytest.approx(np.radians(test_omega), abs=1e-4)

def test_using_north_rotation():
  #test setup
  euler_pole = OC_NA_Pole 
  euler_n_pole = {"lat" : 0.0,  "long": OC_NA_Pole['long'] + 90, "omega" : 1.43 }
  sample_count = 400
  sample_dist = 50000 # m
  test_omega = 0.001
  crop = 0.5 # 50% cropped out

  R = 6371.0E3 # Earth radius in m
  north_v_for_pole_rot= np.sin(np.radians(test_omega)) * R

  sample_n, sample_e, v_east, v_north = \
    tek.create_random_sample_ring(euler_pole, 
                              sample_count, 
                              sample_dist, 
                              test_omega, 
                              crop,
                              euler_n_pole)
  
  x = gn.solve_gauss_newton_2D_transform_geo(sample_e, sample_n, v_east, v_north, euler_pole)
  
  # print(f"North translate X: {x['t_y']}\n")
  assert x['t_y'] == pytest.approx(north_v_for_pole_rot, abs=1e-2)

def test_against_pnw_GPS_data():
  euler_pole = {"lat" : 45.0,  "long" : -118, "omega" : 1.23 }
  center_lat = 45.0
  center_long = -118
  max_distance = 550000 # m
  sample_lats, sample_lons, sample_v_east, sample_v_north = tb.get_GPS_rotation_data(center_lat, center_long, max_distance)

  #pole_result = epr.fit_euler_pole_linear(sample_lats, sample_lons, sample_v_east, sample_v_north)
  x = gn.solve_gauss_newton_2D_transform_geo(sample_lons, sample_lats, sample_v_east, sample_v_north, euler_pole)
  print(f"samples: {len(sample_lats)}")
  gn.print_x(x)
  #ek.print_result ("test_GPS_pole_extraction", pole_result)

