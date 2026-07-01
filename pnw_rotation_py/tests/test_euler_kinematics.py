import numpy as np
import pytest
import test_base as tb
import euler_pole_regression as epr
import euler_kinematics as ek
import test_utils as tu

OC_NA_Pole = {"lat" : 45.54,  "long" : -119.60, "omega" : 1.32 }

def test_euler_pole_from_quad():
  #test setup
  euler_pole = OC_NA_Pole #{"lat" : 45.0,  "long" : -90, "omega" : 1.23 }
  bearings  = [45.0, 135.0, 225.0, 315.0]
  sample_dist = 50000 # m

  sample_lons, sample_lats, sample_v_east, sample_v_north = tu.create_simple_sample_quad(euler_pole, bearings, sample_dist)
  pole_result = epr.fit_euler_pole_linear(sample_lats, sample_lons, sample_v_east, sample_v_north, True)
  # ek.print_result ("test_euler_pole_from_quad", pole_result)

  assert pole_result['omega'] == pytest.approx(euler_pole['omega'])
  assert pole_result['lat'] == pytest.approx(euler_pole['lat'])
  assert pole_result['long'] == pytest.approx(euler_pole['long'])

def test_euler_pole_from_random_disk():
  #test setup
  euler_pole = {"lat" : 45.0,  "long" : -90, "omega" : 1.23 }
  sample_count = 400
  diam = 550 # km
  crop = 1.0 # no crop
  test_omega = 1.23

  sample_lats, sample_lons, sample_v_east, sample_v_north = \
    tu.create_random_sample_ring(euler_pole, sample_count, diam * 1000, test_omega, crop)
  
  pole_result = epr.fit_euler_pole_linear(sample_lats, sample_lons, sample_v_east, sample_v_north, True)
  # ek.print_result ("test_euler_pole_from_random_disk", pole_result)

  assert pole_result['omega'] == pytest.approx(test_omega)
  assert pole_result['lat'] == pytest.approx(euler_pole['lat'])
  assert pole_result['long'] == pytest.approx(euler_pole['long'])

def test_euler_pole_from_random_cropped_disk():
  #test setup
  euler_pole = OC_NA_Pole #{"lat" : 45.0,  "long" : -90, "omega" : 1.23 }
  sample_count = 400
  sample_dist = 50000 # m
  test_omega = 1.23
  crop = 0.5 # 50% cropped out

  sample_lats, sample_lons, sample_v_east, sample_v_north = \
    tu.create_random_sample_ring(
      euler_pole, 
      sample_count, 
      sample_dist, 
      test_omega, 
      crop,
      None)
  pole_result = epr.fit_euler_pole_linear(sample_lats, sample_lons, sample_v_east, sample_v_north, True)
  # ek.print_result ("test_euler_pole_from_random_cropped_disk", pole_result)

  assert pole_result['omega'] == pytest.approx(test_omega)
  assert pole_result['lat'] == pytest.approx(euler_pole['lat'])
  assert pole_result['long'] == pytest.approx(euler_pole['long'])

def test_euler_pole_using_north_rotation():
    #test setup
  euler_pole = OC_NA_Pole 
  euler_n_pole = {"lat" : 0.0,  "long": OC_NA_Pole['long'] + 90, "omega" : 1.43 }
  sample_count = 400
  sample_dist = 50000 # m
  test_omega = 1.23
  crop = 0.5 # 50% cropped out
  sample_lats, sample_lons, sample_v_east, sample_v_north = \
    tu.create_random_sample_ring(euler_pole, 
                              sample_count, 
                              sample_dist, 
                              test_omega, 
                              crop,
                              euler_n_pole)

  pole_result = epr.fit_euler_pole_linear(sample_lats, sample_lons, sample_v_east, sample_v_north, True)
  #ek.print_result ("test_euler_pole_using_north_rotation", pole_result)

  assert pole_result['omega'] == pytest.approx(test_omega)
  assert pole_result['lat'] == pytest.approx(euler_n_pole['lat'])
  assert pole_result['long'] == pytest.approx(euler_n_pole['long'])

def test_GPS_pole_extraction():
  center_lat = 45.0
  center_long = -118
  max_distance = 550000 # m
  lats, lons, v_easts, v_norths, s_e, s_n = tu.get_GPS_rotation_data(center_lat, center_long, max_distance)

  pole_result = epr.fit_euler_pole_linear(lats, lons, v_easts, v_norths)
  epr.print_result ("test_GPS_pole_extraction", pole_result, len(lats))
  
  