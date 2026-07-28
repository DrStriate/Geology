import numpy as np
import pytest
import test_base as tb
import euler_pole_regression as epr
import euler_kinematics as ek
import test_utils as tu
import geo_helper as gh
from geo_helper import PLoc, PAvel, EulerPole

OC_NA_Pole = EulerPole(-119.60, 45.54, 1.32)
SeattlePloc = PLoc(-122.3321, 47.6062)

def test_euler_pole_from_quad():
  #test setup
  euler_pole = OC_NA_Pole #{"lat" : 45.0,  "long" : -90, "omega" : 1.23 }
  azimuths  = [45.0, 135.0, 225.0, 315.0]
  sample_dist = 50000 # m

  sample_lons, sample_lats, sample_v_east, sample_v_north = tu.create_simple_sample_quad(euler_pole, azimuths, sample_dist)
  pole_result = epr.fit_euler_pole_linear(sample_lats, sample_lons, sample_v_east, sample_v_north, True)
  # ek.print_result ("test_euler_pole_from_quad", pole_result)

  assert pole_result.omega == pytest.approx(euler_pole.omega)
  assert pole_result.lat == pytest.approx(euler_pole.lat)
  assert pole_result.long == pytest.approx(euler_pole.long)

def test_euler_pole_from_random_disk():
  #test setup
  euler_pole = EulerPole(-90, 45.0, 1.23)
  sample_count = 400
  diam = 550 # km
  crop = 1.0 # no crop
  test_omega = 1.23

  sample_lats, sample_lons, sample_v_east, sample_v_north = \
    tu.create_random_sample_ring(euler_pole, sample_count, diam * 1000, test_omega, crop)
  
  pole_result = epr.fit_euler_pole_linear(sample_lats, sample_lons, sample_v_east, sample_v_north, True)
  # ek.print_result ("test_euler_pole_from_random_disk", pole_result)

  assert pole_result.omega == pytest.approx(test_omega)
  assert pole_result.long == pytest.approx(euler_pole.long)
  assert pole_result.lat == pytest.approx(euler_pole.lat)

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

  assert pole_result.omega == pytest.approx(test_omega)
  assert pole_result.long == pytest.approx(euler_pole.long)
  assert pole_result.lat == pytest.approx(euler_pole.lat)

def test_euler_pole_using_north_rotation():
    #test setup
  euler_pole = OC_NA_Pole 
  euler_n_pole = EulerPole( OC_NA_Pole.long + 90, 0.0, 1.43)
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

  assert pole_result.long == pytest.approx(euler_n_pole.long)
  assert pole_result.lat == pytest.approx(euler_n_pole.lat)
  assert pole_result.omega == pytest.approx(test_omega)

# errors showing up in breaking up rotations
def test_getPoleRotationOfPoint():
  test_pole = OC_NA_Pole
  test_loc = SeattlePloc

  target = ek.getPoleRotationOfPoint(test_pole, test_loc, 20.0)[0]
  target_midpoint = ek.getPoleRotationOfPoint(test_pole, test_loc, 10.0)[0]
  target2 = ek.getPoleRotationOfPoint(test_pole, target_midpoint, 10.0)[0]
  assert target == pytest.approx(target2)

def test_euler_GPS_pole_extraction():
  center_lat = 45.0
  center_long = -118
  max_distance = 600000 # m
  lats, lons, v_easts, v_norths, s_e, s_n = tu.get_GPS_rotation_data(center_long, center_lat, max_distance)

  pole_result = epr.fit_euler_pole_linear_wtd(lats, lons, v_easts, v_norths, s_e, s_n)
  #epr.print_result ("test_GPS_pole_extraction", pole_result, len(lats))

  pole_result_sb = EulerPole(-115.41237, 43.63921,  0.5512292)
  assert pole_result.long == pytest.approx(pole_result_sb.long)
  assert pole_result.lat == pytest.approx(pole_result_sb.lat)
  assert pole_result.omega == pytest.approx(pole_result_sb.omega)

def test_combined_GPS_pole_extraction():
  center_lat = 45.0
  center_long = -119
  max_distance = 600000 # m
  lats, lons, v_easts, v_norths, s_e, s_n = tu.get_GPS_rotation_data(center_long, center_lat, max_distance)

  pole_result, pAvel_result = epr.extractEulerPoleUsingCombinedRegressions(lats, lons, v_easts, v_norths, s_e, s_n)
  #epr.print_result ("test_GPS_pole_extraction", pole_result, len(lats))
  #pAvel_result.print("PAvel_result: ")

  pole_result_sb = EulerPole(-120.151301, 44.523748, 0.595350533)
  assert pole_result.long == pytest.approx(pole_result_sb.long)
  assert pole_result.lat == pytest.approx(pole_result_sb.lat)
  assert pole_result.omega == pytest.approx(pole_result_sb.omega)
  paVel_result_sb = PAvel(13.769589, 3.694107487)
  assert pAvel_result.azimuth == pytest.approx(paVel_result_sb.azimuth)
  assert pAvel_result.vel == pytest.approx(paVel_result_sb.vel)

def test_euler_pole_from_pLoc(): #test pnw scenario with northerly motion on pole
  ploc = PLoc(OC_NA_Pole.long, OC_NA_Pole.lat) #sample point loc (arbitrary)
  pAzvel = PAvel(0.0, (gh.metersPerDegree()/1000.0))    #point motion north (azimuth, speed in km/ma)
  pole = ek.getEulerPoleFromPlocAndPavel(ploc, pAzvel)
 
  ploc.print("\nploc:")
  pole.print("pole: ")

  assert ploc.long - 90 + 360 == pytest.approx(pole.long)  # translation pole 90 degrees off reference at equator
  assert pole.lat == pytest.approx(0.0, abs = 2e-6)                    # translation north on meridian has a pone on the equator
  assert pole.omega == pytest.approx(1.0, abs=1e-6)

def test_movement_from_Euler_pole(): #test inverse: map above pole back to point
  pole = gh.EulerPole(150.4, 0.0, 1.0)  
  point = PLoc(OC_NA_Pole.long, OC_NA_Pole.lat) #sample point loc
  new_point, vel = ek.getPoleRotationOfPoint(pole, point, 1.0)
  
  point.print("\npoint: ")
  new_point.print("new_point: ")
  print(f"vel: {vel:0.1f}")

  assert new_point.long == pytest.approx(point.long)  # translation pole 90 degrees off reference
  assert new_point.lat == pytest.approx(point.lat + 1, abs = 2e-6) # 1 degree shift north
  assert vel == pytest.approx(gh.metersPerDegree()) # distance (meters) for 1 degree lat movement

def test_3_pole_50ma_yhs_movement():
  yhsLoc0Ma = PLoc(-110.67, 44.43 )
  pnwRotPole, pnwVPavel = tu.getPnwGpsRotPoleAndVelocity()
  naPAvel = PAvel(241.0, 23.0) # degrees, mm / yr
  ma = -50.0 
  
  yhsLoc50Ma = ek.getPlocFromPoleData(naPAvel, pnwRotPole, pnwVPavel, yhsLoc0Ma, ma)

  ploc50maSB = PLoc(-124.34635456003335, 41.52218800247805)
  assert yhsLoc50Ma.lat == pytest.approx(ploc50maSB.lat)
  assert yhsLoc50Ma.long == pytest.approx(ploc50maSB.long)