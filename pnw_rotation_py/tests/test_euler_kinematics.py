import numpy as np
import pytest
import test_base as tb
import euler_pole_regression as epr
import euler_kinematics as ek
import test_utils as tu
import geo_helper as gh
from geo_helper import PLoc, PAvel, EulerPole, R
import gauss_newton as gn
import gauss_newton_old as gno
from pyproj import Geod

OC_NA_Pole = EulerPole(-119.60, 45.54, 1.32)
SeattlePloc = PLoc(-122.3321, 47.6062)

def test_rot_pole_from_quad():
  #test setup
  gh.setGeod(realWorld = False) 
  euler_pole = OC_NA_Pole #{"lat" : 45.0,  "long" : -90, "omega" : 1.32 }
  azimuths  = [45.0, 135.0, 225.0, 315.0]
  sample_dist = 50 # m

  sample_lons, sample_lats, sample_v_east, sample_v_north = tu.create_simple_sample_quad(euler_pole, azimuths, sample_dist)
  pole_result = epr.fit_euler_pole_linear(sample_lats, sample_lons, sample_v_east, sample_v_north, True)
  # ek.print_result ("test_euler_pole_from_quad", pole_result)

  # BUG - The new geo-correct model that (now) works with GPS is pretty bad with quad test. See second version for legacy
  gn_out = gn.solve_gauss_newton_2D_transform(sample_lons, sample_lats, sample_v_east, sample_v_north, pole_result.ploc())
  gn_out2 = gn.solve_gauss_newton_translation(sample_lons, sample_lats, sample_v_east, sample_v_north, pole_result)

  assert pole_result.omega == pytest.approx(euler_pole.omega, abs=2e-3)
  assert pole_result.lat == pytest.approx(euler_pole.lat, abs=0.002)
  assert pole_result.long == pytest.approx(euler_pole.long)

def test_translation_from_quad():
  #test setup
  gh.setGeod(realWorld = False) 
  euler_pole = EulerPole(-90, 45, 0.0)
  v_trans = [1.0, 5.0] # mm/Yr
  azimuths  = [45.0, 135.0, 225.0, 315.0]
  sample_dist = 50 # m

  sample_lons, sample_lats, sample_v_east, sample_v_north = \
    tu.create_simple_sample_quad_w_trans(euler_pole, v_trans, azimuths, sample_dist)

  gn_out = gn.solve_gauss_newton_2D_transform(sample_lons, sample_lats, sample_v_east, sample_v_north, euler_pole.ploc())

  # BUG - The new geo-correct model that (now) works with GPS is pretty bad with quad test. See second version for legacy
  gn_out2 = gn.solve_gauss_newton_translation(sample_lons, sample_lats, sample_v_east, sample_v_north, euler_pole)

  assert gn_out[0] == pytest.approx(v_trans[0])
  assert gn_out[1] == pytest.approx(v_trans[1])

# Debugging test
# def test_translation_w_rot_from_quad():
#   # Now add rotation and test translation
#   gh.setGeod(realWorld = False) 
#   euler_pole = EulerPole(-90, 45, 1.0)
#   v_trans = [1.0, 5.0] # mm/Yr
#   azimuths  = [45.0, 135.0, 225.0, 315.0]
#   sample_dist = 50000 # m

#   sample_lons, sample_lats, sample_v_east, sample_v_north = \
#     tu.create_simple_sample_quad_w_trans(euler_pole, v_trans, azimuths, sample_dist)

#   gn_out = gn.solve_gauss_newton_2D_transform(sample_lons, sample_lats, sample_v_east, sample_v_north, euler_pole.ploc())

#   # BUG - The new geo-correct model that (now) works with GPS is pretty bad with quad test. See second version for legacy
#   gn_out2 = gn.getWeightedAverageValocity_geo(sample_lons, sample_lats, sample_v_east, sample_v_north, euler_pole)

#   assert gn_out[0] == pytest.approx(v_trans[0])
#   assert gn_out[1] == pytest.approx(v_trans[1])

def test_euler_pole_from_random_disk():
  #test setup
  gh.setGeod(realWorld = False)
  euler_pole = EulerPole(-90, 45.0, 1.23)
  sample_count = 400
  diam = 550 # km
  crop = 1.0 # no crop
  test_omega = 1.23

  sample_lats, sample_lons, sample_v_east, sample_v_north = \
    tu.create_random_sample_ring(euler_pole, euler_pole.ploc(), sample_count, diam, test_omega, crop)
  
  pole_result = epr.fit_euler_pole_linear(sample_lats, sample_lons, sample_v_east, sample_v_north, True)
  # ek.print_result ("test_euler_pole_from_random_disk", pole_result)

   # BUG - The new geo-correct model that (now) works with GPS is pretty bad with quad test. See second version for legacy
  gn_out = gn.solve_gauss_newton_2D_transform(sample_lons, sample_lats, sample_v_east, sample_v_north, pole_result.ploc())
  gn_out2 = gn.solve_gauss_newton_translation(sample_lons, sample_lats, sample_v_east, sample_v_north, pole_result)
  #gn_old = gno.solve_gauss_newton_2D_transform(sample_lons, sample_lats, sample_v_east, sample_v_north, pole_result.ploc())

  assert pole_result.omega == pytest.approx(test_omega)
  assert pole_result.long == pytest.approx(euler_pole.long)
  assert pole_result.lat == pytest.approx(euler_pole.lat)

# show that getting the v-pole from a point and a motion can be reversed 
def test_v_pole_from_sample_point():
  gh.setGeod(realWorld = False)
  rotPole = EulerPole(-120.1, 44.427, 0.595) # typical pnw rot pole setup 
  v_in = [0.767, 3.545] # v pavel from typical calibration
  v_pavel = PAvel.from_V(v_in) 
  pnwVPole = ek.getEulerPoleFromPlocAndPavel(rotPole.ploc(), v_pavel)

  sample_ploc = rotPole.ploc()
  v_out = ek.calculate_v_from_EulerPole(pnwVPole, sample_ploc)
  assert v_out[0] == pytest.approx(v_in[0])
  assert v_out[1] == pytest.approx(v_in[1])

def test_v_pole_from_random_disk():
  #test setup
  gh.setGeod(realWorld = False)
  sample_count = 400
  diam = 550 # km
  crop = 1.0 # no crop
  rotPole = EulerPole(-120.1, 44.427, 0.595) # typical pnw rot pole setup 
  v_in = [0.767, 3.545] # v pavel from typical calibration
  v_pavel = PAvel.from_V(v_in) 
  pnwVPole = ek.getEulerPoleFromPlocAndPavel(rotPole.ploc(), v_pavel)
  
  sample_lats, sample_lons, sample_v_east, sample_v_north = \
    tu.create_random_sample_ring(pnwVPole, rotPole.ploc(), sample_count, diam, pnwVPole.omega, crop)
  
  # BUG - The new geo-correct model that (now) works with GPS is pretty bad with quad test. See second version for legacy
  gn_out = gn.solve_gauss_newton_2D_transform(sample_lons, sample_lats, sample_v_east, sample_v_north, rotPole.ploc())
  gn_out2 = gn.solve_gauss_newton_translation(sample_lons, sample_lats, sample_v_east, sample_v_north, rotPole)

  assert gn_out == pytest.approx(v_in, abs=0.005)
  assert gn_out2 == pytest.approx(v_in, abs=0.065)

def test_euler_pole_from_random_cropped_disk():
  #test setup
  gh.setGeod(realWorld = False)
  euler_pole = OC_NA_Pole #{"lat" : 45.0,  "long" : -90, "omega" : 1.23 }
  sample_count = 400
  sample_dist = 500 # km
  test_omega = 1.23
  crop = 0.5 # 50% cropped out

  sample_lats, sample_lons, sample_v_east, sample_v_north = \
    tu.create_random_sample_ring(
      euler_pole, 
      euler_pole.ploc(),
      sample_count, 
      sample_dist, 
      test_omega, 
      crop)
  pole_result = epr.fit_euler_pole_linear(sample_lats, sample_lons, sample_v_east, sample_v_north, True)
  # ek.print_result ("test_euler_pole_from_random_cropped_disk", pole_result)

  assert pole_result.omega == pytest.approx(test_omega)
  assert pole_result.long == pytest.approx(euler_pole.long)
  assert pole_result.lat == pytest.approx(euler_pole.lat)

def test_euler_pole_using_north_rotation():
    #test setup
  gh.setGeod(realWorld = False)
  euler_pole = OC_NA_Pole 
  euler_n_pole = EulerPole( OC_NA_Pole.long + 90, 0.0, 1.43)
  sample_count = 400
  sample_dist = 5000 # km
  test_omega = 1.23
  crop = 0.5 # 50% cropped out
  sample_lats, sample_lons, sample_v_east, sample_v_north = \
    tu.create_random_sample_ring(euler_n_pole, 
                              euler_pole.ploc(),
                              sample_count, 
                              sample_dist, 
                              test_omega, 
                              crop)

  pole_result = epr.fit_euler_pole_linear(sample_lats, sample_lons, sample_v_east, sample_v_north, True)
  #ek.print_result ("test_euler_pole_using_north_rotation", pole_result)

  if pole_result.long < 0:
    assert pole_result.long == pytest.approx(-29.600000000000488)  # euler_n_pole.long +/- 180
  else:
    assert pole_result.long == pytest.approx(180.0 -29.600000000000488)
  assert pole_result.lat == pytest.approx(euler_n_pole.lat)
  assert pole_result.omega == pytest.approx(test_omega)

# errors showing up in breaking up rotations
def test_getPoleRotationOfPoint():
  gh.setGeod(realWorld = False)
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

  pole_result = epr.fit_euler_pole_linear(lats, lons, v_easts, v_norths, s_e, s_n)
  # epr.print_result ("test_GPS_pole_extraction", pole_result, len(lats))

  pole_result_sb = EulerPole(-115.41237, 43.63921,  0.5512292)
  assert pole_result.long == pytest.approx(pole_result_sb.long)
  assert pole_result.lat == pytest.approx(pole_result_sb.lat)
  assert pole_result.omega == pytest.approx(pole_result_sb.omega)

def test_combined_GPS_pole_extraction():
  gh.setGeod(realWorld = True)
  center_lat = 45.0
  center_long = -119
  max_distance = 600000 # m
  lats, lons, v_easts, v_norths, s_e, s_n = tu.get_GPS_rotation_data(center_long, center_lat, max_distance)

  pole_result, pAvel_result = epr.extractEulerPoleUsingCombinedRegressions(lats, lons, v_easts, v_norths, s_e, s_n)
  #epr.print_result ("test_GPS_pole_extraction", pole_result, len(lats))
  #pAvel_result.print("PAvel_result: ")

  pole_result_sb = EulerPole(-120.0989074197468, 44.42671931880976, 0.5949532453109426)
  assert pole_result.long == pytest.approx(pole_result_sb.long)
  assert pole_result.lat == pytest.approx(pole_result_sb.lat)
  assert pole_result.omega == pytest.approx(pole_result_sb.omega)
  paVel_result_sb = PAvel(12.210188380575712, 3.6270103641193825)
  assert pAvel_result.azimuth == pytest.approx(paVel_result_sb.azimuth)
  assert pAvel_result.vel == pytest.approx(paVel_result_sb.vel)

def test_euler_pole_from_pLoc(): #test pnw scenario with northerly motion on pole
  gh.setGeod(realWorld = False)
  ploc = PLoc(OC_NA_Pole.long, OC_NA_Pole.lat) #sample point loc (arbitrary)
  pAzvel = PAvel(0.0, gh.kmPerDegree())    #point motion north (azimuth, speed in km/ma)
  pole = ek.getEulerPoleFromPlocAndPavel(ploc, pAzvel)
 
  # ploc.print("\nploc:")
  # pole.print("pole: ")

  assert ploc.long - 90 + 360 == pytest.approx(pole.long, abs=0.01) # translation pole 90 degrees off reference at equator
  assert pole.lat == pytest.approx(0.0, abs = 0.01)                    # translation north on meridian has a pone on the equator
  assert pole.omega == pytest.approx(1.0, abs=1e-6)

def test_movement_from_Euler_pole(): #test inverse: map above pole back to point
  gh.setGeod(realWorld = False)
  pole = gh.EulerPole(150.4, 0.0, 1.0)  
  point = PLoc(OC_NA_Pole.long, OC_NA_Pole.lat) #sample point loc
  new_point, vel = ek.getPoleRotationOfPoint(pole, point, 1.0)
  
  assert new_point.long == pytest.approx(point.long)  # translation pole 90 degrees off reference
  assert new_point.lat == pytest.approx(point.lat + 1, abs = 2e-6) # 1 degree shift north
  assert vel == pytest.approx(gh.kmPerDegree()) # distance (km) for 1 degree lat movement

def test_v_pole(): #take avg velocity to create pole and then re-create the velocity from the pole
  gh.setGeod(realWorld = False)
  rotPole = EulerPole(-120.1, 44.427, 0.595) # typocal pnw rot pole setup 
  v_in = [0.767, 3.545] # v pavel from typical calibration
  v_pavel = PAvel.from_V(v_in) 

  vPole = ek.getEulerPoleFromPlocAndPavel(rotPole.ploc(), v_pavel)

  v_out = tu.calculate_v_from_Euler_pole2(vPole, rotPole.ploc(), vPole.omega, False)
  assert v_out["v_e"] == pytest.approx(v_in[0])
  assert v_out["v_n"] == pytest.approx(v_in[1])

def test_3_pole_50ma_yhs_movement():
  gh.setGeod(realWorld = True)
  yhsLoc0Ma = PLoc(-110.67, 44.43 )
  pnwRotPole, pnwVPavel = epr.getPnwGpsRotPoleAndVelocity()
  naPAvel = PAvel(241.0, 23.0) # degrees, mm / yr
  ma = -50.0 
  
  yhsLoc50Ma = ek.getPlocFromPoleData(naPAvel, pnwRotPole, pnwVPavel, yhsLoc0Ma, ma)

  ploc50maSB = PLoc(-124.3471338274746, 41.53563788930763)
  assert yhsLoc50Ma.lat == pytest.approx(ploc50maSB.lat)
  assert yhsLoc50Ma.long == pytest.approx(ploc50maSB.long)