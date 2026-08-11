import numpy as np
import pytest
import gauss_newton as gn
import test_utils as tu
from pathlib import Path
from geo_helper import EulerPole, geod, R

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

def test_against_pnw_GPS_data():
  # Absolute path of the script
  script_path = Path(__file__).resolve()
  # print(script_path)

  euler_pole = EulerPole(-118.5, 45, 0)
  center_ploc =  euler_pole.ploc()
  max_distance = 550000 # m

  lats, lons, v_easts, v_norths, s_e, s_n = \
    tu.get_GPS_rotation_data(center_ploc.long, center_ploc.lat, max_distance)

  # x = gn.solve_gauss_newton_2D_transform_geo_wtd(lons, lats, v_easts, v_norths, s_e, s_n, center_ploc)
  x = gn.solve_gauss_newton_translation_wtd(lons, lats, v_easts, v_norths, s_e, s_n, euler_pole)

  #print(f"samples: {len(lats)}")
  #gn.print_x(x)

  x_sb = np.array([1.07832394265, 3.36762990407])
  assert(x[0] == pytest.approx(x_sb[0]))
  assert(x[1] == pytest.approx(x_sb[1]))
