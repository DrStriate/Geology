import numpy as np
import pytest
import geopandas as gpd   
import euler_pole.euler_pole_regression as epr
import euler_pole.euler_kinematics as ek
import src.geo_helper as gh
import src.gauss_newton_2d.gauss_newton as gn
import geopandas as gpd
from shapely.geometry import box

OC_NA_Pole = {"lat" : 45.54,  "long" : -119.60, "omega" : 1.32 }

def create_simple_sample_quad(euler_pole, sample_bearings, sample_dist, return_locs_in_meters = False):
  sample_e = []
  sample_n = []
  sample_v_east = []
  sample_v_north = [] # mm/ yr

  Omega = {"omega": euler_pole['omega'], "phi": np.radians(euler_pole['lat']), "lamb": np.radians(euler_pole['long'])}
  # print("")
  for i in range(len(sample_bearings)):
    sample = ek.create_sample(euler_pole['long'], euler_pole['lat'], sample_bearings[i], sample_dist)
    p = {"phi": np.radians(sample['lat']), "lamb": np.radians(sample['lon'])}
    v = ek.calculate_v_from_Eigen_pole(Omega, p, Omega['omega']);
    if return_locs_in_meters:
      # convert sample points to meters to match v
      s_n, s_e = ek.get_northerly_easterly_from_lat_long_pts(sample['lon'], sample['lat'], euler_pole['long'], euler_pole['lat'])
      sample_e.append(s_e)
      sample_n.append(s_n) 
    else:
      sample_e.append(sample['lon'])
      sample_n.append(sample['lat'])
      
    ### BUG = not sure why I have to negate the rot result from solve_gauss.. given clockwise omega rotation 
    sample_v_east.append(v['v_e'])
    sample_v_north.append(v['v_n'])
    #print(f"{i}: sample['lat']: {sample['lat']:.3f}, sample['lon']: {sample['lon']:.3f}, v_e: {v['v_e']:.2f}  v_n: {v['v_n']:.2f}")

  if return_locs_in_meters:
    return sample_e, sample_n, sample_v_east, sample_v_north
  else:
    return sample_e, sample_n, sample_v_east, sample_v_north
  
def create_random_sample_ring(euler_pole, count, 
                              max_dist, 
                              test_omega, 
                              crop = 1.0, 
                              source_poll = None, 
                              return_locs_in_meters = False):
  rng = np.random.default_rng(seed=42)
  rands = rng.random(size=(count, 2))

  sample_n = []
  sample_e = []
  sample_v_east = []
  sample_v_north = [] # mm/ yr

  Omega = {"omega": euler_pole['omega'], "phi": np.radians(euler_pole['lat']), "lamb": np.radians(euler_pole['long'])}
  if source_poll:
    Omega_source = {"omega": source_poll['omega'], "phi": np.radians(source_poll['lat']), "lamb": np.radians(source_poll['long'])}
  else:
    Omega_source = Omega
  #print("")
  max_long =  ek.create_sample(euler_pole['long'], euler_pole['lat'], 90.0, max_dist)['lon']
  min_long =  ek.create_sample(euler_pole['long'], euler_pole['lat'], 270.0, max_dist)['lon']
  crop_long = min_long + (max_long - min_long) * crop
  cropped_samples = 0;
  for i in range(len(rands)):
    sample = ek.create_sample(euler_pole['long'], euler_pole['lat'], 360.0 * rands[i][0], max_dist * rands[i][1])
    p = {"phi": np.radians(sample['lat']), "lamb": np.radians(sample['lon'])}
    v = ek.calculate_v_from_Eigen_pole(Omega_source, p, test_omega); 

    if sample['lon'] < crop_long:
      if return_locs_in_meters:
        # convert sample points to meters to match v
        s_n, s_e = ek.get_northerly_easterly_from_lat_long_pts(sample['lon'], sample['lat'], euler_pole['long'], euler_pole['lat'])
        sample_e.append(s_e)
        sample_n.append(s_n) 
      else:
        sample_e.append(sample['lon'])
        sample_n.append(sample['lat'])
      
      sample_v_east.append(v['v_e'])
      sample_v_north.append(v['v_n'])

    # print(f"{i}: sample['lat']: {sample['lat']:.3f}, sample['lon']: {sample['lon']:.3f}, v_e: {v['v_e']:.2f}  v_n: {v['v_n']:.2f}")
    cropped_samples += 1
  #print(f"samples = {cropped_samples} out of {count}")
  
  return sample_n, sample_e, sample_v_east, sample_v_north

def get_GPS_rotation_data (center_lat, center_long, max_distance):
  gdf = gpd.read_file("zip://data/NSHM2023_GPS_velocity.zip")
  list_lats = gdf['geometry'].y.values
  list_lons = gdf['geometry'].x.values
  list_v_east = gdf['Ve'].values       
  list_v_north = gdf['Vn'].values    
  
  sample_lats = []
  sample_lons = []
  sample_v_east = []
  sample_v_north = [] # mm/ yr
  for i in range(len(list_lats)):
    dist = gh.GeoHelper.DistanceFromLatLong((list_lats[i], list_lats[i]), (center_lat, center_long))
    if dist < max_distance:
      sample_lats.append(list_lats[i])
      sample_lons.append(list_lons[i])
      sample_v_east.append(list_v_east[i]) 
      sample_v_north.append(list_v_north[i]) 
  return sample_lats, sample_lons, sample_v_east, sample_v_north

def test_euler_pole_from_quad():
  #test setup
  euler_pole = OC_NA_Pole #{"lat" : 45.0,  "long" : -90, "omega" : 1.23 }
  bearings  = [45.0, 135.0, 225.0, 315.0]
  sample_dist = 50000 # m

  sample_lons, sample_lats, sample_v_east, sample_v_north = create_simple_sample_quad(euler_pole, bearings, sample_dist)
  pole_result = epr.fit_euler_pole_linear(sample_lats, sample_lons, sample_v_east, sample_v_north, True)
  # ek.print_result ("test_euler_pole_from_quad", pole_result)

  assert pole_result['omega'] == pytest.approx(euler_pole['omega'])
  assert pole_result['lat'] == pytest.approx(euler_pole['lat'])
  assert pole_result['long'] == pytest.approx(euler_pole['long'])

def test_euler_pole_from_random_disk():
  #test setup
  euler_pole = {"lat" : 45.0,  "long" : -90, "omega" : 1.23 }
  sample_count = 400
  sample_dist = 50000 # m
  crop = 1.0 # no crop
  test_omega = 1.23

  sample_lats, sample_lons, sample_v_east, sample_v_north = \
    create_random_sample_ring(euler_pole, sample_count, sample_dist, test_omega, crop)
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
    create_random_sample_ring(
      euler_pole, 
      sample_count, 
      sample_dist, 
      test_omega, 
      crop,
      None,
      False)
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
    create_random_sample_ring(euler_pole, 
                              sample_count, 
                              sample_dist, 
                              test_omega, 
                              crop,
                              euler_n_pole,
                              False)

  pole_result = epr.fit_euler_pole_linear(sample_lats, sample_lons, sample_v_east, sample_v_north, True)
  #ek.print_result ("test_euler_pole_using_north_rotation", pole_result)

  assert pole_result['omega'] == pytest.approx(test_omega)
  assert pole_result['lat'] == pytest.approx(euler_n_pole['lat'])
  assert pole_result['long'] == pytest.approx(euler_n_pole['long'])

def test_GPS_pole_extraction():
  center_lat = 45.0
  center_long = -118
  max_distance = 350000 # m
  sample_lats, sample_lons, sample_v_east, sample_v_north = get_GPS_rotation_data(center_lat, center_long, max_distance)

  pole_result = epr.fit_euler_pole_linear(sample_lats, sample_lons, sample_v_east, sample_v_north)
  ek.print_result ("test_GPS_pole_extraction", pole_result)
  
  