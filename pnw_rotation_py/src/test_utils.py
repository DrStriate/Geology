import os
import numpy as np
import euler_kinematics as ek
import geopandas as gpd
import geo_helper as gh

OC_NA_Pole = {"lat" : 45.54,  "long" : -119.60, "omega" : 1.32 }

def get_data_file_path(name):
  current_dir = os.path.dirname(os.path.abspath(__file__))
  plugin_root = os.path.dirname(current_dir)
  data_folder_path = os.path.join(plugin_root, "data", name)
  data_folder_path = os.path.normpath(data_folder_path).replace("\\", "/")
  return data_folder_path

def get_test_data():
  return get_GPS_rotation_data(OC_NA_Pole['lat'], OC_NA_Pole['long'], 6e5)

MM_PER_YEAR_TO_M_PER_MA = 1000.0

def get_GPS_rotation_data (center_lat, center_long, max_distance):
  file_path = get_data_file_path("NSHM2023_GPS_velocity.zip")
  gdf = gpd.read_file(f"/vsizip/{file_path}")
  list_lats = gdf['geometry'].y.values
  list_lons = gdf['geometry'].x.values
  list_v_east = gdf['Ve'].values       
  list_v_north = gdf['Vn'].values   
  list_s_east = gdf['Se'].values       
  list_s_north = gdf['Sn'].values 
  "NSHM2023_GPS_velocity.zip"
  lats = []
  lons = []
  v_east = [] # mm/ yr
  v_north = [] # mm/ yr
  s_east = []
  s_north = []

  for i in range(len(list_lats)):
    dist = gh.DistanceFromLatLong((list_lats[i], list_lons[i]), (center_lat, center_long))
    if dist < max_distance:
      lats.append(list_lats[i])
      lons.append(list_lons[i])
      v_east.append(list_v_east[i] * MM_PER_YEAR_TO_M_PER_MA)
      v_north.append(list_v_north[i] * MM_PER_YEAR_TO_M_PER_MA)
      s_east.append(list_s_east[i] * MM_PER_YEAR_TO_M_PER_MA)
      s_north.append(list_s_north[i] * MM_PER_YEAR_TO_M_PER_MA)
  return lats, lons, v_east, v_north, s_east, s_north

def create_simple_sample_quad(euler_pole, bearings, dist):
  longs = []
  lats = []
  v_easts = []  # mm/ yr
  v_norths = [] # mm/ yr

  Omega = {"omega": euler_pole['omega'], "phi": np.radians(euler_pole['lat']), "lamb": np.radians(euler_pole['long'])}
  # print("")
  for i in range(len(bearings)):
    sample = ek.create_sample(euler_pole['long'], euler_pole['lat'], bearings[i], dist)
    p = {"phi": np.radians(sample['lat']), "lamb": np.radians(sample['lon'])}
    v = ek.calculate_v_from_Eigen_pole(Omega, p, Omega['omega']);
    longs.append(sample['lon'])
    lats.append(sample['lat'])
    v_easts.append(v['v_e'])
    v_norths.append(v['v_n'])
    #print(f"{i}: sample['lat']: {sample['lat']:.3f}, sample['lon']: {sample['lon']:.3f}, v_e: {v['v_e']:.2f}  v_n: {v['v_n']:.2f}")
  return longs, lats, v_easts, v_norths
  
def create_random_sample_ring(euler_pole, count, 
                              max_dist, 
                              test_omega, 
                              crop = 1.0, 
                              source_poll = None):
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

  max_long =  ek.create_sample(euler_pole['long'], euler_pole['lat'], 90.0, max_dist)['lon']
  min_long =  ek.create_sample(euler_pole['long'], euler_pole['lat'], 270.0, max_dist)['lon']
  crop_long = min_long + (max_long - min_long) * crop
  cropped_samples = 0;
  for i in range(len(rands)):
    sample = ek.create_sample(euler_pole['long'], euler_pole['lat'], 360.0 * rands[i][0], max_dist * rands[i][1])
    p = {"phi": np.radians(sample['lat']), "lamb": np.radians(sample['lon'])}
    v = ek.calculate_v_from_Eigen_pole(Omega_source, p, test_omega); 

    if sample['lon'] < crop_long:
      sample_e.append(sample['lon'])
      sample_n.append(sample['lat'])
      sample_v_east.append(v['v_e'])
      sample_v_north.append(v['v_n'])

    # print(f"{i}: sample['lat']: {sample['lat']:.3f}, sample['lon']: {sample['lon']:.3f}, v_e: {v['v_e']:.2f}  v_n: {v['v_n']:.2f}")
    cropped_samples += 1
  #print(f"samples = {cropped_samples} out of {count}")
  
  return sample_n, sample_e, sample_v_east, sample_v_north
