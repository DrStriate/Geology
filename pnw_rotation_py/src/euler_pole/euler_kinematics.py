import numpy as np
from pyproj import Geod
# import euler_pole_regression as epr

R = 6371.0E3 # Earth radius in m
  
def create_sample (start_lon, start_lat, bearing, distance):
    geod = Geod(ellps="WGS84")
    # Calculate the terminus point
    end_lon, end_lat, back_azimuth = geod.fwd(
        start_lon, 
        start_lat, 
        bearing, 
        distance)
    return {"lon": end_lon, "lat": end_lat}

def project_V_to_v (V, p): #V is 3D cartesion velocity, phi and lamb in radians
  e_hat = np.array([-np.sin(p['lamb']), np.cos(p['lamb']), 0 ])
  n_hat = np.array([-np.sin(p['phi']) * np.cos(p['lamb']), -np.sin(p['phi']) * np.sin(p['lamb']), np.cos(p['phi'])])
  v_e = np.dot(V, e_hat)
  v_n = np.dot(V, n_hat)
  return {"v_e" : v_e, "v_n" : v_n}

def getHat(p): # returns a normal to the lat,long point
  return np.array([ 
    np.cos(p['phi']) * np.cos(p['lamb']),
    np.cos(p['phi']) * np.sin(p['lamb']),
    np.sin(p['phi'])
    ])

def calculate_v_from_Eigen_pole(Omega, p): # p in {phi, lamb}, Omega in {phi, lamb, omega} radians
  P = R * getHat(p)
  O = Omega['omega'] * getHat(Omega)
  V = np.cross(P, O)
  v = project_V_to_v(V, p)
  return v

def print_result(name, pole_result):
  print(name)
  print(f"Latitude:  {pole_result['latitude']:.4f}° N")
  print(f"Longitude: {pole_result['longitude']:.4f}° E")
  print(f"Rate:      {pole_result['rate (deg/Ma)']:.4f}° / Myr")

# # Test: we have 4 GPS stations tracking a rigid block rotation around the (actual Euler pole)

# #test setup
# euler_point = {"lat" : 45.0,  "long" : -90}
# omega = 1.23 # degrees / Ma (or mm/yr)
# sample_bearings  = [45.0, 135.0, 225.0, 315.0]
# sample_dist = 50000 # m

# sample_lats = []
# sample_lons = []
# sample_v_east = []
# sample_v_north = [] # mm/ yr

# Omega = {"omega": np.radians(omega), "phi": np.radians(euler_point['lat']), "lamb": np.radians(euler_point['long'])}
# for i in range(len(sample_bearings)):
#   sample = create_sample(euler_point['long'], euler_point['lat'], sample_bearings[i], sample_dist)
#   print(f"{i}: sample: {sample}")
#   sample_lats.append(sample['lat'])
#   sample_lons.append(sample['lon'])
#   p = {"phi": np.radians(sample_lats[i]), "lamb": np.radians(sample_lons[i])}
#   v = calculate_v_from_Eigen_pole(Omega, p); 
#   sample_v_east.append(v['v_e'])
#   sample_v_north.append(v['v_n'])
#   print(f"{i}: v: {v}")

# pole_result5 = epr.fit_euler_pole_linear(sample_lats, sample_lons, sample_v_east, sample_v_north)
# print_result ("fit_euler_pole_linear", pole_result5)
