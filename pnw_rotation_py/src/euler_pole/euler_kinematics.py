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

def calculate_v_from_Eigen_pole(Omega, p, omega): # p in {phi, lamb}, Omega in {phi, lamb, omega} radians
  P = R * getHat(p)
  O = np.radians(omega) * getHat(Omega)
  V = np.cross(P, O)
  v = project_V_to_v(V, p)
  return v

def print_result(name, pole_result):
  print(name)
  print(f"Latitude:  {pole_result['lat']:.4f}° N")
  print(f"Longitude: {pole_result['long']:.4f}° E")
  print(f"Rate:      {pole_result['omega']:.4f}° / Myr")
