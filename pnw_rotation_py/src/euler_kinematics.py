import numpy as np
from pyproj import Geod

geod = Geod(ellps="WGS84")
R = 6371.0E3 # Earth radius in m
  
def create_sample (start_lon, start_lat, bearing, distance):
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

def get_hat_p(p): # returns a normal to the phi,lamb point
  return np.array([ 
    np.cos(p['phi']) * np.cos(p['lamb']),
    np.cos(p['phi']) * np.sin(p['lamb']),
    np.sin(p['phi'])
    ])

def get_hat(lat, long):
   return get_hat_p({'lamb': np.radians(long), 'phi': np.radians(lat)})


def calculate_v_from_Eigen_pole(Omega, p, omega): # p in {phi, lamb}, Omega in {phi, lamb, omega} radians
  P = R * get_hat_p(p)
  O = np.radians(omega) * get_hat_p(Omega)
  V = np.cross(P, O)
  v = project_V_to_v(V, p)
  return v

