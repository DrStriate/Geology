import numpy as np
from pyproj import Geod
import geo_helper as gh
from geo_helper import PLoc, EulerPole

geod = Geod(ellps="WGS84")
R = 6371.0E3 # Earth radius in m

# get cartesian value of PLoc lat/long
def getRVector(ploc):
  phi = np.radians(ploc.lat)
  lam = np.radians(ploc.long)
  return np.array([
      np.cos(phi) * np.cos(lam),  
      np.cos(phi) * np.sin(lam),
      np.sin(phi)
   ])

#get omega vector from pole (normalized - unscaled by pole omega)
def getWVector(pole):
  phi = np.radians(pole.lat)
  lam = np.radians(pole.long)
  return np.array([
      np.cos(phi) * np.cos(lam),  
      np.cos(phi) * np.sin(lam),
      np.sin(phi)
   ]) 

def getPoleRotationOfPoint(pole, ploc, ma):
    # Apply Rodrigues' rotation formula
    theta = np.radians(pole.omega) * ma
    v = getRVector(ploc) # v
    k = getWVector(pole) # k
    # 3. Apply Rodrigues' rotation formula
    k_cross_v = np.cross(k, v)
    k_dot_v = np.dot(k, v)
    v_new = v * np.cos(theta) + k_cross_v * np.sin(theta) + k * k_dot_v * (1 - np.cos(theta))
    
    # Convert the rotated Cartesian vector back to lat/long
    new_lon_rad = np.arctan2(v_new[1], v_new[0])
    new_lat_rad = np.arcsin(v_new[2])
    
    new_lon = np.degrees(new_lon_rad)
    new_lat = np.degrees(new_lat_rad)
    
    # Calculate the great-circle displacement on Earth's surface (R = 6371 km)
    angle_rad = np.arccos(np.clip(np.dot(v, v_new), -1.0, 1.0))
    displacement_km = gh.R * angle_rad
    
    return PLoc(new_lon, new_lat), displacement_km

  # Combo pole emulation 
def getCompoundRotationTranslationOfPoint(vPole, rPole, ploc, ma):
    # move the rot pole to the proper loc for ma
    rot_pole_ma_ploc = getPoleRotationOfPoint(vPole, rPole.ploc(), ma)[0]
    ma_rot_pole = gh.EulerPole(rot_pole_ma_ploc.long, rot_pole_ma_ploc.lat, rPole.omega)

    # Rotate by ma scaled rot pole omega 
    loc_2 = getPoleRotationOfPoint(ma_rot_pole, ploc, ma)[0]

    # Move rotated point up according to vPole ma
    loc_3 = getPoleRotationOfPoint(vPole, loc_2, -ma)[0]
    return loc_3

 # Big circle pole for given loc and velocity vector. pAvel is azimuth and speed (e.f. km/ma or mm/yr)
def getEulerPoleFromPlocAndPavel(ploc, pAvel):
    MetersPerDegree = 2 * np.pi * R / 360.0
    KmPerMaToDegreesPerMa = 1000.0 / MetersPerDegree
    degreesPerMa = pAvel.vel * KmPerMaToDegreesPerMa
    omega = degreesPerMa
    
    P = gh.getCartesianFromLanLong(ploc)
    V_e, V_n = gh.getVeVnFromAzvel(ploc, pAvel)
    V = V_e + V_n
    pe_hat = gh.normalize(np.cross(P, V)) # epipolar unit direction vector
    epiPoleLoc = gh.getPlocFromLocNormal(pe_hat)

    return EulerPole(epiPoleLoc.long, epiPoleLoc.lat, omega)

def testIfBigCircleCoplanarity(ploc1, ploc2, ploc3): # (r1 x r2) dot r3 == 0
    r1 = getRVector(ploc1)
    r2 = getRVector(ploc2)
    r3 = getRVector(ploc3)

    test = np.linalg.cross(r1, r2).dot(r3)
    return test

def getFwdAzimuth(ploc1, ploc2):
    # Initialize the WGS84 ellipsoid model
    geod = Geod(ellps='WGS84')
    
    # inv() expects longitude first, then latitude
    fwd_azimuth, back_azimuth, distance = geod.inv(ploc1.long, ploc1.lat, ploc2.lon, ploc2.lat)
    
    # Normalize azimuth to a 0-360 degree scale
    return fwd_azimuth % 360


# Code below probably needs to be refactored to use code/methods above which are more accurate
def create_sample (start_lon, start_lat, azimuth, distance):
    # Calculate the terminus point
    end_lon, end_lat, back_azimuth = geod.fwd(
        start_lon, 
        start_lat, 
        azimuth, 
        distance)
    return gh.PLoc(end_lon, end_lat)

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


def calculate_v_from_Euler_pole(Omega, p, omega): # p in {phi, lamb}, Omega in {phi, lamb, omega} radians
  P = R * get_hat_p(p)
  O = np.radians(omega) * get_hat_p(Omega)
  V = np.cross(P, O)
  v = project_V_to_v(V, p)
  return v

