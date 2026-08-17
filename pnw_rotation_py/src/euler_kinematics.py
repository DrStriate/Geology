import numpy as np
import geo_helper as gh
from geo_helper import PLoc, EulerPole, R

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

# Track yhs from ploc using NA plate motion and PNW rotation and plate motion info
def getPlocFromPoleData(naPAvel, pnwRotPole, pnwVPavel, ploc, ma):
  pnwVPole = getEulerPoleFromPlocAndPavel(pnwRotPole.ploc(), pnwVPavel)
  naPole = getEulerPoleFromPlocAndPavel(ploc, naPAvel)

  # move NA over yhs then move by both pnw poles to its ma location
  loc_2 = getPoleRotationOfPoint(naPole, ploc, -ma)[0]
  loc_3 = getCompoundRotationTranslationOfPoint(pnwVPole, pnwRotPole, loc_2, ma)
  return loc_3

 # Big circle pole for given loc and velocity vector. pAvel is azimuth and speed (e.f. km/ma or mm/yr)
def getEulerPoleFromPlocAndPavel(ploc, pAvel):
    MetersPerDegree = 2 * np.pi * R / 360.0
    KmPerMaToDegreesPerMa = 1.0 / MetersPerDegree
    degreesPerMa = pAvel.vel * KmPerMaToDegreesPerMa
    omega = degreesPerMa
    
    P = gh.getCartesianFromLatLong(ploc)
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

def project_V_to_v (V, p): #V is 3D cartesion velocity, phi and lamb in radians
  e_hat = np.array([-np.sin(p['lamb']), np.cos(p['lamb']), 0 ])
  n_hat = np.array([-np.sin(p['phi']) * np.cos(p['lamb']), -np.sin(p['phi']) * np.sin(p['lamb']), np.cos(p['phi'])])
  v_e = np.dot(V, e_hat)
  v_n = np.dot(V, n_hat)
  return np.array([v_e, v_n])

def get_hat_p(p): # returns a normal to the phi,lamb point
  return np.array([ 
    np.cos(p['phi']) * np.cos(p['lamb']),
    np.cos(p['phi']) * np.sin(p['lamb']),
    np.sin(p['phi'])
    ])

def get_hat(lat, long):
   return get_hat_p({'lamb': np.radians(long), 'phi': np.radians(lat)})

# omega not passed makes v based on euler_pole.omega, otherwise we use argument omega
def calculate_v_from_EulerPole(euler_pole, ploc, omega = None):
  Omega = {"omega": euler_pole.omega, "phi": np.radians(euler_pole.lat), "lamb": np.radians(euler_pole.long)}
  p = {"phi": np.radians(ploc.lat), "lamb": np.radians(ploc.long)}
  v = calculate_v_from_Euler_pole(Omega, p, euler_pole.omega if omega is None else omega)
  # BUG - not clear why the sign is flipped on V (see test_v_pole_from_sample_point)
  return -v

def calculate_v_from_Euler_pole(Omega, p, omega): # p in {phi, lamb}, Omega in {phi, lamb, omega} radians
  P = R * get_hat_p(p)
  O = np.radians(omega) * get_hat_p(Omega)
  V = np.cross(P, O)
  v = project_V_to_v(V, p)
  return v
