from geo_helper import EulerPole, R, PLoc, geod # R is earth radius
import euler_pole_regression as epr
import gauss_newton as gn
import test_utils as tu
import pytest
import numpy as np
from pyproj import Geod

# Using new combined euler pole and offset regression
def test_quad_pole():
  geod = Geod(a=R, b=R) 
  euler_pole  = EulerPole(-99, 45.0, 1.32) # long, lat, omega
  azimuths  = [45.0, 135.0, 225.0, 315.0] # directions  to test lat/long
  sample_dist = 50000 # m
  realWorld = False
  # print("")

  #  Create samples, regress to pole
  sample_lons, sample_lats, sample_v_east, sample_v_north = create_simple_sample_quad1(euler_pole, azimuths, sample_dist, realWorld)
  pole_result = epr.fit_euler_pole_linear(sample_lats, sample_lons, sample_v_east, sample_v_north)
  #pole_result.print("pole_result")
   
  v_offset = gn.solve_gauss_newton_translation(sample_lats, sample_lons, sample_v_east, sample_v_north, pole_result)
  #print(f"v_offset1: {v_offset}")

  assert pole_result.omega == pytest.approx(euler_pole.omega, abs=1e-8)
  assert pole_result.lat == pytest.approx(euler_pole.lat, abs=1e-8)
  assert pole_result.long == pytest.approx(euler_pole.long, abs=1e-8)
  assert v_offset[0] == pytest.approx (0.0, abs=1e-8)
  assert v_offset[1] == pytest.approx (0.0, abs=1e-8)

# demonstrating new pole-offset results in same pole when v's moved by constant and good extraction
def test_offset_quad_pole():
  geod = Geod(a=R, b=R) 
  euler_pole  = EulerPole(-99, 45.0, 1.32) # long, lat, omega
  azimuths  = [45.0, 135.0, 225.0, 315.0] # directions  to test lat/long
  sample_dist = 50000 # m
  realWorld = False

  #  Create samples, regress to pole
  sample_lons, sample_lats, sample_v_east, sample_v_north = create_simple_sample_quad1(euler_pole, azimuths, sample_dist, realWorld)

  # relatively inaccurate global v change: using a pole would be more accurate
  sample_vs = [1.0, 2.0]
  sample_v_east += sample_vs[0]
  sample_v_north += sample_vs[1]

  pole_result = epr.fit_euler_pole_linear(sample_lats, sample_lons, sample_v_east, sample_v_north)
  pole_result.print("pole_result")

  v_offset = gn.solve_gauss_newton_translation(sample_lats, sample_lons, sample_v_east, sample_v_north, pole_result)
  # print(f"v_offset1: {v_offset}")

  assert pole_result.omega == pytest.approx(euler_pole.omega, abs=1e-5)
  assert pole_result.lat == pytest.approx(euler_pole.lat, abs=1e-3)
  assert pole_result.long == pytest.approx(euler_pole.long, abs=0.002)
  assert v_offset[0] == pytest.approx (sample_vs[0], rel=0.002) 
  assert v_offset[1] == pytest.approx (sample_vs[1], rel=0.002)

# repo of 'decomposed' regression we used to iterate to get offset. Note sb pole 
def test_euler_GPS_pole_extraction_legacy():
  geod = Geod(ellps='WGS84')
  center_lat = 45.0
  center_long = -118
  max_distance = 600000 # m
  lats, lons, v_easts, v_norths, s_e, s_n = tu.get_GPS_rotation_data(center_long, center_lat, max_distance)

  pole_result = epr.fit_euler_pole_linear(lats, lons, v_easts, v_norths, s_e, s_n)
  #epr.print_result ("test_GPS_pole_extraction", pole_result, len(lats))

  pole_result_sb = EulerPole(-115.41237, 43.63921,  0.5512292)
  assert pole_result.long == pytest.approx(pole_result_sb.long)
  assert pole_result.lat == pytest.approx(pole_result_sb.lat)
  assert pole_result.omega == pytest.approx(pole_result_sb.omega)

def test_euler_GPS_pole_extraction2():
  geod = Geod(ellps='WGS84')
  center_lat = 45.0
  center_long = -118
  max_distance = 600000 
  
  # 1. Elements loaded cleanly in unscaled mm/yr
  lats, lons, v_easts, v_norths, s_e, s_n = tu.get_GPS_rotation_data(center_long, center_lat, max_distance)

  # 2. Use your pristine baseline legacy pole calculator (which maps rad/yr perfectly if R is in mm)
  # Legacy expected mm/yr and R in meters? -> R_km scale matches your legacy wz extraction.
  pole_result = epr.fit_euler_pole_linear(lats, lons, v_easts, v_norths, s_e, s_n)
  #pole_result.print("Pristine Euler Pole")

  # 3. Pass parameters into the robust translation solver to extract the net sheet motion
  v_offset_mm_yr = gn.solve_gauss_newton_translation_wtd(lats, lons, v_easts, v_norths, s_e, s_n, pole_result)
  
  #print(f"Isolated Average Background Velocity (mm/yr East, North): {v_offset_mm_yr}")
  sb = np.array([ 837.27191367, 2943.58827792])
  assert v_offset_mm_yr == pytest.approx(sb)

def create_simple_sample_quad1(euler_pole, azimuths, dist, realWorld = False):
  longs = np.zeros(4)
  lats = np.zeros(4)
  v_easts = np.zeros(4)
  v_norths = np.zeros(4)

  Omega = {"omega": euler_pole.omega, "phi": np.radians(euler_pole.lat), "lamb": np.radians(euler_pole.long)}
  
  # Calculate angular radius in radians (distance / Earth radius)
  angular_radius = dist / R  

  for i in range(len(azimuths)):
    az_rad = np.radians(azimuths[i])
    
    # Spherically calculate points strictly relative to the pole coordinates
    # This ensures perfect small-circle symmetry
    sample_lat_rad = np.arcsin(np.sin(Omega['phi']) * np.cos(angular_radius) + 
                               np.cos(Omega['phi']) * np.sin(angular_radius) * np.cos(az_rad))
    
    sample_lon_rad = Omega['lamb'] + np.arctan2(np.sin(az_rad) * np.sin(angular_radius) * np.cos(Omega['phi']),
                                                np.cos(angular_radius) - np.sin(Omega['phi']) * np.sin(sample_lat_rad))
    
    sample_lat = np.degrees(sample_lat_rad)
    sample_lon = np.degrees(sample_lon_rad)
    
    p = {"phi": sample_lat_rad, "lamb": sample_lon_rad}
    v = calculate_v_from_Euler_pole(Omega, p, Omega['omega'], realWorld)

    # print(f"{i}: sample.long: {sample_lon:.3f}, sample.lat: {sample_lat:.3f}, v_e: {v['v_e']:.2f}  v_n: {v['v_n']:.2f}")


    longs[i] = sample_lon
    lats[i] = sample_lat
    v_easts [i] = v['v_e']
    v_norths[i] = v['v_n']
    
  return longs, lats, v_easts, v_norths


def create_sample (start_lon, start_lat, azimuth, distance):
    # Calculate the terminus point
    end_lon, end_lat, back_azimuth = geod.fwd(
        start_lon, 
        start_lat, 
        azimuth, 
        distance)
    return PLoc(end_lon, end_lat)


def project_V_to_v (V, p): #V is 3D cartesion velocity, phi and lamb in radians
  e_hat = np.array([-np.sin(p['lamb']), np.cos(p['lamb']), 0 ])
  n_hat = np.array([-np.sin(p['phi']) * np.cos(p['lamb']), -np.sin(p['phi']) * np.sin(p['lamb']), np.cos(p['phi'])])
  v_e = np.dot(V, e_hat)
  v_n = np.dot(V, n_hat)
  return {"v_e" : v_e, "v_n" : v_n}

def get_hat(lat, long):
   return get_hat_p({'lamb': np.radians(long), 'phi': np.radians(lat)})

# WGS84 Constants required for geodetic-to-cartesian conversions
WGS84_A = 6378137.0         # Semi-major axis (meters)
WGS84_E2 = 0.00669437999014  # First eccentricity squared

def get_hat_p(p): 
    # Returns a unit normal vector to the geodetic phi, lamb point
    return np.array([ 
        np.cos(p['phi']) * np.cos(p['lamb']),
        np.cos(p['phi']) * np.sin(p['lamb']),
        np.sin(p['phi'])
    ])


def calculate_v_from_Euler_pole(Omega, p, omega, realWorld):
    if realWorld: 
      # p: dict with {'phi': lat_rad, 'lamb': lon_rad}
      # Omega: dict with {'phi': pole_lat_rad, 'lamb': pole_lon_rad}
      # omega: scalar rotation rate (degrees/Myr or rad/yr depending on your scaling)
      
      # 1. Compute the local radius of curvature along the prime vertical (R_N)
      # This accounts for the ellipsoidal bulge at the station's exact latitude
      R_N = WGS84_A / np.sqrt(1.0 - WGS84_E2 * np.sin(p['phi'])**2)
      
      # 2. Build the 3D position vector using the local radius
      # For horizontal velocities, assume ellipsoidal height h = 0
      P = R_N * get_hat_p(p)
      
      # 3. Convert angular velocity to radians per unit time and build rotation vector
      omega_rad = np.radians(omega)
      O = omega_rad * get_hat_p(Omega)
      
      # 4. Standard rigid cross product 
      V = np.cross(O, P)  # Note: O x P yields the standard right-hand velocity vector
      
      # 5. Project 3D vector to local East and North ellipsoidal vectors
      v = project_V_to_v(V, p)
      return v
    else: # Idealized: usually functional test
      P = R * get_hat_p(p)
      O = np.radians(omega) * get_hat_p(Omega)
      V = np.cross(O, P)  # Standard kinematic rotation vector (Omega x P)
      v = project_V_to_v(V, p)
      return v

# Legacy 'decomposed' (separate and iterated pole and offset regressions)
def fit_euler_pole_linear_legacy(lats, lons, v_east_obs, v_north_obs, s_e = None, s_n = None):
    """
    Finds the exact best-fitting Euler pole using linear least squares.
    
    Args:
      lats (list/array): Latitudes of stations in decimal degrees
      lons (list/array): Longitudes of stations in decimal degrees
      v_east (list/array): East velocity components in mm/yr
      v_north (list/array): North velocity components in mm/yr

      align_pole tests for euler pole pointing to incoming n/s hemisphere (i.e Omega pole flip)
    """
    R = 6371.0E3 # Earth's radius in m
    
    num_stations = len(lats)
    
    # Initialize design matrix A and observation vector B
    A = np.zeros((2 * num_stations, 3))
    B = np.zeros(2 * num_stations)
    
    sum_lats = 0
    for i in range(num_stations):
        # Convert input coordinates to radians
        phi = np.radians(lats[i])
        lam = np.radians(lons[i])
        sum_lats += lats[i]
        
        # root weights for this station
        if s_n is not None and s_e is not None:
            sw_e = 1.0 / s_e[i]
            sw_n = 1.0 / s_n[i]
        else:
            sw_e = 1.0
            sw_n = 1.0
        
        # Weighted East velocity row equations (even rows: 2*i)
        A[2*i, 0] = -R * np.sin(phi) * np.cos(lam)  * sw_e
        A[2*i, 1] = -R * np.sin(phi) * np.sin(lam)  * sw_e
        A[2*i, 2] = R * np.cos(phi)                 * sw_e
        B[2*i]    = v_east_obs[i]                   * sw_e
        
        # Weighted North velocity row equations (odd rows: 2*i+1)
        A[2*i+1, 0] = R * np.sin(lam)               * sw_n
        A[2*i+1, 1] = -R * np.cos(lam)              * sw_n
        A[2*i+1, 2] = 0.0                           * sw_n
        B[2*i+1]    = v_north_obs[i]                * sw_n
        
    north_hemisphere = (sum_lats > 0.0)
    
    # Solves the weighted normal equations: A^T * W * A * omega = A^T * W * B
    omega_cartesian, residuals, rank, s = np.linalg.lstsq(A, B, rcond=None)
    
    wx, wy, wz = omega_cartesian

    if (wz > 0) != north_hemisphere: # if w and incoming data not in the same N/S hemisphere
        wx = -wx
        wy = -wy
        wz = -wz
    
    # Convert the Cartesian angular velocity vector back into Euler Pole parameters
    # 1. Total angular rotation magnitude (rad/yr converted back to deg/Myr)
    # Factor: (1e6 years * 180 degrees) / (pi radians * 1e9 mm to km conversion scale)
    # Since velocities are in mm/yr and R is in km, scaling matches naturally:
    omega_mag_rad = np.sqrt(wx**2 + wy**2 + wz**2) # rad per million years / 1000

    omega_deg_myr = np.degrees(omega_mag_rad) 

    # 2. Latitude and Longitude of the Pole
    lat_pole = np.degrees(np.arcsin(wz / omega_mag_rad))
    lon_pole = np.degrees(np.arctan2(wy, wx))
    
    return EulerPole(lon_pole, lat_pole, omega_deg_myr)

def fit_euler_pole_linear(lats, lons, v_east_obs, v_north_obs, s_e=None, s_n=None):
    """
    Finds the best-fitting Euler pole and localized horizontal translation
    simultaneously, accounting for legacy m/Ma inputs safely.
    """
    num_stations = len(lats)
    
    A_joint = np.zeros((2 * num_stations, 5))
    B = np.zeros(2 * num_stations)
    
    sum_lats = 0
    # Use Earth's radius in kilometers to keep the matrix columns well-conditioned
    R_km = 6371.0  

    for i in range(num_stations):
        phi = np.radians(lats[i])
        lam = np.radians(lons[i])
        sum_lats += lats[i]
        
        if s_n is not None and s_e is not None:
            sw_e = 1.0 / s_e[i]
            sw_n = 1.0 / s_n[i]
        else:
            sw_e = 1.0
            sw_n = 1.0

        e_hat = np.array([-np.sin(lam), np.cos(lam), 0.0])
        n_hat = np.array([-np.sin(phi) * np.cos(lam), -np.sin(phi) * np.sin(lam), np.cos(phi)])
        
        # 1. Build position vector in KILOMETERS 
        P_km = R_km * np.array([np.cos(phi) * np.cos(lam), np.cos(phi) * np.sin(lam), np.sin(phi)])
        
        row_east_pole = np.cross(P_km, e_hat)
        row_north_pole = np.cross(P_km, n_hat)
        
        idx_e = 2 * i
        idx_n = 2 * i + 1
        
        A_joint[idx_e, 0:3] = row_east_pole * sw_e
        A_joint[idx_n, 0:3] = row_north_pole * sw_n
        
        # 2. Convert translation columns to kilometer-scale scaling to match P_km matrix weight
        A_joint[idx_e, 3] = 1.0 * sw_e   
        A_joint[idx_e, 4] = 0.0
        A_joint[idx_n, 3] = 0.0
        A_joint[idx_n, 4] = 1.0 * sw_n   
        
        # 3. Scale input velocities down from m/Ma to mm/yr (which equals km/Ma)
        # This completely strips out the legacy x1000 multiplier during the inversion
        B[idx_e] = (v_east_obs[i] / 1000.0) * sw_e
        B[idx_n] = (v_north_obs[i] / 1000.0) * sw_n
        
    north_hemisphere = (sum_lats > 0.0)
    
    Sol_joint, residuals, rank, s = np.linalg.lstsq(A_joint, B, rcond=None)
    
    Omega_c = Sol_joint[0:3]
    Offset_raw = Sol_joint[3:5]  # Extracted cleanly in mm/yr (km/Ma)
    
    if (Omega_c[2] > 0) != north_hemisphere:
        Omega_c = -Omega_c
        
    Omega_mag = np.linalg.norm(Omega_c)
    omega_deg_myr = np.degrees(Omega_mag)  
    
    lat_pole = np.degrees(np.arcsin(Omega_c[2] / Omega_mag))
    lon_pole = np.degrees(np.arctan2(Omega_c[1], Omega_c[0]))
    
    # 4. Convert your output offset BACK to your legacy application's expected m/Ma scaling
    Offset_legacy = Offset_raw * 1000.0
    
    return EulerPole(lon_pole, lat_pole, omega_deg_myr), Offset_legacy



def solve_gauss_newton_translation_wtd(lats, lons, v_e_obs, v_n_obs, s_e, s_n, euler_pole):
    """
    Computes the true average horizontal/spatial translation offset vector 
    over a regional geodetic network using a weighted Gauss-Newton optimization.
    
    All inputs and outputs are processed in pure mm/yr (or km/Ma).
    """
    num_stations = len(lats)
    
    # Target parameter vector: 3D ECEF translational drift [Tx, Ty, Tz]
    # This prevents the flat-plane geometry distortion over wide networks
    T = np.zeros(3) 
    
    # Constants
    R_km = 6371.0 # Standard Earth radius in km
    Omega = {
        "omega": np.radians(euler_pole.omega), # Convert deg/Myr to rad/Myr (rad/Ma)
        "phi": np.radians(euler_pole.lat),
        "lamb": np.radians(euler_pole.long)
    }
    
    # Build the rotation vector O
    # unit normal to the pole * magnitude in rad/Ma
    O = Omega['omega'] * np.array([
        np.cos(Omega['phi']) * np.cos(Omega['lamb']),
        np.cos(Omega['phi']) * np.sin(Omega['lamb']),
        np.sin(Omega['phi'])
    ])

    # Gauss-Newton Iteration Loop (Linear project converges instantly in 1-2 passes)
    for iteration in range(3):
        J = np.zeros((2 * num_stations, 3))
        r = np.zeros(2 * num_stations)
        
        for i in range(num_stations):
            phi = np.radians(lats[i])
            lam = np.radians(lons[i])
            
            # Extract standard station weights
            sw_e = 1.0 / s_e[i] if s_e is not None else 1.0
            sw_n = 1.0 / s_n[i] if s_n is not None else 1.0
            
            # 1. Local geographic basis vectors
            e_hat = np.array([-np.sin(lam), np.cos(lam), 0.0])
            n_hat = np.array([-np.sin(phi) * np.cos(lam), -np.sin(phi) * np.sin(lam), np.cos(phi)])
            
            # 2. Position vector P in km (since velocity is mm/yr = km/Ma, R must be in km)
            P = R_km * np.array([np.cos(phi) * np.cos(lam), np.cos(phi) * np.sin(lam), np.sin(phi)])
            
            # 3. Rigid body rotation velocity vector (mm/yr)
            V_rot = np.cross(O, P)
            v_e_rot = np.dot(V_rot, e_hat)
            v_n_rot = np.dot(V_rot, n_hat)
            
            # 4. Projected current translation model parameters
            v_e_trans = np.dot(T, e_hat)
            v_n_trans = np.dot(T, n_hat)
            
            # 5. Combined predicted velocity components
            v_e_pred = v_e_rot + v_e_trans
            v_n_pred = v_n_rot + v_n_trans
            
            idx_e = 2 * i
            idx_n = 2 * i + 1
            
            # 6. Calculate Weighted Residuals (Observed - Predicted)
            r[idx_e] = (v_e_obs[i] - v_e_pred) * sw_e
            r[idx_n] = (v_n_obs[i] - v_n_pred) * sw_n
            
            # 7. Populate Jacobian Matrix (Partial derivatives scaled by weights)
            J[idx_e, 0:3] = e_hat * sw_e
            J[idx_n, 0:3] = n_hat * sw_n
            
        # Compute Gauss-Newton Parameter Update step: delta = (J^T * J)^-1 * J^T * r
        try:
            delta, _, _, _ = np.linalg.lstsq(J, r, rcond=None)
        except np.linalg.LinAlgError:
            break
            
        T += delta
        if np.linalg.norm(delta) < 1e-6:
            break
            
    # Convert the optimal 3D spatial drift vector back into a readable local horizontal 2D offset 
    # evaluated at the geometric center coordinates of your dataset for QGIS display
    center_phi = np.radians(np.mean(lats))
    center_lam = np.radians(np.mean(lons))
    
    c_e_hat = np.array([-np.sin(center_lam), np.cos(center_lam), 0.0])
    c_n_hat = np.array([-np.sin(center_phi) * np.cos(center_lam), -np.sin(center_phi) * np.sin(center_lam), np.cos(center_phi)])
    
    t_east_avg = np.dot(T, c_e_hat)
    t_north_avg = np.dot(T, c_n_hat)
    
    return np.array([t_east_avg, t_north_avg])
