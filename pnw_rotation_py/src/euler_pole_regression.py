import numpy as np
from geo_helper import EulerPole, PLoc, getPAvel,R
import gauss_newton as gn
import test_utils as tu

# second best model (so far) but small errors in offset (0.002) but also works with GPS data 
def fit_euler_pole_linear(lats, lons, v_east_obs, v_north_obs, s_e = None, s_n = None):
    """
    Finds the exact best-fitting Euler pole using linear least squares.
    
    Args:
      lats (list/array): Latitudes of stations in decimal degrees
      lons (list/array): Longitudes of stations in decimal degrees
      v_east (list/array): East velocity components in mm/yr
      v_north (list/array): North velocity components in mm/yr

      align_pole tests for euler pole pointing to incoming n/s hemisphere (i.e Omega pole flip)
    """
    
    num_stations = len(lats)

    if num_stations < 3:
        return EulerPole(0, 0, 0)
    
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
    tu.test_regression_stats(omega_cartesian, A, B, residuals, False)

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
    
    return EulerPole(lon_pole, lat_pole, omega_deg_myr, is_clockwise=True)


# Best results so far in tests (with or without offset) but horrible in GPS - almost certainly due to crazy scaling
def fit_euler_pole_linear2(lats, lons, v_east_obs, v_north_obs, s_e=None, s_n=None):
    """
    Finds the best-fitting Euler pole and localized horizontal translation
    simultaneously, accounting for legacy m/Ma inputs safely.
    """
    num_stations = len(lats)
    
    A_joint = np.zeros((2 * num_stations, 5))
    B = np.zeros(2 * num_stations)
    
    sum_lats = 0

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
        P = R * np.array([np.cos(phi) * np.cos(lam), np.cos(phi) * np.sin(lam), np.sin(phi)])
        
        row_east_pole = np.cross(P, e_hat)
        row_north_pole = np.cross(P, n_hat)
        
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
    Offset_legacy = Offset_raw
    
    return EulerPole(lon_pole, lat_pole, omega_deg_myr, is_clockwise=True), Offset_legacy

# Early attempt at a unified regression (rot, offset). Gemini left it at T for offset, but it still seems incorrect in testing
def fit_euler_pole_linear3(lats, lons, v_east_obs, v_north_obs, s_e=None, s_n=None):
    """
    Finds the best-fitting Euler pole and global 3D translational offset
    simultaneously, preventing mathematical cross-talk across wide regions.
    """
    num_stations = len(lats)
    
    # 6 Columns: 3 for Rotation [Omega_x, Y, Z], 3 for 3D Translation [T_x, Y, Z]
    A_joint = np.zeros((2 * num_stations, 6))
    B = np.zeros(2 * num_stations)
    
    sum_lats = 0
    R_val = 6371000.0  # Spherical Earth radius in meters

    for i in range(num_stations):
        phi = np.radians(lats[i])
        lam = np.radians(lons[i])
        sum_lats += lats[i]
        
        # Apply standard data weights
        if s_n is not None and s_e is not None:
            sw_e = 1.0 / s_e[i]
            sw_n = 1.0 / s_n[i]
        else:
            sw_e = 1.0
            sw_n = 1.0

        # Construct local basis frames
        e_hat = np.array([-np.sin(lam), np.cos(lam), 0.0])
        n_hat = np.array([-np.sin(phi) * np.cos(lam), -np.sin(phi) * np.sin(lam), np.cos(phi)])
        
        # Build 3D Geocentric station position vector (meters)
        P = R_val * np.array([np.cos(phi) * np.cos(lam), np.cos(phi) * np.sin(lam), np.sin(phi)])
        
        # 1. Rigorous rigid block rotation row entries (O x P)
        row_east_pole = np.cross(P, e_hat)
        row_north_pole = np.cross(P, n_hat)
        
        idx_e = 2 * i
        idx_n = 2 * i + 1
        
        # Populate rotation matrix rows (0, 1, 2)
        A_joint[idx_e, 0:3] = row_east_pole * sw_e
        A_joint[idx_n, 0:3] = row_north_pole * sw_n
        
        # 2. Rigorous 3D Translation row entries (Projecting T onto local station frames)
        A_joint[idx_e, 3:6] = e_hat * sw_e
        A_joint[idx_n, 3:6] = n_hat * sw_n
        
        # 3. Populate observation target array
        B[idx_e] = v_east_obs[i] * sw_e
        B[idx_n] = v_north_obs[i] * sw_n
        
    north_hemisphere = (sum_lats > 0.0)
    
    # Run the joint linear regression pass
    Sol_joint, residuals, rank, s = np.linalg.lstsq(A_joint, B, rcond=None)
    
    # Isolate parameters cleanly
    Omega_c = Sol_joint[0:3]
    T_cartesian = Sol_joint[3:6]  # Uniform 3D ECEF translation vector (m/Ma)
    
    if (Omega_c[2] > 0) != north_hemisphere:
        Omega_c = -Omega_c
        
    Omega_mag = np.linalg.norm(Omega_c)
    omega_deg_myr = np.degrees(Omega_mag)  # Outputs cleanly as °/Myr
    
    lat_pole = np.degrees(np.arcsin(Omega_c[2] / Omega_mag))
    lon_pole = np.degrees(np.arctan2(Omega_c[1], Omega_c[0]))
    
    # Return the translation vector alongside the Euler Pole definition
    offset_e = np.dot(e_hat, T_cartesian)
    offset_n = np.dot(n_hat, T_cartesian)
    return EulerPole(lon_pole, lat_pole, omega_deg_myr, is_clockwise=True), np.array([offset_e, offset_n])

def extractEulerPoleUsingCombinedRegressions(lat_list, long_list, ve_list, vn_list, we_list = None, wn_list = None):

    # get data Euler pole from the raw data set
    raw_pole = fit_euler_pole_linear(lat_list, long_list, ve_list, vn_list, we_list, wn_list)
    # raw_pole.print("1: rotPole: ")

    # apply Gauss-Newton analysis to get any translation (non-rotation) components
    # offset = gn.solve_gauss_newton_translation_wtd(long_list, lat_list, ve_list, vn_list, we_list, wn_list, raw_pole)
    offset = gn.solve_gauss_newton_2D_transform_geo_wtd(long_list, lat_list, ve_list, vn_list, we_list, wn_list, raw_pole)
    # print(f"offset: {offset}")    

    # print(f"solve_gauss_newton_2D_transform_geo_wtd: {offset}")

    for iteration in range(3):
        # strip any translation element to get rot-only 
        rot_ve_list = np.array(ve_list) - offset[0]
        rot_vn_list = np.array(vn_list) - offset[1] 

        # get data Euler pole from the raw data set
        rot_pole = fit_euler_pole_linear(lat_list, long_list, rot_ve_list, rot_vn_list, we_list, wn_list)
        # rot_pole.print("2: rotPole: ")

        pass_offset = gn.solve_gauss_newton_2D_transform_geo_wtd(long_list, lat_list, rot_ve_list, rot_vn_list, we_list, wn_list, rot_pole)
        offset += pass_offset
        # print(f"offset: {offset}")
    # get Velocity pole PAVel info 
    pnwVPAVel = getPAvel(offset[0], offset[1])

    # offset is in meters per ma and we want a rate (km/ma or mm/yr) so we need to scale
    return rot_pole, pnwVPAVel

def getPnwGpsRotPoleAndVelocity(sample_center, sample_radius): # radius km 
  lats, longs, ves, vns, wes, wns=\
    tu.get_GPS_rotation_data(sample_center.long, sample_center.lat, sample_radius)   
  rot_pole, pnwVPAVel = extractEulerPoleUsingCombinedRegressions(lats, longs, ves, vns, wes, wns)
  return rot_pole, pnwVPAVel

def print_result(name, pole_result, point_count = 0):
    print(f"{name} count: {point_count}")
    print(f"Longitude: {pole_result.long:.5f}° E")   
    print(f"Latitude:  {pole_result.lat:.5f}° N")
    print(f"Rate:      {pole_result.omega:.5f}° / Myr")
