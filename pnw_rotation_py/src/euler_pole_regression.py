import numpy as np
from geo_helper import EulerPole, getPAvel, PLoc
import test_utils as tu
from gauss_newton import solve_gauss_newton_2D_transform_geo_wtd

def fit_euler_pole_linear(lats, lons, v_east_obs, v_north_obs, align_pole = True):
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
    # Each station provides 2 equations (East and North)
    A = np.zeros((2 * num_stations, 3))
    B = np.zeros(2 * num_stations)
    sum_lats = 0
    
    for i in range(num_stations):
        # Convert input coordinates to radians
        phi = np.radians(lats[i])
        lam = np.radians(lons[i])
        sum_lats += lats[i]
        
        # East velocity row equations
        A[2*i, 0] = -R * np.sin(phi) * np.cos(lam)  # R * n-hat[0] 
        A[2*i, 1] = -R * np.sin(phi) * np.sin(lam)  # R * n_hat[1]
        A[2*i, 2] = R * np.cos(phi)                 # R * n_hat[2]
        B[2*i] = v_east_obs[i]
        
        # North velocity row equations
        A[2*i+1, 0] = R * np.sin(lam)               # -R * e-hat[0]
        A[2*i+1, 1] = -R * np.cos(lam)              # -R * e-hat[1]
        A[2*i+1, 2] = 0.0                           # -R * e-hat[2]
        B[2*i+1] = v_north_obs[i]
    
    north_hemisphere = (sum_lats > 0.0)
        
    # Solve the linear system A * omega = B using standard least squares
    # This solves the normal equation: omega = (A^T * A)^(-1) * A^T * B
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

def fit_euler_pole_linear_wtd(lats, lons, v_east_obs, v_north_obs, s_e, s_n, align_pole = True):
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
        
        # Current root weights for this station
        sw_e = 1.0 / s_e[i]
        sw_n = 1.0 / s_n[i]
        
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

def extractEulerPoleUsingCombinedRegressions(lat_list, long_list, ve_list, vn_list, we_list, wn_list):
    # get data Euler pole from the raw data set
    raw_pole = fit_euler_pole_linear_wtd(lat_list, long_list, ve_list, vn_list, we_list, wn_list)

    # apply Gauss-Newton analysis to get any translation (non-rotation) components
    gn_out = solve_gauss_newton_2D_transform_geo_wtd(long_list, lat_list, ve_list, vn_list, we_list, wn_list, raw_pole.ploc())

    # strip any translation element to get rot-only 
    rot_ve_list = np.array(ve_list) - gn_out['t_x']
    rot_vn_list = np.array(vn_list) - gn_out['t_y']

    # get data Euler pole from the raw data set
    rot_pole = fit_euler_pole_linear_wtd(lat_list, long_list, rot_ve_list, rot_vn_list, we_list, wn_list)

    # get Velocity pole PAVel info 
    pnwVPAVel = getPAvel(gn_out['t_x'] * 1e-3, gn_out['t_y'] * 1e-3)

    # offset is in meters per ma and we want a rate (km/ma or mm/yr) so we need to scale
    return rot_pole, pnwVPAVel

def print_result(name, pole_result, point_count = 0):
    print(f"{name} count: {point_count}")
    print(f"Longitude: {pole_result.long:.5f}° E")   
    print(f"Latitude:  {pole_result.lat:.5f}° N")
    print(f"Rate:      {pole_result.omega:.5f}° / Myr")

def getPnwGpsRotPoleAndVelocity(self, raw_data_center, distance):
    # get rot data
    sample_radius = 600.0 # km
    sample_center = PLoc(-119.0, 45.0)
    lats, longs, ves, vns, wes, wns=\
      tu.get_GPS_rotation_data(sample_center.long, sample_center.lat, sample_radius * 1000)   
    rot_pole, pnwVPAVel = extractEulerPoleUsingCombinedRegressions(lats, longs, ves, vns, wes, wns)
    return rot_pole, pnwVPAVel