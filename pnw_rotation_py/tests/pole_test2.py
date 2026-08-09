from geo_helper import EulerPole, R, PLoc # R is earth radius
import euler_pole_regression as epr
import test_utils as tu
import pytest
import numpy as np
from pyproj import Geod

def test_euler_GPS_pole_extractionX():
  center_lat = 45.0
  center_long = -118
  max_distance = 600000 
  
  # 1. Elements loaded cleanly in unscaled mm/yr
  lats, lons, v_easts, v_norths, s_e, s_n = tu.get_GPS_rotation_data(center_long, center_lat, max_distance)

  # 2. Use your pristine baseline legacy pole calculator (which maps rad/yr perfectly if R is in mm)
  # Legacy expected mm/yr and R in meters? -> R_km scale matches your legacy wz extraction.
  pole_result = fit_euler_pole_linear_legacy(lats, lons, v_easts, v_norths, s_e, s_n)
  pole_result.print("Pristine Euler Pole")

  # 3. Pass parameters into the robust translation solver to extract the net sheet motion
  v_offset_mm_yr = solve_gauss_newton_translation_wtd(lats, lons, v_easts, v_norths, s_e, s_n, pole_result)
  
  print(f"Isolated Average Background Velocity (mm/yr East, North): {v_offset_mm_yr}")

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
