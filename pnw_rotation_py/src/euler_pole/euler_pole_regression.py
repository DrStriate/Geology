import numpy as np

def fit_euler_pole_linear(lats, lons, v_east_obs, v_north_obs):
    """
    Finds the exact best-fitting Euler pole using linear least squares.
    
    Args:
      lats (list/array): Latitudes of stations in decimal degrees
      lons (list/array): Longitudes of stations in decimal degrees
      v_east (list/array): East velocity components in mm/yr
      v_north (list/array): North velocity components in mm/yr
    """
    R = 6371.0E03 # Earth's radius in m
    
    # Convert input coordinates to radians
    num_stations = len(lats)
    
    # Initialize design matrix A and observation vector B
    # Each station provides 2 equations (East and North)
    A = np.zeros((2 * num_stations, 3))
    B = np.zeros(2 * num_stations)
    
    for i in range(num_stations):
        phi = np.radians(lats[i])
        lam = np.radians(lons[i])
        
        # East velocity row equations
        A[2*i, 0] = -R * np.sin(phi) * np.cos(lam)
        A[2*i, 1] = -R * np.sin(phi) * np.sin(lam)
        A[2*i, 2] = R * np.cos(phi)
        B[2*i] = v_east_obs[i]
        
        # North velocity row equations
        A[2*i+1, 0] = R * np.sin(lam)
        A[2*i+1, 1] = -R * np.cos(lam)
        A[2*i+1, 2] = 0.0
        B[2*i+1] = v_north_obs[i]
        
    # Solve the linear system A * omega = B using standard least squares
    # This solves the normal equation: omega = (A^T * A)^(-1) * A^T * B
    omega_cartesian, residuals, rank, s = np.linalg.lstsq(A, B, rcond=None)
    
    wx, wy, wz = omega_cartesian
    
    # Convert the Cartesian angular velocity vector back into Euler Pole parameters
    # 1. Total angular rotation magnitude (rad/yr converted back to deg/Myr)
    # Factor: (1e6 years * 180 degrees) / (pi radians * 1e9 mm to km conversion scale)
    # Since velocities are in mm/yr and R is in km, scaling matches naturally:
    omega_mag_rad = np.sqrt(wx**2 + wy**2 + wz**2) # rad per million years / 1000

    omega_deg_myr = np.degrees(omega_mag_rad) 

    # 2. Latitude and Longitude of the Pole
    lat_pole = np.degrees(np.arcsin(wz / omega_mag_rad))
    lon_pole = np.degrees(np.arctan2(wy, wx))
    
    return {
        "latitude": lat_pole,
        "longitude": lon_pole,
        "rate (deg/Ma)": omega_deg_myr
    }

# --- Examplwe Usage ---
sample_lons  = np.array([-90.0, -180.0, 90.0, 0.0])
# sample_lats  = np.array([0.0, 0.0, 0.0, 0.0])
# obs_ve = np.array([10.0, 10.0, 10.0, 10.0])
# obs_vn = np.array([0.0, 0.0, 0.0, 0.0])

sample_lats = np.array([48.23, 45.54, 42.844, 45.54])
sample_Lons = np.array([-119.6, -115.748, -119.6, -123.452])
obs_ve = np.array([62.256, 1.5, -66.227, 1.5])
obs_vn = np.array([0.0, -62.470, 0.0, 62.47]) 

found_pole = fit_euler_pole_linear(sample_lats, sample_lons, obs_ve, obs_vn)

print(f"Latitude: {found_pole['latitude']:.2f}°N")
print(f"Pole Longitude: {found_pole['longitude']:.2f}°E")
print(f"Angular Velocity: {found_pole['rate (deg/Ma)']:.4f} °/Myr")
