import numpy as np
import geopandas as gpd    

def calculate_euler_pole_c(lats, lons, v_east_obs, v_north_obs):
    """
    Finds the exact best-fitting Euler pole using linear least squares.
    
    Args:
      lats (list/array): Latitudes of stations in decimal degrees
      lons (list/array): Longitudes of stations in decimal degrees
      v_east (list/array): East velocity components in mm/yr
      v_north (list/array): North velocity components in mm/yr
    """
    R = 6371.0 # Earth's radius in km
    
    # Convert input coordinates to radians
    num_stations = len(lats)
    
    # Initialize design matrix A and observation vector B
    # Each station provides 2 equations (East and North)
    A = np.zeros((2 * num_stations, 3))
    B = np.zeros(2 * num_stations)
    
    for i in range(num_stations):
        phi = lats[i]
        lam = lons[i]
        
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
    
    # Scale factor to output standard deg/Myr
    scale = (180 / np.pi)

    omega_deg_myr = omega_mag_rad * scale
    
    # 2. Latitude and Longitude of the Pole
    lat_pole = np.degrees(np.arcsin(wz / omega_mag_rad))
    lon_pole = np.degrees(np.arctan2(wy, wx))
    
    return {
        "pole_latitude_deg": lat_pole,
        "pole_longitude_deg": lon_pole,
        "rotation_rate_deg_per_myr": omega_deg_myr
    }

# def calculate_euler_pole_c(lats, lons, v_east, v_north):
#     """
#     Calculates the Euler Pole from lists of GPS points and velocities.
#     """
#     # Earth Radius in mm (to match velocity units of mm/yr)
#     R = 6371000.0 * 1000.0 

#     # Convert inputs to numpy arrays and convert degrees to radians
#     phi = np.radians(np.array(lats))
#     lam = np.radians(np.array(lons))

#     # Convert velocities from mm/yr to rad/yr relative to Earth's surface
#     ve = np.array(v_east) / R
#     vn = np.array(v_north) / R

#     num_stations = len(lats)

#     # Initialize the design matrix A and data vector b
#     A = np.zeros((2 * num_stations, 3))
#     b = np.zeros(2 * num_stations)

#     for i in range(num_stations):
#         # CORRECTED: East velocity equation index
#         A[2*i, 0] = -np.sin(phi[i]) * np.cos(lam[i])
#         A[2*i, 1] = -np.sin(phi[i]) * np.sin(lam[i])
#         A[2*i, 2] = np.cos(phi[i])   # This term cannot be zero!
#         b[2*i] = ve[i]

#         # CORRECTED: North velocity equation index
#         A[2*i+1, 0] = np.sin(lam[i])
#         A[2*i+1, 1] = -np.cos(lam[i])
#         A[2*i+1, 2] = 0.0            # Rotation around Z produces no local North motion
#         b[2*i+1] = vn[i]

#     # Solve the overdetermined linear system using Ordinary Least Squares
#     Omega, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)
#     print("Omega: ", Omega)

#     # Extract Cartesian rotation vector components (rad/yr)
#     wx, wy, wz = Omega

#     # Calculate the total magnitude of rotation (rad/yr)
#     omega_total = np.sqrt(wx**2 + wy**2 + wz**2)

#     # Convert total rotation rate to standard geodetic units (degrees / Million Years)
#     omega_deg_per_myr = omega_total * (180.0 / np.pi) * 1e6

#     # Calculate Euler Pole Location (Latitude and Longitude in degrees)
#     pole_lat = np.degrees(np.arcsin(wz / omega_total))
#     pole_lon = np.degrees(np.arctan2(wy, wx))

#     return {
#         "pole_latitude_deg": pole_lat,
#         "pole_longitude_deg": pole_lon,
#         "rotation_rate_deg_per_myr": omega_deg_per_myr
#     }

# print(calculate_euler_pole(lats, lons, v_east, v_north))

# def calculate_euler_pole(lats, lons, v_east, v_north):
#     """
#     Calculates the Euler Pole from lists of GPS points and velocities.
    
#     Parameters:
#     lats (list/array): Latitudes of stations in decimal degrees
#     lons (list/array): Longitudes of stations in decimal degrees
#     v_east (list/array): East velocity components in mm/yr
#     v_north (list/array): North velocity components in mm/yr
    
#     Returns:
#     dict: Dictionary containing pole latitude, longitude, and rotation rate.
#     """
#     # Earth Radius in mm (to match velocity units of mm/yr)
#     R = 6371000.0 * 1000.0 
    
#     # Convert inputs to numpy arrays and convert degrees to radians
#     phi = np.radians(np.array(lats))
#     lam = np.radians(np.array(lons))
    
#     # Convert velocities from mm/yr to rad/yr relative to Earth's surface
#     ve = np.array(v_east) / R
#     vn = np.array(v_north) / R
    
#     num_stations = len(lats)
    
#     # Initialize the design matrix A (2 equations per station) and data vector b
#     A = np.zeros((2 * num_stations, 3))
#     b = np.zeros(2 * num_stations)
    
#     for i in range(num_stations):
#         # East velocity equation index
#         A[2*i, 0] = -np.sin(lam[i])
#         A[2*i, 1] = np.cos(lam[i])
#         A[2*i, 2] = 0.0
#         b[2*i] = ve[i]
        
#         # North velocity equation index
#         A[2*i+1, 0] = -np.sin(phi[i]) * np.cos(lam[i])
#         A[2*i+1, 1] = -np.sin(phi[i]) * np.sin(lam[i])
#         A[2*i+1, 2] = np.cos(phi[i])
#         b[2*i+1] = vn[i]
        
#     # Solve the overdetermined linear system using Ordinary Least Squares
#     # A * Omega = b  ->  Omega = (A_T * A)^-1 * A_T * b
#     Omega, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)
#     print ("num_stations: ", num_stations)
#     print("A: ", A)
#     print("b: ", b)
#     print ("Omega ", Omega)

#     # Extract Cartesian rotation vector components (rad/yr)
#     wx, wy, wz = Omega
    
#     # Calculate the total magnitude of rotation (rad/yr)
#     omega_total = np.sqrt(wx**2 + wy**2 + wz**2)
    
#     # Convert total rotation rate to standard geodetic units (degrees / Million Years)
#     # (rad/yr) * (180/pi deg/rad) * (1,000,000 yr/Myr)
#     omega_deg_per_myr = omega_total * (180.0 / np.pi) * 1e6
    
#     # Calculate Euler Pole Location (Latitude and Longitude in degrees)
#     pole_lat = np.degrees(np.arcsin(wz / omega_total))
#     pole_lon = np.degrees(np.arctan2(wy, wx))
    
#     return {
#         "pole_latitude_deg": pole_lat,
#         "pole_longitude_deg": pole_lon,
#         "rotation_rate_deg_per_myr": omega_deg_per_myr
#     }


# # 1. Define physical constants and input parameters
# R = 6378137.0  # Earth's radius in meters (or standard mm, units just scale Omega)
# lat_val = 45.0
# longitudes_deg = [0.0, 120.0, -120.0]

# # Observed local velocities (East, North, Up)
# ve = 10.0
# vn = 0.0
# vu = 0.0

# A_list = []
# V_list = []

# # 2. Loop through each station to build the linear system
# for lon_deg in longitudes_deg:
#     lat = np.radians(lat_val)
#     lon = np.radians(lon_deg)
    
#     # Global Geocentric Cartesian Coordinates (X, Y, Z)
#     x = R * np.cos(lat) * np.cos(lon)
#     y = R * np.cos(lat) * np.sin(lon)
#     z = R * np.sin(lat)
    
#     # Transform local ENU velocity to global XYZ velocity (Crucial Step!)
#     # Formulated using standard geodetic rotation matrices:
#     vx = -np.sin(lon) * ve - np.sin(lat) * np.cos(lon) * vn + np.cos(lat) * np.cos(lon) * vu
#     vy =  np.cos(lon) * ve - np.sin(lat) * np.sin(lon) * vn + np.cos(lat) * np.sin(lon) * vu
#     vz =  0.0 * ve         + np.cos(lat) * vn               + np.sin(lat) * vu
    
#     # Skew-symmetric design matrix block derived from: v = Omega x r
#     A_i = np.array([
#         [ 0,  z, -y],
#         [-z,  0,  x],
#         [ y, -x,  0]
#     ])
    
#     A_list.append(A_i)
#     V_list.append([vx, vy, vz])

# # Combine all stations into global matrices
# A = np.vstack(A_list)
# V = np.array(V_list).flatten()

# # 3. Solve the Least Squares problem: Omega = (A^T * A)^-1 * A^T * V
# Omega, residuals, rank, s = np.linalg.lstsq(A, V, rcond=None)

# # 4. Convert global angular velocity vector back to Euler Pole parameters
# omega_mag = np.linalg.norm(Omega)
# pole_lat = np.degrees(np.arctan2(Omega[2], np.sqrt(Omega[0]**2 + Omega[1]**2)))
# pole_lon = np.degrees(np.arctan2(Omega[1], Omega[0]))

# # 5. Output results
# print("--- Regression Results ---")
# print(f"Global Omega Vector : [{Omega[0]:.4e}, {Omega[1]:.4e}, {Omega[2]:.4e}]")
# print(f"Angular Rotation Rate: {omega_mag:.6e} rad/time")
# print(f"Euler Pole Latitude  : {pole_lat:.1f}° N")
# print(f"Euler Pole Longitude : {pole_lon:.1f}° E")

def calculate_euler_pole(lats, lons, v_east, v_north):
#     """
#     Calculates the Euler Pole from lists of GPS points and velocities.
#     => Coded from writeup rather than asking Gemini to code it (yet) again...
    
#     Parameters:
#     lats (list/array): Latitudes of stations in decimal degrees
#     lons (list/array): Longitudes of stations in decimal degrees
#     v_east (list/array): East velocity components in mm/yr
#     v_north (list/array): North velocity components in mm/yr

    # 1. Define physical constants and input parameters
    R = 6.371E09  # Earth's radius in mm

    # Convert inputs to numpy arrays and convert degrees to radians
    phi = np.radians(np.array(lats))
    lam = np.radians(np.array(lons))

    # Convert velocities from mm/yr to rad/yr relative to Earth's surface
    ve = np.array(v_east) / R
    vn = np.array(v_north) / R

    num_stations = len(lats)

    # Initialize the design matrix A and data vector b
    A = np.zeros((3 * num_stations, 3))
    b = np.zeros(3 * num_stations)
        
    for i in range(num_stations):
        # define 3D position vectors
        x_i = R * np.cos(phi[i]) * np.cos(lam[i])
        y_i = R * np.cos(phi[i]) * np.sin(lam[i])
        z_i = R * np.sin(phi[i])

        #define 3D vecocity vectors
        v_x = -np.sin(lam[i]) * ve[i] - np.sin(phi[i]) * np.cos(lam[i]) * vn[i]
        v_y = np.cos(lam[i]) * ve[i] - np.sin(phi[i]) * np.sin(lam[i]) * vn[i]
        v_z = np.cos(phi[i]) * vn[i]
        v_mag = np.sqrt(v_x**2 + v_y**2 + v_z**2)

        # set up linear system
        A[3*i] =    [0.0, z_i, -y_i]
        A[3*i + 1] = [-z_i, 0, x_i]
        A[3*i + 2] = [y_i, -x_i, 0]

        b[3*i] = v_x
        b[3*i + 1] = v_y
        b[3*i + 2] = v_z

    # Solve the overdetermined linear system using Ordinary Least Squares
    Omega, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)
    print("Omega: ", Omega)

    # Extract Cartesian rotation vector components (rad/yr)
    wx, wy, wz = Omega

    # Calculate the total magnitude of rotation (rad/yr)
    omega_total = np.sqrt(wx**2 + wy**2 + wz**2)

    # Convert total rotation rate to standard geodetic units (degrees / Million Years)
    omega_deg_per_myr = omega_total * (180.0 / np.pi) * 1e6

    # Calculate Euler Pole Location (Latitude and Longitude in degrees)
    pole_lat = np.degrees(np.arcsin(wz / omega_total))
    pole_lon = np.degrees(np.arctan2(wy, wx))

    return {
        "pole_latitude_deg": pole_lat,
        "pole_longitude_deg": pole_lon,
        "rotation_rate_deg_per_myr": omega_deg_per_myr
    }
