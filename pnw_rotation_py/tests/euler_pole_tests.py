# The following code stub will make the 'run Pyrthon file from a subfolder' work,
# in VS Code but it is not recommended to use this in production code.
from config_root import setRoot
setRoot()

import src.euler_pole.euler_pole_regression as epr
import numpy as np
import geopandas as gpd    

# # ==========================================
# EXAMPLE USAGE WITH SYNTHETIC DATA
# ==========================================
# # Suppose we have 3 GPS stations tracking a rigid block rotation
# sample_lats  = [34.0, 35.5, 36.2]
# sample_lons  = [-118.0, -119.2, -117.5]
# sample_v_east = [15.2, 14.8, 16.1]     # mm/yr
# sample_v_north = [5.1, 4.8, 5.5]       # mm/yr

# Suppose we have 4 GPS stations tracking a rigid block rotation around polaris at 45 deg N
sample_lons  = [-45.0, 45.0, 45.0, -45.0]
sample_lats  = [45.0, 45.0, -45.0, -45.0]
sample_v_east = [10.0, 10.0, -10.0, -10.0]     # mm/yr
sample_v_north = [10.0, -10.0, -10.0, 10.0]    # mm/yr

# Suppose we have 4 GPS stations tracking a rigid block rotation around polaris at equator (0 deg lat)
# 10mm at equator => 8.99 x 10E-8 degrees / year
# sample_lons  = [-90.0, -180.0, 90.0, 0.0]
# sample_lats  = [0.0, 0.0, 0.0, 0.0]
# sample_v_east = [10.0, 10.0, 10.0, 10.0]     # mm/yr
# sample_v_north = [0.0, 0.0, 0.0, 0.0]    # mm/yr

# Suppose we have 3 GPS stations tracking a rigid block rotation
# sample_lats  = [45.0, 45.0, 45.0]
# sample_lons  = [0, 120.0, -120.0]
# sample_v_east = [10, 10, 10]     # mm/yr
# sample_v_north = [0.0, 0.0, 0.0] # mm/yr

# ==========================================
# EXAMPLE USAGE WITH GPS DATA
# ==========================================
# Load the shapefile data directly from the zip
gdf = gpd.read_file("zip://data/NSHM2023_GPS_velocity.zip")
# sample_lats = gdf['geometry'].y.values
# sample_lons = gdf['geometry'].x.values
# sample_v_east = gdf['Ve'].values       
# sample_v_north = gdf['Vn'].values    

# Calculate the pole
pole_result = epr.calculate_euler_pole_c(sample_lats, sample_lons, sample_v_east, sample_v_north)

# Print results cleanly
print("--- Calculated Euler Pole ---")
print(f"Latitude:  {pole_result['pole_latitude_deg']:.4f}° N")
print(f"Longitude: {pole_result['pole_longitude_deg']:.4f}° E")
print(f"Rate:      {pole_result['rotation_rate_deg_per_myr']:.4f}° / Myr")