import geo_helper as gh
import geopandas as gpd   
import pytest

NA_Speed = 35e-3    # m / yr (Current)
NA_Bearing = 240.0  # degrees azimuth
YHS_lat = 44.43     # Yellowstone hotspot caldera
YHS_long = -110.67
JdF_Lat = 48.25     # Center of Strait of Juan de Fuca
JdF_Long = -124.0

def get_GPS_rotation_data (center_lat, center_long, max_distance):
  gdf = gpd.read_file("zip://data/NSHM2023_GPS_velocity.zip")
  list_lats = gdf['geometry'].y.values
  list_lons = gdf['geometry'].x.values
  list_v_east = gdf['Ve'].values       
  list_v_north = gdf['Vn'].values    
  
  sample_lats = []
  sample_lons = []
  sample_v_east = []
  sample_v_north = [] # mm/ yr
  for i in range(len(list_lats)):
    dist = gh.DistanceFromLatLong((list_lats[i], list_lats[i]), (center_lat, center_long))
    if dist < max_distance:
      sample_lats.append(list_lats[i])
      sample_lons.append(list_lons[i])
      sample_v_east.append(list_v_east[i]) 
      sample_v_north.append(list_v_north[i]) 
  return sample_lats, sample_lons, sample_v_east, sample_v_north

def test_latitudeFromDistance():
    distN = 1000 # 1 km
    latSb = 0.00898 # degrees
    lat = gh.latitudeFromDistN(distN)
    assert lat == pytest.approx(latSb, abs=0.001)

def test_lonitudeFromDistance():
    distN = 1000 # 1 km
    lat = 45
    lonSb = 0.0127 # degrees
    lon = gh.longitudeFromDistE(lat, distN)
    assert lon == pytest.approx(lonSb, abs=0.001)

def test_pole_rotation():
  yhs = {"lat" : 44.43, "lon" :-110.67}
  V = gh.getPoleRotationV(yhs["lat"], yhs["lon"])
  assert V[0] == -0.0019364865124427936
  assert V[1] == -0.016349206364427868
  return
