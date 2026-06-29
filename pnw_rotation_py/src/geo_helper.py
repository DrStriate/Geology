import math
import numpy as np
from haversine import haversine, Unit
from pyproj import Geod

geod = Geod(ellps="WGS84")
R = 6371.0E3 # Earth radius in m

def get_northerly_easterly_from_lat_long_pts(lon1, lat1, lon2, lat2):
    # inv() expects (lon1, lat1, lon2, lat2)
    # forward_azimuth is the angle from point 1 to point 2 (degrees clockwise from North)
    forward_azimuth, back_azimuth, distance_meters = geod.inv(lon2, lat2, lon1, lat1)
    
    # Convert azimuth to radians
    azimuth_rad = np.radians(forward_azimuth)
    
    # Calculate components
    northerly = distance_meters * np.cos(azimuth_rad)
    easterly = distance_meters * np.sin(azimuth_rad)
    
    return northerly, easterly

# gets lat and long converted to coordinate distances from pole
def get_sample_pts(long_list, lat_list, pole):
  pe_list = []
  pn_list = []
  for i in range(len(long_list)):
    # convert sample points to meters
    p_n, p_e = get_northerly_easterly_from_lat_long_pts(long_list[i], lat_list[i], pole['long'], pole['lat'])
    pe_list.append(p_e)
    pn_list.append(p_n)
  return pe_list, pn_list 

def find_moments(long_list, lat_list, ve_list, vn_list, pole):
    sum_alpha = 0.0
    avg_alpha = 0.0
    count = len(ve_list) 
    pe_list, pn_list = get_sample_pts(long_list, lat_list, pole)
    for i in range(count):
        v = np.array([ve_list[i], vn_list[i]])
        p = np.array([pe_list[i], pn_list[i]])
        s = p + v
        norm_s = np.linalg.norm(s)
        norm_p = np.linalg.norm(p)
        dot_vp = np.dot(s, p)/(norm_p * norm_s)
        angle_vp = np.degrees(np.acos(dot_vp))
        # print(f"angle_vp: {angle_vp:.4f}")
        sum_alpha += angle_vp
    if count > 0:
        avg_alpha = sum_alpha / count
    return avg_alpha


def clamp(value, minimum, maximum):
    return max(minimum, min(value, maximum))

def latitudeFromDistN(dist): # dist in meters North
    lat = math.atan2(dist, R) * 180.0 / math.pi
    return lat

def longitudeFromDistE(latitude, dist): # meters East
    latitudeRadians = math.radians(latitude)
    radiusOfParallel = R * math.cos(latitudeRadians) # m
    longitudeDeltaRadians = dist / radiusOfParallel
    return math.degrees(longitudeDeltaRadians)

def LatLongForDeDn(latitude, longitude, de, dn): #de, dn in meters
    latOut = latitudeFromDistN(dn) + latitude
    longOut = longitudeFromDistE(latOut, de) + longitude
    return [latOut, longOut]

def DistanceFromLatLong(point1, point2): # (lat, lon)
    if point1[0] < -90 or point1[0] > 90 or point2[0] < -90 or point1[0] > 90:
        print("hello")
    return haversine(point1, point2, unit=Unit.METERS)

# OC_NA Eigen pole from Wells & Simpson 2001
def getPoleRotationV(lat, lon): # angles in degrees, motion per Ma
    OC_NAws = {"lat" : 45.54, "lon" : -119.6, "AngVel": 1.32 * 1E-6} 
    g = Geod(ellps='WGS84') 
    fwd_azimuth, back_azimuth, R = g.inv( OC_NAws["lon"],OC_NAws["lat"], lon, lat)
    p_hat = (np.sin(np.radians(fwd_azimuth)), np.cos(np.radians(fwd_azimuth))) # unit vector from lat lon
    d_hat = (p_hat[1], -p_hat[0]) # transpose is the rotation direction
    d = R * np.tan(np.radians(OC_NAws["AngVel"]))
    V = (d_hat[0] * d, d_hat[1] * d)
    return V

def get_pole_rotation_velocity(pole_lat, pole_lon, rotation_rate_deg_per_myr, lat, lon):
    g = Geod(ellps='WGS84') 
    fwd_azimuth, back_azimuth, R = g.inv( pole_lon, pole_lat, lon, lat)
    p_hat = (np.sin(np.radians(fwd_azimuth)), np.cos(np.radians(fwd_azimuth))) # unit vector from lat lon
    d_hat = (p_hat[1], -p_hat[0]) # transpose is the rotation direction
    d = R * np.tan(np.radians(rotation_rate_deg_per_myr)) * 10E-3 # m / Ma
    V = (d_hat[0] * d, d_hat[1] * d)  # mm / Yr
    return V





