import math
import numpy as np
from haversine import haversine
from pyproj import Geod

EarthRadius = 6371 * 1000  # m

class GeoHelper:
    @staticmethod
    def clamp(value, minimum, maximum):
        return max(minimum, min(value, maximum))

    @staticmethod
    def latitudeFromDistN(dist): # dist in meters North
        lat = math.atan2(dist, EarthRadius) * 180.0 / math.pi
        return lat

    @staticmethod
    def longitudeFromDistE(latitude, dist): # meters East
        latitudeRadians = math.radians(latitude)
        radiusOfParallel = EarthRadius * math.cos(latitudeRadians) # m
        longitudeDeltaRadians = dist / radiusOfParallel
        return math.degrees(longitudeDeltaRadians)

    @staticmethod
    def LatLongForDeDn(latitude, longitude, de, dn): #de, dn in meters
        latOut = GeoHelper.latitudeFromDistN(dn) + latitude
        longOut = GeoHelper.longitudeFromDistE(latOut, de) + longitude
        return [latOut, longOut]

    def DistanceFromLatLong(point1, point2): # (lat, lon)
        return haversine(point1, point2)
    
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





