import numpy as np
from haversine import haversine, Unit
from pyproj import Geod
from dataclasses import dataclass


R = 6371.0 # Earth radius in km

geod = Geod(ellps="WGS84")

def kmPerDegree():
    return R * np.pi / 180.0

def getMagnitude(d_e, d_n):    # only use for very small distances
    return np.sqrt(d_e * d_e + d_n * d_n)

def getAzimuth (v_e, v_n): 
    azimuth = np.degrees(np.arctan2(v_e, v_n))
    if azimuth < 0.0 :
        azimuth += 360.0 # arctan2 returns -180 ... 180 wheras azimuth needs to be 0 ... 360
    return azimuth

def getPAvel (v_e, v_n):
    return PAvel(getAzimuth(v_e, v_n), getMagnitude(v_e, v_n))

@dataclass
class PLoc:
    long: float
    lat: float
    def print(self, label = ""): 
        print (f"{label} long: {self.long:0.3f}, lat:  {self.lat:0.3f}")
    
    def __iadd__(self, other):  # Overriding the += operator
        if isinstance(other, PLoc):
            self.long += other.long
            self.lat += other.lat
            return self 
        return NotImplemented

@dataclass 
# PVel should only be used for (linear) Velocity. For V the standard is mm/Yr (= km/ma) 
# distance measurements in plate kinematics should be in degrees, not kilometers.
class PVel:
    east: float
    north: float
    def azimuth (self):
        return getAzimuth(self.east, self.north)
    def magnitude (self): # only valid for very small magnitudes
        return getMagnitude(self.east, self.north)
    def print(self, label = ""): 
        print (f"{label} east: {self.east:0.3f}, north:  {self.north:0.3f}")

@dataclass
# PAvel azimuth is in degrees (cw from N) and vel is in  mm/yr or km/Ma (equiv)
# PAvel is an Azimuth & Speed proxy for PVel
# as with PVel, distance measurements in plate kinematics should be in degrees, not kilometers.
class PAvel:
    azimuth: float
    vel: float
    def print(self, label = ""): 
        print (f"{label} azimuth: {self.azimuth:0.3f}, vel:  {self.vel:0.3f}")
    #def getAzimuth (ploc1, ploc2): #useful to track forward azimuth on 
       
@dataclass
class EulerPole:
    long: float
    lat: float
    omega: float
    def ploc(self):
        return PLoc(self.long, self.lat)
    # def setPloc(self, ploc):
    #     self.lat = ploc.lat
    #     self.long = ploc.long
    def print(self, label = ""):
        print(f"{label} long: {self.long:0.3f}, lat: {self.lat:0.3f}, omega: {self.omega:.6f}")
    # def pointRotateForMa(self, ploc, ma):
    #     return ek.getPoleRotationOfPoint(self, ploc, ma)[0]
###

# Get point from PAvel new point pavel * ma distant
def getPointFromPavel(start_point, pAVel, ma):
    """
    Calculates the destination latitude/longitude using a spherical Earth model.
    """
    # Convert degrees to radians
    lat1 = np.radians(start_point.lat)
    lon1 = np.radians(start_point.long)
    azimuth = np.radians(pAVel.azimuth)
    
    # Angular distance covered
    angular_vel = ( pAVel.vel * ma ) / R
    
    # Calculate destination latitude
    lat2 = np.arcsin(np.sin(lat1) * np.cos(angular_vel) +
                     np.cos(lat1) * np.sin(angular_vel) * np.cos(azimuth))
    
    # Calculate destination longitude
    lon2 = lon1 + np.arctan2(np.sin(azimuth) * np.sin(angular_vel) * np.cos(lat1),
                             np.cos(angular_vel) - np.sin(lat1) * np.sin(lat2))
    
    # Convert back from radians to degrees
    destination_lat = np.degrees(lat2)
    destination_lon = np.degrees(lon2)
    
    # Normalize longitude to be between -180 and +180
    destination_lon = (destination_lon + 540) % 360 - 180
    
    return PLoc(destination_lon, destination_lat)

# This is for angular easterly and northerly rotation angles, not dists
def getNortherlyEasterlyFromLatLongPoints(lon1, lat1, lon2, lat2):
    # forward_azimuth is the angle from point 1 to point 2 (degrees clockwise from North)
    forward_azimuth, back_azimuth, distance_meters = geod.inv(lon1, lat1, lon2, lat2)

    # Convert azimuth to radians
    azimuth_rad = np.radians(forward_azimuth)
    
    # Calculate components - actually calculating 
    northerly = distance_meters * np.cos(azimuth_rad)
    easterly = distance_meters * np.sin(azimuth_rad)
    return northerly, easterly

def getFwdAzimuthFromLocations (point1, point2):
   # forward_azimuth is the angle from point 1 to point 2 (degrees clockwise from North)
    forward_azimuth, back_azimuth, distance_meters = geod.inv(point1.long, point1.lat, point2.long, point2.lat)
    return forward_azimuth

# Code below probably needs to be refactored to use code/methods above which are more accurate
def create_sample (start_lon, start_lat, azimuth, distance):
    # Calculate the terminus point
    end_lon, end_lat, back_azimuth = geod.fwd(
        start_lon, 
        start_lat, 
        azimuth, 
        distance)
    return PLoc(end_lon, end_lat)

# using geo_helper for geod to keep models consistent
# def getFwdAzimuth(ploc1, ploc2):
#     # Initialize the WGS84 ellipsoid model
#     geod = Geod(ellps='WGS84')
    
#     # inv() expects longitude first, then latitude
#     fwd_azimuth, back_azimuth, distance = geod.inv(ploc1.long, ploc1.lat, ploc2.lon, ploc2.lat)
    
#     # Normalize azimuth to a 0-360 degree scale
#     return fwd_azimuth % 360

# gets lat and long converted to coordinate distances from pole. This is an approximation 
def getSamplePoints(long_list, lat_list, center_ploc):
  pe_list = []
  pn_list = []
  for i in range(len(long_list)):
    # convert sample points to meters
    p_n, p_e = getNortherlyEasterlyFromLatLongPoints(center_ploc.long, center_ploc.lat, long_list[i], lat_list[i])
    pe_list.append(p_e)
    pn_list.append(p_n)
  return pe_list, pn_list 

def clamp(value, minimum, maximum):
    return max(minimum, min(value, maximum))

# Be wary of use of these distance metrics. 
def latitudeFromDistN(dist): # dist in meters North
    lat = np.arctan2(dist, R) * 180.0 / np.pi
    return lat

def longitudeFromDistE(latitude, dist): # meters East
    latitudeRadians = np.radians(latitude)
    radiusOfParallel = R * np.cos(latitudeRadians) # m
    longitudeDeltaRadians = dist / radiusOfParallel
    return np.degrees(longitudeDeltaRadians)

# Great circle distance
def getDistanceBetweenPoints(point1, point2): #both PLocs
    if point1.lat < -90 or point1.lat > 90 or point2.lat < -90 or point1.lat > 90:
        print("getDistanceBetweenPoints: bounding error")
    return haversine((point1.lat, point1.long), (point2.lat, point2.long), unit=Unit.METERS)

### Epipolar calculations
def locToRadians(pLoc):
    lam = np.radians(pLoc.long)
    phi = np.radians(pLoc.lat) 
    return lam, phi

def normalize(vect):
    mag = np.linalg.norm(vect)
    if mag > 0:
        return vect / mag
    return vect

def getCartesianFromLatLong (pLoc):
    lam, phi = locToRadians(pLoc)
    P = np.array([0, 0, 0])
    P[0] = R * np.cos(lam) * np.cos(phi)
    P[1] = R * np.sin(lam) * np.cos(phi)
    P[2] = R * np.sin(phi)
    return P

def getPlocFromLocNormal(p_hat):
    phi = np.arcsin(p_hat[2])
    lam = np.arctan2(p_hat[1], p_hat[0])
    return PLoc(np.degrees(lam), np.degrees(phi))

def getVeVnFromAzvel(pLoc, pAzvel): #cartesian Ve and Vn for point, and motion azimuth and magnitude (mm/Y)
    lam, phi = locToRadians(pLoc)
    # unit vectors for 'easterly' and 'northerly' at P
    e_hat = np.array([-np.sin(lam), np.cos(lam), 0.0])
    n_hat = np.array([-np.sin(phi) * np.cos(lam), -np.sin(phi) * np.sin(lam), np.cos(phi)])
    # 2D motion vector at point
    V = np.array([np.sin(np.radians(pAzvel.azimuth)) * pAzvel.vel,
                    np.cos(np.radians(pAzvel.azimuth)) * pAzvel.vel])
    # return scaled velocity in easterly and northerly directions
    return e_hat * V[0], n_hat * V[1]

def find_moments(long_list, lat_list, ve_list, vn_list, pole):
    sum_alpha = 0.0
    avg_alpha = 0.0
    count = len(ve_list) 
    pe_list, pn_list = getSamplePoints(long_list, lat_list, pole)
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
