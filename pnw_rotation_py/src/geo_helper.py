import numpy as np
from pyproj import Geod
from dataclasses import dataclass, field

R = 6371.0 # Earth radius in km

geod = Geod(ellps="WGS84")

def setGeod(realWorld):
    global geod  
    if realWorld:
        geod = Geod(ellps="WGS84")
    else:
        geod = Geod(a=R*1e3, b=R*1e3)

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
    phi: float = field(init=False)
    lam: float = field(init=False)

    def __post_init__(self):
        self.phi = np.radians(self.lat)
        self.lam = np.radians(self.long)

    def print(self, label = ""): 
        print (f"{label} long: {self.long:0.3f}, lat:  {self.lat:0.3f}")
    
    def __iadd__(self, other):  # Overriding the += operator
        if isinstance(other, PLoc):
            self.long += other.long
            self.lat += other.lat
            self.lam += other.lam
            self.phi += other.phi 
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
# Note that distance measurements in plate kinematics should be in degrees, not kilometers.
class PAvel:
    azimuth: float
    vel: float
    def print(self, label = ""): 
        print (f"{label} azimuth: {self.azimuth:0.3f}, vel:  {self.vel:0.3f}")
    @classmethod
    def from_V(cls, v2) -> 'PAvel': # V2 is [east_v, north_v]
        return cls(
            np.degrees(np.arctan2(v2[0], v2[1])),
            np.hypot(v2[0], v2[1]))

@dataclass
class EulerPole:
    long: float
    lat: float
    omega: float
    ploc: PLoc = field(init=False)

    def __post_init__(self):
        self.ploc = PLoc(self.long, self.lat)

    is_clockwise: bool = False  # Set to True for clockwise poles needing antipodal shift

    def print(self, label: str = ""):
        print(f"{label} long: {self.long:0.3f}, lat: {self.lat:0.3f}, omega: {self.omega:.6f}")

    def normalize(self) -> 'EulerPole':
        """
        Ensures right-hand rule math produces the correct surface motion.
        If rotation is clockwise, returns a NEW EulerPole flipped to its 
        antipodal position with is_clockwise=False so standard CCW vector math works.
        """
        if self.is_clockwise:
            # Shift longitude by 180 and keep within [-180, 180]
            new_long = (self.long + 180) % 360
            if new_long > 180:
                new_long -= 360

            return EulerPole(
                long=new_long,
                lat=-self.lat,
                omega=self.omega,
                is_clockwise=False
            )
        return self        

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

# returns km
def getNortherlyEasterlyFromLatLongPoints(lon1, lat1, lon2, lat2):
    # forward_azimuth is the angle from point 1 to point 2 (degrees clockwise from North)
    forward_azimuth, back_azimuth, distance_meters = geod.inv(lon1, lat1, lon2, lat2)

    # Convert azimuth to radians
    azimuth_rad = np.radians(forward_azimuth)
    
    # Calculate components - actually calculating 
    northerly_km = distance_meters * np.cos(azimuth_rad) * 0.001
    easterly_km = distance_meters * np.sin(azimuth_rad) * 0.001
    return northerly_km, easterly_km

def getFwdAzimuthFromLocations (point1, point2): #both PLocs
   # forward_azimuth is the angle from point 1 to point 2 (degrees clockwise from North)
    forward_azimuth, back_azimuth, distance_meters = geod.inv(point1.long, point1.lat, point2.long, point2.lat)
    return forward_azimuth % 360

# Great circle distance (km)
def getDistanceBetweenPoints(point1, point2): #both PLocs
    forward_azimuth, back_azimuth, distance_meters = geod.inv(point1.long, point1.lat, point2.long, point2.lat)
    return distance_meters / 1000.0

def create_sample (start_lon, start_lat, azimuth, distance): # distance in km!
    # Calculate the terminus point
    end_lon, end_lat, back_azimuth = geod.fwd(
        start_lon, 
        start_lat, 
        azimuth, 
        distance * 1000.0) 
    return PLoc(end_lon, end_lat)

def clamp(value, minimum, maximum):
    return max(minimum, min(value, maximum))

def normalize(vect):
    mag = np.linalg.norm(vect)
    if mag > 0:
        return vect / mag
    return vect

# Code below probably needs to be refactored to use code/methods above which are more accurate

# gets lat and long converted to coordinate distances from pole. This is an approximation 
def getSamplePoints(long_list, lat_list, center_ploc):
  p_e = np.zeros(len(long_list))
  p_n = np.zeros(len(long_list))
  for i in range(len(long_list)):
    # convert sample points to meters
    p_n[i], p_e[i] = getNortherlyEasterlyFromLatLongPoints(center_ploc.long, center_ploc.lat, long_list[i], lat_list[i])
  return p_e, p_n 

# Be wary of use of these distance metrics. 

def latitudeFromDistN(dist): # dist in km North
    lat = np.arctan2(dist, R) * 180.0 / np.pi
    return lat

def longitudeFromDistE(latitude, dist): # km East
    latitudeRadians = np.radians(latitude)
    radiusOfParallel = R * np.cos(latitudeRadians) # m
    longitudeDeltaRadians = dist / radiusOfParallel
    return np.degrees(longitudeDeltaRadians)



