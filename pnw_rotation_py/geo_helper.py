import math

EarthRadius = 6371.01 * 1000  # m

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



