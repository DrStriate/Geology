import math

EarthRadius = 6371.01 * 1000  # m

class GeoHelper:

    def clamp(value, minimum, maximum):
        return max(minimum, min(value, maximum))

    def latutideFromDistN(dist): # dist in meters North
        lat = math.atan2(dist, EarthRadius) * 180.0 / math.pi
        return lat

    def longitudeFromDist(latitude, dist): # meters East
        latitudeRadians = math.radians(latitude)
        radiusOfParallel = EarthRadius * math.cos(latitudeRadians) # m
        longitudeDeltaRadians = dist / radiusOfParallel
        return math.degrees(longitudeDeltaRadians)

