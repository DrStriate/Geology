import math
earthRadius = 6371.01 * 1000 # m
class GeoHelper:

    def latutideFromDistN(dist): # dist in meters North
        lat = math.atan2(dist, earthRadius) * 180.0 / math.pi
        return lat

    def longitudeFromDist(latitude, dist): # meters East
        latitudeRadians = math.radians(latitude)
        radiusOfParallel = earthRadius * math.cos(latitudeRadians) # m
        longitudeDeltaRadians = dist / radiusOfParallel
        return math.degrees(longitudeDeltaRadians)

