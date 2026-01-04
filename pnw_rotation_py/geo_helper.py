import math
earthRadius = 6371.01 * 1000 # m
class GeoHelper:

    def latutideFromDistN(distN): # dist in meters North
       lat = math.atan2(distN, earthRadius) * 180.0 / math.pi
       return lat

    def longitudeFromDist(latitude, distE): # meters East
        latitudeRadians = math.radians(latitude)
        radiusOfParallel = earthRadius * math.cos(latitudeRadians) # m
        longitudeDeltaRadians = distE / radiusOfParallel
        return math.degrees(longitudeDeltaRadians)

