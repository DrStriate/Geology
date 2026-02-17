# First order simulator of Juan de Fuca plate occlusion of YHS in mantle
# takes a date (number of Ma ago) and a longitude and determines ii occluded
import math
from .geo_helper import GeoHelper as gh

class JFP:
    InitialSubductionDate = -45e6  # Years ca
    maxSubductionDepth = 500.0e3   # m
    plateAngle = 45.0              # degrees
    trackingLatitude = 42.0        # degrees N
    trenchlongitude = -125.0       # degrees E
    subductionRate = 28.0          # mm / year (East)

    clampedSubductionEast = maxSubductionDepth * math.tan(plateAngle / 180 * math.pi)

    @staticmethod
    def leadingEdgeLongitude (date):   # date of YHS
        if date < JFP.InitialSubductionDate:
            return JFP.trenchlongitude

        plateMotion = (date - JFP.InitialSubductionDate) * JFP.subductionRate / 1e3 # m
        plateMotionEast = plateMotion * math.cos(JFP.plateAngle / 180 * math.pi)
        plateMotionDepth = plateMotion * math.sin(JFP.plateAngle / 180 * math.pi)

        # clamp at max depth
        if plateMotionDepth > JFP.maxSubductionDepth:
            plateMotionEast = JFP.clampedSubductionEast

        plateLongitudeCoverage = gh.longitudeFromDistE(JFP.trackingLatitude, plateMotionEast)
        plateEdgeLongitude = JFP.trenchlongitude + plateLongitudeCoverage

        # print("date: " + str(date))
        # print("YHSlongitude: " + str(longitude))
        # print("plateEdgeLongitude: " + str(plateEdgeLongitude))

        return plateEdgeLongitude # plate edge west pf hotspot
