# Calculates estimated NA plate motion based on combined plate velocity and rotation data
from dataclasses import dataclass
import math
import numpy as np
from scipy.spatial.distance import pdist

from .geo_helper import GeoHelper as gh

# qGis versoion
from .geo_helper import GeoHelper
# test version
#from geo_helper import geoHelper

# Plate motion is measuring the change in position of an inertial reference point (e.g. the YHS) on the surface
# of the NA Plate (lat/long). So a motion of the plate will result in the opposite motion of that point
#
# Velocities of the inertial points on the plate are measured in meters per year while the Lat/Lon are scaled by DeltaT

@dataclass
class PLoc:
    long: float
    lat: float

@dataclass
class PDist:
    East: float
    North: float

class PlateMotion:
    # Historic estimates of rates compared to current (0 Ma)
    # so data is Ma vs Rate Scaling in steps of 5 Ma
    # need to validate and refine these data using current studies
    ScalingMa = [
        1.0, # 0 Ma (current)
        2,0, # 5 Ma
        4.0, # 10 Ma...
        6.0, # 15 Ma
        9.0, # 20 Ma
        10.0,# 25 Ma
        6.0, # 30 Ma
        5.5, # 35 Ma
        5.0, # 40 Ma
        4.5, # 45 Ma
        4.0] # 50 Ma

    def __init__(self):
        self.locYhs = PLoc(0, 0)
        self.locNa = PLoc(0,0)
        self.distRot = PDist(0,0)
        self.naPlateVe = 0
        self.naPlateVn = 0
        self.dataFile = None
        self.blockRotationV = []
        self.NAPlateV = []

    def initialize(self, startT, initLat, initLong, naSpeed, naBearing, interpFunction, dataFile = None): # speed in m/yr, bearing is azimuth degrees
        self.naPlateVn = math.cos(math.radians(naBearing)) * naSpeed # math.cos(247.5) * 46  mm / Y
        self.naPlateVe = math.sin(math.radians(naBearing)) * naSpeed # math.sin(247.5) * 46  mm / Y
        self.currentYr = startT
        self.locYhs = PLoc(initLong, initLat)
        self.locNa = PLoc(initLong, initLat)
        self.interpFunction = interpFunction
        self.dataFile = dataFile
        # if self.dataFile:
        #     dataFile.write("long, lat, Na-e, Na-n, Rot-e, Rot-n, Rot-idx, Delta-e, Delta-n, Delta-long, Delta-lat\n")
        return self.locYhs

    def getNextState(self, deltaT, rotData, applyNaScaling, applyRotation, verbose=False):
        deltaDistNa = PDist(0,0)        # NA Plate component of delta dist
        deltaLocNa = PLoc(0,0)
        deltaDistYhs = PDist(0,0)       # combined delta vector
        deltaLocYhs = PLoc(0,0)

        self.currentYr += deltaT
        motionSense = -1.0 if deltaT > 0 else 1.0

        # first compute the unrotated update for Na (the unrotated YHS path)
        deltaDistNa.North += motionSense * self.naPlateVn * abs(deltaT)
        deltaDistNa.East += motionSense * self.naPlateVe * abs(deltaT)

        deltaLocNa.lat = gh.latitudeFromDistN(deltaDistNa.North)
        deltaLocNa.long = gh.longitudeFromDistE(self.locNa.lat + deltaLocNa.lat, deltaDistNa.East)

        self.locNa.lat +=  deltaLocNa.lat
        self.locNa.long +=  deltaLocNa.long

        # next compute the block rotation change
        if (applyRotation and rotData):
            rotV = rotData.getRotationV(self.locYhs.long, self.locYhs.lat,
                                        self.interpFunction != "ClosestEntry")
            appliedMaScaling = 1.0
            if applyNaScaling and self.currentYr < 0.0:
                scaleIdx = min(-self.currentYr / 5.0e6,len(self.ScalingMa) -1.0)
                appliedMaScaling = np.interp(scaleIdx, np.arange(len(self.ScalingMa)), self.ScalingMa)

            self.distRot.East  = motionSense * rotV[0] * appliedMaScaling  * abs(deltaT)# m / yr
            self.distRot.North = motionSense * rotV[1] * appliedMaScaling  * abs(deltaT)

            if verbose:
                azimuthR = math.atan2(self.distRot.East, self.distRot.North)
                print("Rot De: " + str(self.distRot.East) + ", Dn: " + str(self.distRot.North) + ", az: ", math.degrees(azimuthR))

        # Apply NA Motion
        if verbose:
            azimuthNA = math.atan2(deltaDistNa.East, deltaDistNa.North)
            print("Na Ve: " + str(deltaDistNa.East) + ", Vn: " + str(deltaDistNa.North) + ", az: ", math.degrees(azimuthNA))

        deltaDistYhs.North = deltaDistNa.North + self.distRot.North
        deltaDistYhs.East = deltaDistNa.East + self.distRot.East

        if verbose:
            print ("currentYr: " + str(self.currentYr))
            azimuthS = math.atan2(deltaDistYhs.East, deltaDistYhs.North)
            print("Sum Ve: " + str(deltaDistYhs.East) + ", Vn: " + str(deltaDistYhs.North) + ", az: ", math.degrees(azimuthS))

        #scale motion by time and convert distance to lat/long
        deltaLocYhs.lat = gh.latitudeFromDistN(deltaDistYhs.North)
        deltaLocYhs.long = gh.longitudeFromDistE(self.locYhs.lat + deltaLocYhs.lat,
                                                     deltaDistYhs.East)
        #update current state
        self.locYhs.lat += deltaLocYhs.lat
        self.locYhs.long += deltaLocYhs.long

        #latRange = [47.26, 47.40] # zig-zags in nearest sample run
        #latRange = [46.94, 47.16] # sudden jump right at lat 47.09
        # latRange = [0, 0]
        # if self.dataFile and self.currentState.latitude > latRange[0] and self.currentState.latitude < latRange[1]:
        #     self.dataFile.write(
        #         f"{self.currentState.longitude:.4f}" + ", " +
        #         f"{self.currentState.latitude:.4f}" + ", " +
        #         f"{deltaNa.East:.4f}" + ", " +
        #         f"{deltaNa.North:.4f}" + ", " +
        #         f"{deltaRot.vEast:.4f}" + ", " +
        #         f"{deltaRot.North:.4f}" + ", " +
        #         str(rotEntry.rotIdx) + ", " +
        #         f"{deltaState.vEast:.4f}" + ", " +
        #         f"{deltaState.North:.4f}" + ", " +
        #         f"{deltaState.longitude:.4f}" + ", " +
        #         f"{deltaState.latitude:.4f}" + "\n")
        return self.locYhs







