# Calculates estimated NA plate motion based on combined plate velocity and rotation data
from dataclasses import dataclass
import math
import numpy as np
from .geo_helper import GeoHelper as gh
from .rot_data import PState

# qGis versoion
from .geo_helper import GeoHelper
# test version
#from geo_helper import geoHelper

# Plate motion is measuring the change in position of an inertial reference point (e.g. the YHS) on the surface
# of the NA Plate (lat/long). So a motion of the plate will result in the opposite motion of that point
#
# Velocities of the inertial points on the plate are measured in meters per year while the Lat/Lon are scaled by DeltaT

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
        self.currentState = PState(0,0,0,0,0)
        self.naPlateVe = 0
        self.naPlateVn = 0
        self.dataFile = None

    def initialize(self, startT, initLat, initLong, naSpeed, naBearing, interpFunction, dataFile = None): # speed in m/yr, bearing is azimuth degrees
        self.naPlateVn = math.cos(math.radians(naBearing)) * naSpeed # math.cos(247.5) * 46  mm / Y
        self.naPlateVe = math.sin(math.radians(naBearing)) * naSpeed # math.sin(247.5) * 46  mm / Y
        self.currentYr = startT
        self.currentState = PState(initLong, initLat, 0, 0, 0)
        self.interpFunction = interpFunction
        self.dataFile = dataFile
        if self.dataFile:
            dataFile.write("long, lat, Na-e, Na-n, Rot-e, Rot-n, Rot-idx, Delta-e, Delta-n, Delta-long, Delta-lat\n")
        return self.currentState

    def getNextState(self, deltaT, rotData, applyNaScaling, applyRotation, verbose=False):
        deltaRot = PState(0,0,0,0,0)    # rotation component of delta
        deltaNa = PState(0,0,0,0,0)     # NA Plate component of delta
        deltaState = PState(0,0,0,0,0)  # combined delta vector

        self.currentYr += deltaT
        #verbose = True
        motionSense = -1.0 if deltaT > 0 else 1.0
        if (applyRotation and rotData):
            appliedMaScaling = 1.0
            if applyNaScaling and self.currentYr < 0.0:
                scaleIdx = min(-self.currentYr / 5.0e6,len(self.ScalingMa) -1.0)
                appliedMaScaling = np.interp(scaleIdx, np.arange(len(self.ScalingMa)), self.ScalingMa)

            if self.interpFunction == "ClosestEntry":
                # get the closest rotation entry velocity for current location
                rotEntry = rotData.getClosestRotEntry(self.currentState.longitude, self.currentState.latitude)

            else: # Interpolated
                rotEntry = rotData.getLinearInterpSample(self.currentState.longitude, self.currentState.latitude)

            if rotEntry is None or rotEntry.rotIdx == -1:
                print('No close rotation entry found')
                return

            deltaRot.vEast  = motionSense * rotEntry.vEast / 1000.0 * appliedMaScaling # m / yr
            deltaRot.vNorth = motionSense * rotEntry.vNorth / 1000.0 * appliedMaScaling

            if verbose:
                azimuthR = math.atan2(deltaRot.vEast, deltaRot.vNorth)
                print("Rot Ve: " + str(deltaRot.vEast) + ", Vn: " + str(deltaRot.vNorth) + ", az: ", math.degrees(azimuthR))

        # Apply NA Motion
        deltaNa.vNorth += motionSense * self.naPlateVn
        deltaNa.vEast += motionSense * self.naPlateVe

        if verbose:
            azimuthNA = math.atan2(deltaNa.vEast, deltaNa.vNorth)
            print("Na Ve: " + str(deltaNa.vEast) + ", Vn: " + str(deltaNa.vNorth) + ", az: ", math.degrees(azimuthNA))

        deltaState.vNorth = deltaNa.vNorth + deltaRot.vNorth
        deltaState.vEast = deltaNa.vEast + deltaRot.vEast

        if verbose:
            print ("currentYr: " + str(self.currentYr))
            azimuthS = math.atan2(deltaState.vEast, deltaState.vNorth)
            print("Sum Ve: " + str(deltaState.vEast) + ", Vn: " + str(deltaState.vNorth) + ", az: ", math.degrees(azimuthS))

        #scale motion by time and convert distance to lat/long
        deltaState.latitude = gh.latutideFromDistN(deltaState.vNorth * abs(deltaT))
        deltaState.longitude = gh.longitudeFromDist(self.currentState.latitude + deltaState.latitude,
                                                           deltaState.vEast * abs(deltaT))
        #update current state
        nextState = PState(0,0,0,0,0)

        nextState.latitude = self.currentState.latitude + deltaState.latitude
        nextState.longitude = self.currentState.longitude + deltaState.longitude
        nextState.vEast = deltaRot.vEast #used to show rot influence
        nextState.vNorth = deltaRot.vNorth
        #nextState.rotIdx = rotEntry.rotIdx

        #latRange = [47.26, 47.40] # zig-zags in nearest sample run
        #latRange = [46.94, 47.16] # sudden jump right at lat 47.09
        latRange = [0, 0]
        if self.dataFile and self.currentState.latitude > latRange[0] and self.currentState.latitude < latRange[1]:
            self.dataFile.write(
                f"{self.currentState.longitude:.4f}" + ", " +
                f"{self.currentState.latitude:.4f}" + ", " +
                f"{deltaNa.vEast:.4f}" + ", " +
                f"{deltaNa.vNorth:.4f}" + ", " +
                f"{deltaRot.vEast:.4f}" + ", " +
                f"{deltaRot.vNorth:.4f}" + ", " +
                str(rotEntry.rotIdx) + ", " +
                f"{deltaState.vEast:.4f}" + ", " +
                f"{deltaState.vNorth:.4f}" + ", " +
                f"{deltaState.longitude:.4f}" + ", " +
                f"{deltaState.latitude:.4f}" + "\n")
        self.currentState = nextState
        return nextState







