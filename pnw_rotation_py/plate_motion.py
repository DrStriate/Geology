# Calculates estimated NA plate motion based on combined plate velocity and rotation data
from dataclasses import dataclass
import math

# qGis versoion
from .geo_helper import geoHelper
# test version
#from geo_helper import geoHelper

# Plate motion is measuring the change in position of an inertial reference point (e.g. the YHS) on the surface
# of the NA Plate (lat/long). So a motion of the plate will result in the opposite motion of that point
#
# Velocities of the inertial points on the plate are measured in meters per year while the Lat/Lon are scaled by DeltaT

@dataclass
class PState:
    longitude: float
    latitude: float
    vEast: float
    vNorth: float
    rotIdx: int

    # Historic estimates of rates compared to current (0 Ma)

class PlateMotion:
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

    def initialize(self, initLat, initLong, naSpeed, naBearing): # speed in m/yr, bearing is azimuth degrees
        self.naPlateVn = math.cos(math.radians(naBearing)) * naSpeed # math.cos(247.5) * 46  mm / Y
        self.naPlateVe = math.sin(math.radians(naBearing)) * naSpeed # math.sin(247.5) * 46  mm / Y
        self.totalT = 0.0
        self.currentState = PState(initLong, initLat, 0, 0, 0)
        return self.currentState

    def getNextState(self, deltaT, rotData, applyNaScaling, applyRotation, applyNaMotion):
        deltaState = PState(0,0,0,0,0)

        closestIdx = -1;
        self.totalT += deltaT
        if (applyRotation and rotData):
            # get the closest rotation entry velocity for current location
            closestIdx = rotData.getClosestRotEntry(self.currentState.longitude, self.currentState.latitude)
            if closestIdx == -1:
                print('No close rotation entry found')
                return None
            appliedMaScaling = 1.0
            if applyNaScaling and self.totalT < 0.0:
                scaleIdx = min(int(-self.totalT / 5.0e6),len(self.ScalingMa) -1)
                appliedMaScaling = self.ScalingMa[scaleIdx]

            print ("totalT: " + str(self.totalT))

            rotation = rotData.rotFeatureList[closestIdx]
            deltaState.vEast  = -rotation[2] / 1000.0 * appliedMaScaling # m / yr
            deltaState.vNorth  = -rotation[3] / 1000.0 * appliedMaScaling

        if (applyNaMotion):
            deltaState.vNorth = deltaState.vNorth - self.naPlateVn
            deltaState.vEast = deltaState.vEast - self.naPlateVe

        #scale motion by time and convert distance to lat/long
        deltaState.latitude = geoHelper.latutideFromDistN(deltaState.vNorth * deltaT)
        deltaState.longitude = geoHelper.longitudeFromDist(self.currentState.latitude + deltaState.latitude,
                                                           deltaState.vEast * deltaT)
        #update current state
        nextState = PState(0,0,0,0,0)

        nextState.latitude = self.currentState.latitude + deltaState.latitude
        nextState.longitude = self.currentState.longitude + deltaState.longitude
        nextState.vEast = deltaState.vEast #self.currentState.vEast + deltaState.vEast
        nextState.vNorth = deltaState.vNorth #self.currentState.vNorth + deltaState.vNorth
        nextState.rotIdx = closestIdx;

        self.currentState = nextState
        return nextState





