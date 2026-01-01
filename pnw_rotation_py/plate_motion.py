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

class PlateMotion:
    def __init__(self):
        self.currentState = PState(0,0,0,0,0)
        self.naPlateVe = 0
        self.naPlateVn = 0

    def initialize(self, initLat, initLong, naSpeed, naBearing): # speed in m/yr, bearing is azimuth degrees
        self.naPlateVn = math.cos(math.radians(naBearing)) * naSpeed # math.cos(247.5) * 46  mm / Y
        self.naPlateVe = math.sin(math.radians(naBearing)) * naSpeed # math.sin(247.5) * 46  mm / Y

        self.currentState = PState(initLong, initLat, 0, 0, 0)
        return self.currentState

    def getNextState(self, deltaT, rotData, applyRotation, applyNaMotion):
        deltaState = PState(0,0,0,0,0)

        closestIdx = -1;
        if (applyRotation and rotData):
            # get the closest rotation entry velocity for current location
            closestIdx = rotData.getClosestRotEntry(self.currentState.longitude, self.currentState.latitude)
            if closestIdx == -1:
                print('No close rotation entry found')
                return None
            rotation = rotData.rotFeatureList[closestIdx]
            deltaState.vEast  = -rotation[2] / 1000.0 # m / yr
            deltaState.vNorth  = -rotation[3] / 1000.0

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





