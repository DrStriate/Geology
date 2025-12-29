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

class PlateMotion:
    def __init__(self):
        self.currentState = PState(0,0,0,0)
        self.stateList = []
        self.naPlateVe = 0
        self.naPlateVn = 0

    def initialize(self, initLat, initLong, naSpeed, naBearing): # speed in m/yr, bearing is azimuth degrees
        self.naPlateVn = math.cos(math.radians(naBearing)) * naSpeed # math.cos(247.5) * 46  mm / Y
        self.naPlateVe = math.sin(math.radians(naBearing)) * naSpeed # math.sin(247.5) * 46  mm / Y

        self.currentState = PState(initLong, initLat, 0, 0)
        self.stateList = [self.currentState]
        return self.currentState

    def getNextState(self, deltaT, rotData, applyRotation, applyNaMotion):
        deltaState = PState(0,0,0,0)

        if (applyRotation and rotData):
            # get the closest rotation entry velocity for current location
            rotation = rotData.getClosestRotEntry(self.currentState.longitude, self.currentState.latitude)
            if not rotation:
                print('No close rotation entry found')
                return None
            deltaState.vEast  = -rotation.featureAttribute(2)
            deltaState.vNorth  = -rotation.featureAttribute(3)

        if (applyNaMotion):
            deltaState.vNorth -= self.naPlateVn
            deltaState.vEast -= self.naPlateVe

        #scale motion by time and convert distance to lat/long
        deltaState.latitude = geoHelper.latutideFromDistN(deltaState.vNorth * deltaT)
        deltaState.longitude = geoHelper.longitudeFromDist(self.currentState.latitude + deltaState.latitude,
                                                           deltaState.vEast * deltaT)
        #update current state
        nextState = PState(0,0,0,0)

        nextState.latitude = self.currentState.latitude + deltaState.latitude
        nextState.longitude = self.currentState.longitude + deltaState.longitude
        nextState.vEast = self.currentState.vEast + deltaState.vEast
        nextState.vNorth = self.currentState.vNorth + deltaState.vNorth

        self.stateList.append(nextState)
        self.currentState = nextState

        return nextState





