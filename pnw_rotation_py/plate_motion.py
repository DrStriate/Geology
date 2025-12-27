# Calculates estimated NA plate motion based on combined plate velocity and rotation data
from dataclasses import dataclass
from qgis._core import QgsMessageLog, Qgis
import math
from .geo_helper import longitudeFromDist, latutideFromDistN

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

    def initialize(self, naSpeed, naBearing):
        self.currentState = PState(0,0,0,0)
        self.stateList = []
        self.naPlateVn = math.cos(naBearing) * naSpeed # math.cos(247.5) * 46  mm / Y
        self.naPlateVe = math.sin(naBearing) * naSpeed # math.sin(247.5) * 46  mm / Y

    def getNextState(self, deltaT, rotData, applyRotation, applyNaMotion):
        deltaState = PState(0,0,0,0)

        if (applyRotation):
            # get the closest rotation entry velocity for current location
            rotation = rotData.getClosestRotEntry(self.currentState.longitude, self.currentState.latitude)
            if not rotation:
                QgsMessageLog.logMessage('No close rotation entry found', tag=PlateMotion.__name__, level=Qgis.Info)
                return None
            deltaState.vEast  = rotation.featureAttribute(2)
            deltaState.vNorth  = rotation.featureAttribute(3)

        if (applyNaMotion):
            deltaState.vNorth += self.naPlateVn
            deltaState.vEast += self.naPlateVe

        #scale motion by time and convert distance to lat/long
        deltaState.latitude = latutideFromDistN(deltaState.vNorth * deltaT)
        deltaState.longitude = longitudeFromDist(self.currentState.latitude + deltaState.latitude,
                                                           deltaState.vEast * deltaT)

        #update current state
        self.currentState.latitude = self.currentState.latitude + deltaState.latitude
        self.currentState.longitude = self.currentState.longitude + deltaState.longitude
        self.currentState.vEast = self.currentState.vEast + deltaState.vEast * deltaT
        self.currentState.vNorth = self.currentState.vNorth + deltaState.vNorth * deltaT

        self.stateList.append(self.currentState)

        return self.currentState





