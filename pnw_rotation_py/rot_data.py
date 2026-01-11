from qgis._core import QgsMessageLog, Qgis, QgsProject, QgsFeature
from dataclasses import dataclass
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import LinearNDInterpolator

@dataclass
class PState:
    longitude: float
    latitude: float
    vEast: float
    vNorth: float
    rotIdx: int


class RotData:
    rotDataLayerName = 'nshm2023_GPS_velocity'
    XRange = {"min": -125.0, "max": -110.0}
    YRange = {"min": 40.0, "max": 50.0}
    XYSteps = 50
    XSize = (XRange["max"] - XRange["min"]) / XYSteps
    YSize = (YRange["max"] - YRange["min"]) / XYSteps


    def __init__(self):
        self.rotSourceLayer = None
        self.rotDataLoaded = False
        self.rotFieldList = []
        self.rotFeatureList = []
        self.u_interp = None
        self.v_interp = None
        self.interp_func = None
        self.interpFunction = "ClosestEntry"

    def clearData(self):
        self.rotFieldList = []
        self.rotFeatureList = []
        return

    def load(self):
        if self.rotDataLoaded:
            return True

        print ("Rot_d2ata loading ")
        # get rotation data layer
        dataLayers = QgsProject.instance().mapLayersByName(self.rotDataLayerName)
        if len(dataLayers) == 0 :
            QgsMessageLog.logMessage('Could not access source data layer' +
                                     self.rotDataLayerName, tag=RotData.name, level=Qgis.Info)
            return False

        self.rotSourceLayer = dataLayers[0]

        #get rot data fields
        fields = self.rotSourceLayer.fields()
        for field in fields:
            self.rotFieldList.append(field)

        #get rot features
        featureList = self.rotSourceLayer.getFeatures()
        for feature in featureList:
            self.rotFeatureList.append(feature)

        self.rotDataLoaded = True
        return True

    # Setup routines to prep for sampling depending on interpFunction
    def setupSampling(self, interpFunction):
        self.interp_func = interpFunction
        if self.interp_func == "LinearNDInterpolator":
            self.linearInterpSamplerSetup()
        return

    # Sample functions dependent on interp_function setting
    def sample(self, longitude, latitude):
        if self.interpFunction == "ClosestEntry":
            return self.getClosestRotEntry(longitude, latitude)
        elif self.interpFunction == "LinearNDInterpolator":
            return self.getLinearInterpSample(longitude, latitude)
        return None

    def linearInterpSamplerSetup(self):
        # build x (longitude), y (latitude), u (Ve) and v (vn) arrays
        x, y, u, v = [], [], [], []
        for feature in self.rotFeatureList:
            x.append(feature.attribute(0))
            y.append(feature.attribute(1))
            u.append(feature.attribute(2))
            v.append(feature.attribute(3))

        print("setup data for linear interpolation")
        points = np.column_stack((x, y))  # (N, 2)
        values = np.column_stack((u, v))  # (N, 2)
        self.interp_func = LinearNDInterpolator(points, values)

        print("resampleToGrid")
        return

    def getLinearInterpSample(self, longitude, latitude):
        if longitude < self.XRange["min"] or longitude > self.XRange["max"]:
            print("Longitude out of range: " + str(longitude))
            return None
        if latitude < self.YRange["min"] or latitude > self.YRange["max"]:
            print("latitude out of range: " + str(latitude))
            return None

    # returns PState of closest rot vector
    def getClosestRotEntry(self, longitude, latitude):
        if not self.rotDataLoaded:
            if not self.load():
                QgsMessageLog.logMessage('No rotation data loaded', RotData.name, Qgis.Info)
                return None

        closestIdx = -1
        minDist = 1e10
        idx = 0
        for feature in self.rotFeatureList:
            f_lon = feature.attribute(0)
            f_lat = feature.attribute(1)
            dist = (f_lon - longitude)**2 + (f_lat - latitude)**2
            if (dist < minDist):
                minDist = dist
                closestIdx = idx
            idx += 1

        if closestIdx == -1:
            return None

        # self.compareRotFeatureToInterp(self.rotFeatureList[closestIdx])

        return PState(
            self.rotFeatureList[closestIdx].attribute(0),
            self.rotFeatureList[closestIdx].attribute(1),
            self.rotFeatureList[closestIdx].attribute(2),
            self.rotFeatureList[closestIdx].attribute(3),
            closestIdx)

        u, v = self.interp_func(longitude, latitude)

        # note index set to zero - need to create a new QFeature array for interpolated rot display TO DO
        return PState(longitude, latitude, u, v, 0)

    # samping using griddata - not yet complete

    # Google guidance on 'python math library smooth 2d vector field'
    #
    # import numpy as np
    # from scipy.interpolate import griddata
    # import matplotlib.pyplot as plt
    #
    # # Assume you have scattered data points (x, y, u, v)
    # # x = [...], y = [...], u = [...], v = [...]
    #
    # # Create a new, regular mesh grid for the smoothed field
    # xx, yy = np.meshgrid(np.linspace(min(x), max(x), 50),
    #                      np.linspace(min(y), max(y), 50))
    #
    # points = np.transpose(np.vstack((x, y)))
    #
    # # Interpolate the U and V components onto the new grid
    # u_interp = griddata(points, u, (xx, yy), method='cubic')
    # v_interp = griddata(points, v, (xx, yy), method='cubic')
    #
    # # Plot the smoothed vector field
    # # plt.quiver(xx, yy, u_interp, v_interp)
    # # plt.show()

    def interpTest(self):
        x, y, u, v = [], [], [], []
        # for feature in self.rotFeatureList:
        #     compareRotFeatureToInterp(self, feature)
                # x.append(lon)
                # y.append(lat)
                # u.append(veErr)
                # v.append(vnErr)

        # plt.quiver(x, y, u, v)
        # plt.show()
        return

    def compareRotFeatureToInterp (self, feature):
        lon = feature.attribute(0)
        lat = feature.attribute(1)
        ve = feature.attribute(2)
        vn = feature.attribute(3)

        interpState = self.getInterpRotEntry(lon, lat)
        if interpState is not None:
            lonErr = lon - interpState.longitude
            latErr = lat - interpState.latitude
            veErr = ve - interpState.vEast
            vnErr = vn - interpState.vNorth

            print(
                "veErr: " + str(veErr) + " (" + str(ve) + " raw vs " + str(interpState.vEast) + " interp)\n" +
                "vnErr: " + str(vnErr) + " (" + str(vn) + " raw vs " + str(interpState.vNorth) + " interp)")
        else:
            print("No interp found")
        return