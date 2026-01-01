from qgis._core import QgsMessageLog, Qgis, QgsProject, QgsFeature

class RotData:
    rotDataLayerName = 'nshm2023_GPS_velocity'

    def __init__(self):
        self.rotSourceLayer = None
        self.rotDataLoaded = False
        self.rotFieldList = []
        self.rotFeatureList = []

    def clearData(self):
        self.rotFieldList = []
        self.rotFeatureList = []
        return

    def load(self):
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

        return closestIdx