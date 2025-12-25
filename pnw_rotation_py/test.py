#import geoHelper
import geoHelper
from qgis.core import QgsApplication
#import sys

def main():
    QgsApplication.setPrefixPath("%OSGEO4W_ROOT%/bin/qgis-bin.exe", True)
    qgs = QgsApplication([], False)
    qgs.initQgis()

    print("Running unit tests...\n")
    test_latitudeFromDistance()
    test_lonitudeFromDistance()

    qgs.exitQgis()

def testFunction(val, valSB, tol, name):
    if abs(val - valSB) < tol:
        print ('Test ' + name + ' passed!')
    else:
        print ('Test ' + name + ' failed. Was ' + str(val) + ' sb ' + str(valSB))

def test_latitudeFromDistance():
    distN = 1000 # 1 km
    latSb = 0.00898 # degrees
    lat = geoHelper.latutideFromDistN(distN)
    testFunction(lat, latSb, 0.001, 'test_latitudeFromDistance')

def test_lonitudeFromDistance():
    distN = 1000 # 1 km
    lat = 45
    lonSb = 0.0127 # degrees
    lat = geoHelper.longitudeFromDist(lat, distN)
    testFunction(lat, lonSb, 0.001, 'test_lonitudeFromDistance')

if __name__ == "__main__":
    main()