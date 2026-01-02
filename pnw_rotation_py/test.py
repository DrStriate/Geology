
from geo_helper import latutideFromDistN, longitudeFromDist
#from qgis.core import QgsApplication


def main():
    # QgsApplication.setPrefixPath("%OSGEO4W_ROOT%/bin/qgis-bin.exe", True)
    # qgs = QgsApplication([], False)
    # qgs.initQgis()

    print("Running unit tests...\n")
    test_latitudeFromDistance()
    test_lonitudeFromDistance()

    # qgs.exitQgis()
#>>>>>>> Stashed changes

def testFunction(val, valSB, tol, name):
    if abs(val - valSB) < tol:
        print (name + ' passed!')
    else:
        print (name + ' failed. Was ' + str(val) + ' sb ' + str(valSB))

def test_latitudeFromDistance():
    distN = 1000 # 1 km
    latSb = 0.00898 # degrees
    lat = latutideFromDistN(distN)
    testFunction(lat, latSb, 0.001, 'test_latitudeFromDistance')

def test_lonitudeFromDistance():
    distN = 1000 # 1 km
    lat = 45
    lonSb = 0.0127 # degrees
    lat = longitudeFromDist(lat, distN)
    testFunction(lat, lonSb, 0.001, 'test_lonitudeFromDistance')

def test_plate_motion():
    plateMotion = PlateMotion()
    initialState = plateMotion.initialize(YHS_lat, YHS_long, NA_Speed, NA_Bearing)
    period = 1e6        # 1 My
    currentState = plateMotion.getNextState(period, None, False, True)
    testFunction(currentState.latitude, 44.67164642221743, 0.0001, 'test_plate_motion.latitude')
    testFunction(currentState.longitude, -110.33020007827903, 0.0001, 'test_plate_motion.longitude')
    testFunction(currentState.vEast, 0.02687, 0.000001, 'test_plate_motion.vEast')
    testFunction(currentState.vNorth, 0.02687, 0.000001, 'test_plate_motion.vNorth')

if __name__ == "__main__":
    main()