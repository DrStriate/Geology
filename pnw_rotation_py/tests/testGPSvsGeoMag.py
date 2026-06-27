# The following code stub will make the 'run Pyrthon file from a subfolder' work,
# in VS Code but it is not recommended to use this in production code.
# from config_root import setRoot
# setRoot()

import math
import csv
import os

from rot_data import RotData
from src.geo_helper import GeoHelper as gh

def test():
    rotData = RotData()
    rotData.load()
    rotData.setupSampling("LinearNDInterpolator")

    lastT = 0
    lastR = 0
    plugin_dir = os.path.dirname(__file__)
    data_path = os.path.join(plugin_dir, "testData.csv")
    with open(data_path) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        row_num = 0
        for col in csv_reader:
            #if row_num == 0:
                #print(f'Column names are {", ".join(col)}')
            if row_num > 0:
                print(f'{col[0]} lat: {col[1]}, long: {col[2]}, R {col[3]}, T {col[4]}.')

                lat = float(col[1])
                long = float(col[2])
                R = float(col[3])
                T = float(col[4])

                if (lat > 0 and long < 0): # Skip oaleo_d vs oresent_d comparison if no lat/long
                    rotV = rotData.getRotationV(long, lat, True) # m / yr
                    present_d = math.sqrt(sqr(rotV[0]) + sqr(rotV[1])) # d from v_e and v_n
                    #print("present d", present_d)

                    # calculate movement using R, T and rotation pol
                    paleo_d = calculateGeoMagMove ((lat, long), R, T)
                    #print("paleo d = ", paleo_d)

                    # print("\tRatio paleo_d/present_d = ", paleo_d/present_d)

                omega = (R - lastR)/(T - lastT)
                print("\tomega_", row_num, ": ", omega)

                lastT = T
                lastR = R

            row_num += 1

    omegaAv = lastR / lastT
    print("omega avg: ", omegaAv)
    return

def sqr(x) :
    return x * x

def calculateGeoMagMove (point, R, T): #r in degrees, T in Ma returns m/yr
    OcNaPoint = (45.54, -119.60) # rot axis lat lon
    r = gh.DistanceFromLatLong(point, OcNaPoint) * 1000 # m
    d = r * math.atan(math.radians(R/(T*1e6))) # m / yr
    return d


