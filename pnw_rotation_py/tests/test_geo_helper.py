from geo_helper import *
import pytest as pt

def test_getPoleRotationOfPoint():
  testPAdist = PAdist(0.0, 0.0)
  point = PLoc(long=-121.57708451886487, lat=40.84987174388522)
  ma = -50.0
  pole = {'lat': 44.871330515586074,
          'long':-119.8854828943077,
          'omega': 0.5614672520657451}

  pole_center = PLoc(pole['long'], pole['lat'])
  total_angle =  pole['omega'] * ma
  radius1 = getDistanceBetweenPoints(point, pole_center)
  azimuth1 = getFwdAzimuthFromLocations(point, pole_center)
  azimuth2 = azimuth1 - total_angle
  outPoint = getPointFromAzimuthDistance(pole_center, azimuth2, radius1)
  radius2 = getDistanceBetweenPoints(outPoint, pole_center)
  print(f"radius1: {radius1:0.1f}")
  print(f"radius2: {radius2:0.1f}")
  print(f"azimuth2: {azimuth2:0.2f}")
  assert radius1 == pt.approx(radius2, 1e-4)


# Revise this test to validare azimuths from pole using for precise azimuths
# forward_azimuth, back_azimuth, distance_meters = geod.inv(lon2, lat2, lon1, lat1)