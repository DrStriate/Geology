import sys
from geo_helper import GeoHelper as gh

def main():
  yhs = {"lat" : 44.43, "lon" :-110.67}
  V = gh.getPoleRotationV(yhs["lat"], yhs["lon"])
  print("V: ", V)
  return

if __name__ == "__main__":
    sys.exit(main())