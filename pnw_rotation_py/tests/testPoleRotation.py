import sys
# The following code stub will make the 'run Pyrthon file from a subfolder' work,
# in VS Code but it is not recommended to use this in production code.
from config_root import setRoot
setRoot()

from src.geo_helper import GeoHelper as gh

def main():
  yhs = {"lat" : 44.43, "lon" :-110.67}
  V = gh.getPoleRotationV(yhs["lat"], yhs["lon"])
  print("V: ", V)
  return

if __name__ == "__main__":
    sys.exit(main())