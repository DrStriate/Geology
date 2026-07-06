import numpy as np
from .src import geo_helper as gh

PLoc = gh.PLoc  # 3. Now it is safe to assign your class
PDist = gh.PDist

YHS_lat = 44.43           # Yellowstone caldera yhs @ 0 Ma
YHS_long = -110.67
NA_plate_speed = 23e-3    # m / yr (Current) = Adjusted to Owyhee=Humbolt cauldera ~14Ma
NA_plate_bearing = 241.0  # degrees azimuth

class YhsPath:
  def __init__(self):
    self.yhs_loc = PLoc(YHS_long, YHS_lat)
    self.na_plate_v = {'e': np.sin(np.radians(NA_plate_bearing)) * NA_plate_speed, 
                       'n': np.cos(np.radians(NA_plate_bearing)) * NA_plate_speed}
    self.pnw_rot_pole  = {"lat" : 45.54, "lon" : 119.6, "omega": 1.32} # default NA_OC_WS
    self.pnw_pole_v = {'e': 0, 'n': 0}
  
  def get_yhs_loc(self, ma):
    
    # 1: Move yhs loc by NA speed scaled by aa from 0 Ma location
    delta_d = PDist(-self.na_plate_v['e'] * ma, -self.na_plate_v['n'] * ma)
    ma_yhs_lat, ma_yhs_long = gh.LatLongForDeDn(self.yhs_loc.lat,self.yhs_loc.long, delta_d.east, delta_d.north)
    retval = PLoc(ma_yhs_long, ma_yhs_lat)
    retval.print("get_yhs_loc retval ")
    return retval
  


    