import numpy as np
from .src import geo_helper as gh
from .src import test_utils as tu
from .src import euler_pole_regression as epr
from .src import gauss_newton as gn

PLoc = gh.PLoc  # 3. Now it is safe to assign your class
PDist = gh.PDist

SF = 1000.0 # conversion from units in km/ma and yrs to get meters (km/ma * yrs / SF = m)

YHS_lat = 44.43           # Yellowstone caldera yhs @ 0 Ma
YHS_long = -110.67
NA_plate_speed = 23.0    # mm / yr (Current) = Adjusted to Owyhee=Humbolt cauldera ~14Ma
NA_plate_bearing = 241.0  # degrees 

class YhsPath:
  def __init__(self, parent):

    self.parent = parent
    self.yhs_loc = PLoc(YHS_long, YHS_lat)
    self.na_plate_v = {'e': np.sin(np.radians(NA_plate_bearing)) * NA_plate_speed, 
                       'n': np.cos(np.radians(NA_plate_bearing)) * NA_plate_speed}
    self.pnw_pole_v = {'e': 0, 'n': 0}
    self.sample_radius = 500.0 # km
    self.Rot_Data_Sample_center = PLoc(-118.0, 45.0)
    
    self.pnw_rot_pole, self.pnw_rot_pole_v = self.getRotPoleAndVelocity(self.Rot_Data_Sample_center, self.sample_radius)
    self.setup_graphics(parent)

  def setup_graphics(self, parent):
    parent.rotDestLayer.dataProvider().truncate()
    parent.yhsRotFeatureList = []

  def reset_data(self): # any statefulness that changes with runs should be resettable
    return

  def draw_translation_vector(self, pStart, pdist):
    feature = self.parent.rotData.createRotFeature(pStart, pdist, 0.001) 
    self.parent.yhsRotFeatureList.append(feature)
    
  def getRotPoleAndVelocity(self, raw_data_center, distance):
    # get rot data
    lat_list, long_list, ve_list, vn_list, se, sn =\
      tu.get_GPS_rotation_data(raw_data_center.lat, raw_data_center.long, distance * 1000)   

    # get data Euler pole from the raw data set
    raw_pole = epr.fit_euler_pole_linear(lat_list, long_list, ve_list, vn_list)

    # apply Gauss-Newton analysis to get any translation (non-rotation) components
    gn_out = gn.solve_gauss_newton_2D_transform_geo(long_list, lat_list, ve_list, vn_list, raw_pole)

    # strip any translation element to get rot-only 
    rot_ve_list = np.array(ve_list) - gn_out['t_x']
    rot_vn_list = np.array(vn_list) - gn_out['t_y']

    # get data Euler pole from the raw data set
    rot_pole = epr.fit_euler_pole_linear(lat_list, long_list, rot_ve_list, rot_vn_list)

    # offset is in meters per ma and we want a rate (km/ma or mm/yr) so we need to scale
    return rot_pole, PDist(gn_out['t_x'] / SF, gn_out['t_y'] / SF)
  
  def get_yhs_loc(self, yrs): # yrs is years
    # 1: Move yhs loc by NA speed scaled by ma from 0 Ma location
    delta_d = PDist(-self.na_plate_v['e'] * yrs, -self.na_plate_v['n'] * yrs)
    ma_yhs_lat, ma_yhs_long = gh.LatLongForDeDn(self.yhs_loc.lat,self.yhs_loc.long, delta_d.east, delta_d.north)
    new_yhs_loc1 = PLoc(ma_yhs_long, ma_yhs_lat)
    new_yhs_loc1.print("get_yhs_loc new_yhs_loc1")

    # 2: Move up by pole velocity scaled by ma from 0 Ma location
    dist_moved_by_pole_v = PDist(self.pnw_rot_pole_v.east * yrs / SF, self.pnw_rot_pole_v.north * yrs / SF)
    new_yhs_loc2 = gh.getPlocForPdist(new_yhs_loc1, dist_moved_by_pole_v)
    new_yhs_loc2.print("get_yhs_loc new_yhs_loc2 ")
    self.draw_translation_vector(new_yhs_loc1, dist_moved_by_pole_v)

    return new_yhs_loc1