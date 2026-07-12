import numpy as np
from . import geo_helper as gh
from .src.geo_helper import PAdist, PLoc, PDist
from .src import test_utils as tu
from .src import euler_pole_regression as epr
from .src import gauss_newton as gn
from . import path_layer as pl

SF = 1000.0 # conversion from units in km/ma and yrs to get meters (km/ma * yrs / SF = m)

YHS_lat = 44.43           # Yellowstone caldera yhs @ 0 Ma
YHS_long = -110.67
NA_plate_speed = 23.0    # mm / yr (Current) = Adjusted to Owyhee=Humbolt cauldera ~14Ma
NA_plate_azimuth = 241.0  # degrees 

class YhsPath:
  def __init__(self, parent):
    self.parent = parent
    self.yhs_loc = PLoc(YHS_long, YHS_lat)
    self.path_layer_manager = pl.PathLayerManager()
    self.yhs_pole = {'long' : -76.38, 'lat': 49.60, 'omega': 0.774}

    self.NAPAdist = PAdist(NA_plate_azimuth, NA_plate_speed)
    self.pnw_pole_v = {'e': 0, 'n': 0}

    self.trans_path_layer = None
    self.rot_path_layer = None

    self.useGpuData = False
    self.useGpuModel(True)

    self.path_layer_manager = pl.PathLayerManager()
  
  def useGpuModel (self, setTrue):
    if self.useGpuData == setTrue:
      return
    if setTrue:
      self.sample_radius = 600.0 # km
      self.Rot_Data_Sample_center = PLoc(-119.0, 45.0)
      self.pnw_rot_pole, self.pnw_rot_pole_v = self.getRotPoleAndVelocity(self.Rot_Data_Sample_center, self.sample_radius)
    else: # use published pole / velocity info
      self.useGpuData - False
      self.sample_radius = 0
      self.pnw_rot_pole = {"lat" : 45.54,  "long" : -119.60, "omega" : 1.2 } #OC_NA_Pole
      self.pnw_rot_pole_v = PDist(0.0, 5.0)
      label_text1 = f"{self.pnw_rot_pole['long']:.4f}, {self.pnw_rot_pole['lat']:.4f}, {self.pnw_rot_pole['omega']:.3f} deg, "
      label_text2 = f"e: {(self.pnw_rot_pole_v.east / 1E3):.2f} km, n: {(self.pnw_rot_pole_v.north / 1E3):.2f} km"
      self.parent.geoWhiteboard.draw_target(self.pnw_rot_pole['long'], self.pnw_rot_pole['lat'], label_text1 + label_text2)
    self.useGpuData = setTrue
    self.erase_everything()

  def checkLayersCreated(self): 
    self.yhs_path_layer = self.path_layer_manager.getInstance("YHS Path", "red")
    self.trans_path_layer = self.path_layer_manager.getInstance("Pole Trans", "blue")
    self.rot_path_layer = self.path_layer_manager.getInstance("Pole Rot", "black")

  def erase_everything(self): # any statefulness that changes with runs should be resettable
    if self.path_layer_manager is not None:
      self.path_layer_manager.erase_everything()
  
  def closeLayers(self):
    self.path_layer_manager.close_layers()
  
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
    # Plot and label the Euler rotation pole
    self.checkLayersCreated()
    label_text1 = f"{self.pnw_rot_pole['long']:.4f}, {self.pnw_rot_pole['lat']:.4f}, {self.pnw_rot_pole['omega']:.3f} deg, "
    label_text2 = f"e: {(self.pnw_rot_pole_v.east / 1E3):.2f} km, n: {(self.pnw_rot_pole_v.north / 1E3):.2f} km, {self.sample_radius} km"
    self.parent.geoWhiteboard.draw_target(self.pnw_rot_pole['long'], self.pnw_rot_pole['lat'], label_text1 + label_text2)

    # 1: Move yhs loc by NA speed scaled by ma from 0 Ma location (red line)
    #new_yhs_loc1=self.yhs_path_layer.addAnnotationsForPoleRotationOfPoint(self.yhs_loc, self.yhs_pole, yrs/1e6)
    plateAsp = gh.PAdist(self.NAPAdist.azimuth, -self.NAPAdist.dist * yrs / SF)
    new_yhs_loc1 = gh.getPlocForPlocAndPAdist(self.yhs_loc, plateAsp)
    #new_yhs_loc1.print("get_yhs_loc new_yhs_loc1")

    # 2: Move by ma scaled pole translation v (blue line)
    v_dist = -self.pnw_rot_pole_v.magnitude() * yrs / SF
    new_yhs_loc2 = gh.getPointFromAzimuthDistance(new_yhs_loc1, self.pnw_rot_pole_v.azimuth(), v_dist) 
    yhs_v_paths = [
        [(new_yhs_loc1.long, new_yhs_loc1.lat), (new_yhs_loc2.long, new_yhs_loc2.lat), "translation v"]]
    self.trans_path_layer.add_run_paths_to_path_layer(yhs_v_paths)
    new_yhs_loc2.print("get_yhs_loc new_yhs_loc2 ")

    # 3: Rotate by ma scaled pole omega 
    loc = self.rot_path_layer.addAnnotationsForPoleRotationOfPoint(new_yhs_loc2, self.pnw_rot_pole, yrs / 1e6)
    new_yhs_loc3 = gh.getPoleRotationOfPoint(self.pnw_rot_pole, new_yhs_loc2, yrs / 1e6)
    new_yhs_loc3.print("new_yhs_loc3")
    self.parent.geoWhiteboard.draw_target(new_yhs_loc3.long, new_yhs_loc3.lat, "YHS")

    return new_yhs_loc1 #(renders step 1 (new_yhs_loc1))