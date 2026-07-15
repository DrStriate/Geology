import numpy as np
from . import geo_helper as gh
from .src.geo_helper import PAdist, PLoc, PDist, EulerPole
from .src import test_utils as tu
from .src import euler_pole_regression as epr
from .src import gauss_newton as gn
from . import path_layer as pl
from .src import euler_kinematics as ek

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

    self.NAPAdist = PAdist(NA_plate_azimuth, NA_plate_speed)
    self.yhs_pole = EulerPole( -76.38, 49.60,  0.774 )

    self.pnw_rot_pole = None  # set up by setupEulerPoles
    self.pnw_v_pole = None    # ditto

    self.trans_path_layer = None
    self.rot_path_layer = None

    self.useGpsData = False
    self.setupEulerPoles(True)

    self.path_layer_manager = pl.PathLayerManager()

    self.pole_model = 2 # 1 is NA then Pole-V then Pole-R, 2 is NA - Translated-R - Pole-V
  
  def setupEulerPoles (self, useGpsData):
    if self.useGpsData == useGpsData:
      return
    
    if useGpsData:
      sample_radius = 600.0 # km
      self.Rot_Data_Sample_center = PLoc(-119.0, 45.0)
      self.pnw_rot_pole, pnw_v_dist = self.getRotPoleAndVelocity(self.Rot_Data_Sample_center, sample_radius)

    else: # use published pole / velocity info
      self.pnw_rot_pole = EulerPole(-119.60, 45.54, 1.2 ) #OC_NA_Pole
      pnw_v_dist = PDist(0.0, 5.0)

    self.pnw_v_azdist = gh.getPAdist(pnw_v_dist.east, pnw_v_dist.north)
    self.useGpsData = useGpsData
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

  def choosePoleModel(self, selected):
    if selected:
      self.pole_model = 1
    else:
      self.pole_model = 2
  
  def getRotPoleAndVelocity(self, raw_data_center, distance):
    # get rot data
    lat_list, long_list, ve_list, vn_list, se, sn =\
      tu.get_GPS_rotation_data(raw_data_center.long, raw_data_center.lat, distance * 1000)   

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

  # Plot and label the Euler rotation pole
  def display_rot_pole_info(self):
    label_text1 = f"{self.pnw_rot_pole.long:.4f}, {self.pnw_rot_pole.lat:.4f}, {self.pnw_rot_pole.omega:.3f} deg, "
    label_text2 = f"az: {(self.pnw_v_azdist.azimuth / 1E3):.2f} deg, v: {(self.pnw_v_azdist.dist / 1E3):.2f} km/ma"
    self.parent.geoWhiteboard.draw_target(self.pnw_rot_pole.long, self.pnw_rot_pole.lat, label_text1 + label_text2)
  
  def get_yhs_loc(self, yrs): # yrs is years
    self.checkLayersCreated()

    # 1: Move yhs loc by NA speed scaled by ma from 0 Ma location (red line)
    paDistNA = gh.PAdist(self.NAPAdist.azimuth, self.NAPAdist.dist)
    pole_NA = ek.getEulerPoleFromPlocAndPazdist(self.yhs_loc, paDistNA)
    loc_1 = self.yhs_path_layer.addAnnotationsForPoleRotationOfPoint(self.yhs_loc, pole_NA, -yrs/1e6)
    #loc_1.print("loc_1")

    # get the pnw_v_pole given it varies with loc_1 (based on NA movement)
    self.pnw_v_pole = ek.getEulerPoleFromPlocAndPazdist(loc_1, self.pnw_v_azdist)

    if self.pole_model == 1: # move YHS loc by pnw_v_pole then rotate (most direct model)
      
      # 2: Move by ma scaled pole translation v (blue line) - note pole is position dpendent
      loc_2 = self.trans_path_layer.addAnnotationsForPoleRotationOfPoint(loc_1, self.pnw_v_pole, -yrs / 1e6)
      #loc_2.print("loc_2: ")

      # 3: Rotate by ma scaled pole omega 
      loc_3 = self.rot_path_layer.addAnnotationsForPoleRotationOfPoint(loc_2, self.pnw_rot_pole, yrs / 1e6)

    else: # model 2: move rot pole down by pole_v, rotate loc1_i to loc_2
      # #2: Rotate by pre-translated rot pole 
      new_rot_pole_ploc, dist = ek.getPoleRotationOfPoint(self.pnw_v_pole, self.pnw_rot_pole.ploc(), yrs / 1e6)
      t_RotPole = EulerPole(new_rot_pole_ploc.long, new_rot_pole_ploc.lat, self.pnw_rot_pole.omega)
      self.parent.geoWhiteboard.draw_target(t_RotPole.long, t_RotPole.lat, 
                                            f"rot pole({t_RotPole.long:0.3f}, {t_RotPole.lat:0.3f})")
      loc_2 = self.rot_path_layer.addAnnotationsForPoleRotationOfPoint(loc_1, t_RotPole, yrs / 1e6)
      #loc_2.print("loc_2: ")

      # #3: Translate back up 
      loc_3 = self.trans_path_layer.addAnnotationsForPoleRotationOfPoint(loc_2, self.pnw_v_pole, -yrs / 1e6)

    self.parent.geoWhiteboard.draw_target(loc_3.long, loc_3.lat, f"YHS ({loc_3.long:0.3f}, {loc_3.lat:0.3f})")
    #loc_3.print("loc_3: ")
    return loc_3