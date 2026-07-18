import numpy as np
from dataclasses import dataclass
from .src import geo_helper as gh
from .src.geo_helper import PAvel, PLoc, PVel, EulerPole
from .src import test_utils as tu
from .src import euler_pole_regression as epr
from .src import gauss_newton as gn
from .src import euler_kinematics as ek
from . import path_layer as pl

SF = 1000.0 # conversion from units in km/ma and yrs to get meters (km/ma * yrs / SF = m)

YHS_lat = 44.43           # Yellowstone caldera yhs @ 0 Ma
YHS_long = -110.67
NA_plate_speed = 23.0    # mm / yr (Current) = Adjusted to Owyhee=Humbolt cauldera ~14Ma
NA_plate_azimuth = 241.0  # degrees 

@dataclass
class YhsPropertyBag:
  NAPAvel: PAvel
  PnwVPAvel: PAvel
  PnwRotPole: EulerPole

class YhsPath:
  def __init__(self, parent):
    self.parent = parent
    self.path_layer_manager = pl.PathLayerManager()
    
    self.yhsLoc = PLoc(YHS_long, YHS_lat)

    # NA Plate 
    self.NAPAvel = PAvel(NA_plate_azimuth, NA_plate_speed)
    self.NAPole = None # Set up from above NAPADist in setupEulerPoles()
    #self.yhs_pole = EulerPole( -76.38, 49.60,  0.774 ) ?/

    # PNW Rotation 
    self.PnwRotPole = None # Set up from above NAPADist in setupEulerPoles()

    # PNW Translation 
    self.PnwVPAvel = PAvel(0, 0)
    self.PnwVPole = None # Set up from above NAPADist in setupEulerPoles()

    ###

    self.PnwVPoleLayer = None
    self.PnwRotPoleLayer = None
    self.NaPoleLayer = None

    self.useGpsData = False
    self.pole_model = 2 # 1 is NA then Pole-V then Pole-R, 2 is NA - Translated-R - Pole-V
    self.setupEulerPoles(True)
  
  def getYhsPropertyBag(self):
    propertyBag = YhsPropertyBag(self.NAPAvel, self.PnwVPAvel, self.PnwRotPole)
    return propertyBag
  
  def setupEulerPoles (self, useGpsData):
    if self.useGpsData == useGpsData:
      return
    
    if useGpsData:
      sample_radius = 600.0 # km
      self.Rot_Data_Sample_center = PLoc(-119.0, 45.0)
      self.PnwRotPole, self.PnwVPAvel = self.getRotPoleAndVelocity(self.Rot_Data_Sample_center, sample_radius)

    else: # use published pole / velocity info
      self.PnwRotPole = EulerPole(-119.60, 45.54, 1.2 ) #OC_NA_Pole
      self.PnwVPAvel= PVel(0.0, 5.0)

    self.useGpsData = useGpsData
    self.erase_everything()
    #self.display_rot_pole_info()

  def checkLayersCreated(self): 
    self.NaPoleLayer = self.path_layer_manager.getInstance("NA pole", "red")
    self.PnwVPoleLayer = self.path_layer_manager.getInstance("Pnw V pole ", "blue")
    self.PnwRotPoleLayer = self.path_layer_manager.getInstance("Pnw Rot pole", "black")

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
    lat_list, long_list, ve_list, vn_list, we_list, wn_list =\
      tu.get_GPS_rotation_data(raw_data_center.long, raw_data_center.lat, distance * 1000)   

    # get data Euler pole from the raw data set
    raw_pole = epr.fit_euler_pole_linear(lat_list, long_list, ve_list, vn_list)# we_list, wn_list)

    # apply Gauss-Newton analysis to get any translation (non-rotation) components
    gn_out = gn.solve_gauss_newton_2D_transform_geo_wtd(long_list, lat_list, ve_list, vn_list, we_list, wn_list, raw_pole)

    # strip any translation element to get rot-only 
    rot_ve_list = np.array(ve_list) - gn_out['t_x']
    rot_vn_list = np.array(vn_list) - gn_out['t_y']

    # get data Euler pole from the raw data set
    rot_pole = epr.fit_euler_pole_linear(lat_list, long_list, rot_ve_list, rot_vn_list)

    # get Velocity pole PAVel info 
    pnwVPAVel = gh.getPAvel(gn_out['t_x'] / SF, gn_out['t_y'] / SF)

    # offset is in meters per ma and we want a rate (km/ma or mm/yr) so we need to scale
    return rot_pole, pnwVPAVel

  # Plot and label the Euler rotation pole
  def display_rot_pole_info(self):
    label_text1 = f"{self.PnwRotPole.long:.4f}, {self.PnwRotPole.lat:.4f}, {self.PnwRotPole.omega:.3f} deg, "
    label_text2 = f"az: {self.NAPAvel.azimuth:.2f} deg, v: {(self.NAPAvel.vel):.2f} km/ma"
    self.parent.geoWhiteboard.draw_target(self.PnwRotPole.long, self.PnwRotPole.lat, label_text1 + label_text2)
  
  def get_yhs_loc(self, yrs): # yrs is years
    self.checkLayersCreated()

    # 1: Move yhs loc by NA speed scaled by ma from 0 Ma location (red line)
    self.NAPole = ek.getEulerPoleFromPlocAndPazvel(self.yhsLoc, self.NAPAvel)
    loc_1 = self.NaPoleLayer.addAnnotationsForPoleRotationOfPoint(self.yhsLoc, self.NAPole, -yrs/1e6)
    #loc_1.print("loc_1")

    # get the PnwVPole given it varies with loc_1 (based on NA movement)
    self.PnwVPole = ek.getEulerPoleFromPlocAndPazvel(loc_1, self.PnwVPAvel)

    if self.pole_model == 1: # move YHS loc by PnwVPole then rotate (most direct model)
      
      # 2: Move by ma scaled pole translation v (blue line) - note pole is position dpendent
      loc_2 = self.PnwVPoleLayer.addAnnotationsForPoleRotationOfPoint(loc_1, self.PnwVPole, -yrs / 1e6)
      #loc_2.print("loc_2: ")

      # 3: Rotate by ma scaled pole omega 
      loc_3 = self.PnwRotPoleLayer.addAnnotationsForPoleRotationOfPoint(loc_2, self.PnwRotPole, yrs / 1e6)

    else: # model 2: move rot pole down by pole_v, rotate loc1_i to loc_2
      # #2: Rotate by pre-translated rot pole 
      new_rot_pole_ploc, disp = ek.getPoleRotationOfPoint(self.PnwVPole, self.PnwRotPole.ploc(), yrs / 1e6)
      t_RotPole = EulerPole(new_rot_pole_ploc.long, new_rot_pole_ploc.lat, self.PnwRotPole.omega)
      self.parent.geoWhiteboard.draw_target(t_RotPole.long, t_RotPole.lat, 
                                            f"{((int))(yrs / 1e6)} Ma pole({t_RotPole.long:0.3f}, {t_RotPole.lat:0.3f})")
      loc_2 = self.PnwRotPoleLayer.addAnnotationsForPoleRotationOfPoint(loc_1, t_RotPole, yrs / 1e6)
      #loc_2.print("loc_2: ")

      # #3: Translate back up 
      loc_3 = self.PnwVPoleLayer.addAnnotationsForPoleRotationOfPoint(loc_2, self.PnwVPole, -yrs / 1e6)

    self.parent.geoWhiteboard.draw_target(loc_3.long, loc_3.lat, f"YHS ({loc_3.long:0.3f}, {loc_3.lat:0.3f})")
    #loc_3.print("loc_3: ")
    return loc_3