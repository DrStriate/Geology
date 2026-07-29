import numpy as np
from dataclasses import dataclass, asdict, replace
from .src import euler_pole_regression as epr
import json
import dacite 

from .src.geo_helper import PAvel, PLoc, PVel, EulerPole
from .src import test_utils as tu
from .src import euler_pole_regression as epr
from .src import gauss_newton as gn
from .src import euler_kinematics as ek
from . import path_layer as pl

SF = 1000.0 # conversion from units in km/ma and yrs to get meters (km/ma * yrs / SF = m)

YHS_lat = 44.43           # Yellowstone caldera yhs @ 0 Ma
YHS_long = -110.67

NA_plate_speed = 23.0     # mm / yr (Current) = Adjusted to Owyhee=Humbolt cauldera ~14Ma
NA_plate_azimuth = 241.0  # degrees 

@dataclass
class YhsPropertyBag:
  NAPAvel: PAvel
  PnwVPAvel: PAvel
  PnwRotPole: EulerPole

class YhsPath:
  PrintAlignmentErrors = False
  def __init__(self, parent):
    self.parent = parent
    self.path_layer_manager = pl.PathLayerManager()
    
    self.PrintAlignmentErrors = False
    # self.cross = None # Test PrintAlignmentErrors vector

    ### Data structures needed 
    self.yhs_loc = PLoc(YHS_long, YHS_lat)

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
    self.PnwComboLayer = None

    self.useGpsData = False
    self.pole_model = 2 # 1 is NA then Pole-V then Pole-R, 2 is NA - Translated-R - Pole-V
    self.setupEulerPoles(True)

    self.delta_ve = 0
    self.delta_vn = 0

  # Serialization (Object -> JSON) 
  def get_serialize_bag(self) -> str:
    # asdict recursively handles nesting automatically
    return json.dumps(asdict(self.getYhsPropertyBag()), indent=4)

  # Deserialization (JSON -> Object) 
  def deserialize_and_set_bag(self, json_str: str) -> YhsPropertyBag:
    data_dict = json.loads(json_str)
    
    # dacite automatically looks at the type hints (e.g., NAPAvel: PAvel)
    # and recursively instantiates the inner classes for you.
    propertyBag = dacite.from_dict(data_class=YhsPropertyBag, data=data_dict)
    self.NAPAvel = propertyBag.NAPAvel
    self.PnwVPAvel = propertyBag.PnwVPAvel
    self.PnwRotPole = propertyBag.PnwRotPole  
  
  def getYhsPropertyBag(self):
    propertyBag = YhsPropertyBag(self.NAPAvel, self.PnwVPAvel, self.PnwRotPole)
    return propertyBag
  
  def setYhsPropertyBag(self, propertyBag):
    self.NAPAvel = propertyBag.NAPAvel
    self.PnwVPAvel = propertyBag.PnwVPAvel
    self.PnwRotPole = propertyBag.PnwRotPole
    self.useGpsData = False
    self.setupEulerPoles(False)

  def setupEulerPoles (self, useGpsData):
    if self.useGpsData == useGpsData:
      return
    if useGpsData:
      self.PnwRotPole, self.PnwVPAvel = tu.getPnwGpsRotPoleAndVelocity()

    self.useGpsData = useGpsData
    self.erase_everything()

  def getPnwGpsRotPoleAndVelocity(self):
      self.PnwRotPole, self.PnwVPAvel = tu.getPnwGpsRotPoleAndVelocity()

  def checkLayersCreated(self): 
    self.NaPoleLayer = self.path_layer_manager.getInstance("NA pole", "red")
    self.PnwVPoleLayer = self.path_layer_manager.getInstance("Pnw V pole ", "blue")
    self.PnwRotPoleLayer = self.path_layer_manager.getInstance("Pnw Rot pole", "black")
    self.PnwComboLayer = self.path_layer_manager.getInstance("Pnw R-V pole", "green")

  def erase_everything(self): # any statefulness that changes with runs should be resettable
    if self.path_layer_manager is not None:
      self.path_layer_manager.erase_everything()
  
  def closeLayers(self):
    self.path_layer_manager.close_layers()

  def modeSet(self, poleModel):
    self.pole_model = poleModel

  def setDefaultNAPole(self):
    self.NAPAvel = PAvel(NA_plate_azimuth, NA_plate_speed)

  # Plot and label the NA Euler rotation pole
  def display_NA_pole_info(self):
    label_text1 = f"0 Ma YHS ({self.yhs_loc.long:.3f}, {self.yhs_loc.lat:.3f}), "
    label_text2 = f"az: {self.NAPAvel.azimuth:.1f} deg, v: {(self.NAPAvel.vel):.1f} km/Ma"
    self.parent.geoWhiteboard.draw_target(self.yhs_loc.long, self.yhs_loc.lat, label_text1 + label_text2)
    self.parent.geoWhiteboard.draw_target(self.PnwRotPole.long, self.PnwRotPole.lat, 
                                 f"0 Ma pole ({self.PnwRotPole.long:0.3f}, {self.PnwRotPole.lat:0.3f})")
    
  def get_yhs_loc(self, currentMa, deltaMa, steps): # yrs is years
    self.checkLayersCreated()

    # get a local rot pole we can move as needed 
    runPnwRotPole = replace(self.PnwRotPole)
    # get the PnwVPole which is time invariant 
    self.PnwVPole = ek.getEulerPoleFromPlocAndPavel(runPnwRotPole.ploc(), self.PnwVPAvel)

    # 1: Move yhs loc by NA speed scaled by ma from 0 Ma location (red line)
    self.NAPole = ek.getEulerPoleFromPlocAndPavel(self.yhs_loc, self.NAPAvel)
    loc_1 = self.NaPoleLayer.RenderPoleMotionForMa(self.yhs_loc, self.NAPole, -currentMa)
    #loc_1.print("loc_1")

    if self.pole_model == 1: # move YHS loc by PnwVPole then rotate (most direct model)

      # 2: Move by ma scaled pole translation v (blue line) - note pole is position dpendent
      loc_2 = self.PnwVPoleLayer.RenderPoleMotionForMa(loc_1, self.PnwVPole, -currentMa)
      #loc_2.print("loc_2: ")

      # 3: Rotate by ma scaled pole omega 
      loc_3 = self.PnwRotPoleLayer.RenderPoleMotionForMa(loc_2, runPnwRotPole, currentMa)
      self.parent.geoWhiteboard.draw_target(loc_3.long, loc_3.lat, f"{currentMa} Ma YHS ({loc_3.long:0.3f}, {loc_3.lat:0.3f})")
      #loc_3.print("loc_3: ")

    elif self.pole_model == 2: # model 2: move rot pole down by pole_v, rotate loc_1 to loc_2

      # 2: Translate rot pole  
      new_rot_pole_ploc = ek.getPoleRotationOfPoint(self.PnwVPole, runPnwRotPole.ploc(), currentMa)[0]
      t_RotPole = EulerPole(new_rot_pole_ploc.long, new_rot_pole_ploc.lat, runPnwRotPole.omega)
      self.parent.geoWhiteboard.draw_target(t_RotPole.long, t_RotPole.lat, 
                                            f"{currentMa} Ma pole ({t_RotPole.long:0.3f}, {t_RotPole.lat:0.3f})")
      # Rotate loc by pre-translated rot pole
      loc_2 = self.PnwRotPoleLayer.RenderPoleMotionForMa(loc_1, t_RotPole, currentMa) # WHY NOW POSITIVE Ma?
      #loc_2.print("loc_2: ")

      # 3: Translate back up 
      loc_3 = self.PnwVPoleLayer.RenderPoleMotionForMa(loc_2, self.PnwVPole, -currentMa)
      self.parent.geoWhiteboard.draw_target(loc_3.long, loc_3.lat, f"{currentMa} Ma YHS ({loc_3.long:0.3f}, {loc_3.lat:0.3f})")
      #loc_3.print("loc_3: ")
    
    else: # pole model 3: run pole from start Ma but plot progress points up to final ma

      # self.parent.geoWhiteboard.draw_target(loc_1.long, loc_1.lat, f"{currentMa} Ma YHS ({loc_1.long:0.3f}, {loc_1.lat:0.3f})")
      loc_3 = self.PnwComboLayer.RenderComboPoleMotionForMa(loc_1, self.PnwVPole, runPnwRotPole, currentMa)
      self.parent.geoWhiteboard.draw_target(loc_3.long, loc_3.lat, f"{currentMa} Ma YHS ({loc_3.long:0.3f}, {loc_3.lat:0.3f})")
      self.PnwRotPoleLayer.RenderAzimuthMarkersforMa(loc_1, self.PnwVPole, runPnwRotPole, currentMa)
  
    return loc_3

  def displayGPSDataAndPoles(self):
    # delta_ve = 0
    # delta_vn = 0
    # get rot data
    diam = 600 # km
    center_lat = 45.0
    center_long = -119.0
    lat_list, long_list, ve_list, vn_list, se, sn =\
      tu.get_GPS_rotation_data(center_long, center_lat, diam * 1000)
    mod_ve_list = np.array(ve_list) - self.delta_ve
    mod_vn_list = np.array(vn_list) - self.delta_vn
    self.delta_ve, self.delta_vn = self.finish_test_setup(
    lat_list, long_list, mod_ve_list, mod_vn_list, diam, self.delta_ve, self.delta_vn, se, sn)

  def finish_test_setup(self, lat_list, long_list, ve_list, vn_list, diam, d_ve = None, d_vn = None, s_e = None, s_n = None ):
    # clear and set rot data in Qgis
    self.parent.rotDestLayer.dataProvider().truncate()
    self.parent.yhsRotFeatureList = []
    for i in range(len(lat_list)):
      feature = self.parent.rotData.createRotFeature(
          PLoc(long_list[i], lat_list[i]), PVel(ve_list[i], vn_list[i]), 0.001)
      self.parent.yhsRotFeatureList.append(feature)        

    # get euler pole and gauss newton results and display
    pole = epr.fit_euler_pole_linear_wtd(lat_list, long_list, ve_list, vn_list, s_e, s_n)
    gn_out = gn.solve_gauss_newton_2D_transform_geo_wtd(long_list, lat_list, ve_list, vn_list, s_e, s_n, pole)

    label_text1 = f"{pole.long:.4f}, {pole.lat:.4f}, {pole.omega:.3f} deg, "
    label_text2 = f"e: {(gn_out['t_x'] / 1E3):.2f} km, n: {(gn_out['t_y'] / 1E3):.2f} km, {diam} km"
    self.parent.geoWhiteboard.draw_target(pole.long, pole.lat, label_text1 + label_text2)
    #print(label_text1 + label_text2)

    # if translation correction added, add a delta_V vector to the target to sho that
    if d_ve != 0.0 or d_vn != 0.0:    
      feature = self.parent.rotData.createRotFeature(
        PLoc(pole.long, pole.long), PVel(self.delta_ve + d_ve, self.delta_vn + d_vn), 0.001) 
      self.parent.yhsRotFeatureList.append(feature)
    return gn_out['t_x'], gn_out['t_y']

  def clearGPSparams(self):
    self.delta_ve = 0.0
    self.delta_vn = 0.0  
