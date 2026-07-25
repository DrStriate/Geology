import numpy as np
from dataclasses import dataclass, asdict
import json
import dacite 

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
  PrintAlignmentErrors = False
  def __init__(self, parent):
    self.parent = parent
    self.path_layer_manager = pl.PathLayerManager()
    
    self.PrintAlignmentErrors = False
    self.cross = None # Test PrintAlignmentErrors vector

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
      self.setGpsPoleModel()

    self.useGpsData = useGpsData
    self.erase_everything()

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

  def setGpsPoleModel(self):
      sample_radius = 600.0 # km
      self.Rot_Data_Sample_center = PLoc(-119.0, 45.0)
      self.PnwRotPole, self.PnwVPAvel = self.getRotPoleAndVelocity(self.Rot_Data_Sample_center, sample_radius)

  def getRotPoleAndVelocity(self, raw_data_center, distance):
    # get rot data
    lats, longs, ves, vns, wes, wns=\
      tu.get_GPS_rotation_data(raw_data_center.long, raw_data_center.lat, distance * 1000)   

    rot_pole, pnwVPAVel = epr.extractEulerPoleUsingCombinedRegressions(lats, longs, ves, vns, wes, wns)
    return rot_pole, pnwVPAVel

  # Plot and label the NA Euler rotation pole
  def display_NA_pole_info(self):
    label_text1 = f"0 Ma YHS ({self.yhs_loc.long:.3f}, {self.yhs_loc.lat:.3f}), "
    label_text2 = f"az: {self.NAPAvel.azimuth:.1f} deg, v: {(self.NAPAvel.vel):.1f} km/Ma"
    self.parent.geoWhiteboard.draw_target(self.yhs_loc.long, self.yhs_loc.lat, label_text1 + label_text2)
    self.parent.geoWhiteboard.draw_target(self.PnwRotPole.long, self.PnwRotPole.lat, 
                                 f"0 Ma pole ({self.PnwRotPole.long:0.3f}, {self.PnwRotPole.lat:0.3f})")
    
  def get_yhs_loc(self, currentMa, deltaMa, steps): # yrs is years
    self.checkLayersCreated()

    # 1: Move yhs loc by NA speed scaled by ma from 0 Ma location (red line)
    self.NAPole = ek.getEulerPoleFromPlocAndPavel(self.yhs_loc, self.NAPAvel)
    loc_1 = self.NaPoleLayer.RenderPoleMotionForMa(self.yhs_loc, self.NAPole, -currentMa)
    #loc_1.print("loc_1")

    # get the PnwVPole given its loc covaries with Rot Pole 
    self.PnwVPole = ek.getEulerPoleFromPlocAndPavel(self.PnwRotPole.ploc(), self.PnwVPAvel)

    if self.pole_model == 1: # move YHS loc by PnwVPole then rotate (most direct model)
      
      # 2: Move by ma scaled pole translation v (blue line) - note pole is position dpendent
      loc_2 = self.PnwVPoleLayer.RenderPoleMotionForMa(loc_1, self.PnwVPole, -currentMa)
      #loc_2.print("loc_2: ")

      # 3: Rotate by ma scaled pole omega 
      loc_3 = self.PnwRotPoleLayer.RenderPoleMotionForMa(loc_2, self.PnwRotPole, currentMa)
      self.parent.geoWhiteboard.draw_target(loc_3.long, loc_3.lat, f"{currentMa} Ma YHS ({loc_3.long:0.3f}, {loc_3.lat:0.3f})")
      #loc_3.print("loc_3: ")

    elif self.pole_model == 2: # model 2: move rot pole down by pole_v, rotate loc_1 to loc_2
      # 2: Translate rot pole and rotate loc by pre-translated rot pole 
      new_rot_pole_ploc = ek.getPoleRotationOfPoint(self.PnwVPole, self.PnwRotPole.ploc(), currentMa)[0]

      # get the PnwVPole given its loc covaries with Rot Pole 
      self.PnwVPole = ek.getEulerPoleFromPlocAndPavel(self.PnwRotPole.ploc(), self.PnwVPAvel)

      if (self.PrintAlignmentErrors is True):
        if self.cross is not None:  # r1 dot (r2 cross r3) == 0? Test big circle alignment
          print(f"error m = {(gh.R * np.dot(ek.getRVector(new_rot_pole_ploc), self.cross))}")
        self.cross = np.cross(ek.getRVector(new_rot_pole_ploc), ek.getRVector(self.PnwRotPole.ploc()))

      t_RotPole = EulerPole(new_rot_pole_ploc.long, new_rot_pole_ploc.lat, self.PnwRotPole.omega)
      self.parent.geoWhiteboard.draw_target(t_RotPole.long, t_RotPole.lat, 
                                            f"{currentMa} Ma pole ({t_RotPole.long:0.3f}, {t_RotPole.lat:0.3f})")
      loc_2 = self.PnwRotPoleLayer.RenderPoleMotionForMa(loc_1, t_RotPole, currentMa) # WHY NOW POSITIVE Ma?
      #loc_2.print("loc_2: ")

      # 3: Translate back up 
      loc_3 = self.PnwVPoleLayer.RenderPoleMotionForMa(loc_2, self.PnwVPole, -currentMa)
      self.parent.geoWhiteboard.draw_target(loc_3.long, loc_3.lat, f"{currentMa} Ma YHS ({loc_3.long:0.3f}, {loc_3.lat:0.3f})")
      #loc_3.print("loc_3: ")
    
    else: # pole model 3: run pole from start Ma but progress that point up to final ma
      # self.parent.geoWhiteboard.draw_target(loc_1.long, loc_1.lat, f"{currentMa} Ma YHS ({loc_1.long:0.3f}, {loc_1.lat:0.3f})")
      loc_3 = self.PnwComboLayer.RenderMultiPoleMotionForMa(loc_1, self.PnwVPAvel, self.PnwRotPole, currentMa)
      self.parent.geoWhiteboard.draw_target(loc_3.long, loc_3.lat, f"{currentMa} Ma YHS ({loc_3.long:0.3f}, {loc_3.lat:0.3f})")

        
    return loc_3