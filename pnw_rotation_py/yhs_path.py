import numpy as np
from dataclasses import dataclass, asdict, replace
from .src import euler_pole_regression as epr
import json
import dacite 
import os

from .src.geo_helper import PAvel, PLoc, PVel, EulerPole, getPAvel, setGeod
from .src import test_utils as tu
from .src import euler_pole_regression as epr
from .src import gauss_newton as gn
from .src import euler_kinematics as ek
from . import path_layer as pl
from .na_plate_gplates import read_yhs_location as ry

SF = 1000.0 # conversion from units in km/ma and yrs to get meters (km/ma * yrs / SF = m)

YHS_lat = 44.43           # Yellowstone caldera yhs @ 0 Ma
YHS_long = -110.67

# NA_plate_speed = 23.0     # mm / yr (Current) = Adjusted to Owyhee=Humbolt cauldera ~14Ma
# NA_plate_azimuth = 241.0  # degrees 

@dataclass
class YhsPropertyBag:
  NaPlateDataName: str
  PnwVPAvel: PAvel
  PnwRotPole: EulerPole

class YhsPath:
  def __init__(self, parent):
    self.parent = parent
    self.path_layer_manager = pl.PathLayerManager()

    # Set up deault plate data
    self.setupNAPLateData("yhs_continuous_1ma_Muller2019.geojson")
    self.useGpsData = False
    self.pole_model = 2 # 1 is NA then Pole-V then Pole-R, 2 is NA - Translated-R - Pole-V
    self.yhs_loc = PLoc(YHS_long, YHS_lat)

    self.delta_ve = 0
    self.delta_vn = 0

  # Serialization (Object -> JSON) 
  def get_serialize_bag(self) -> str:
    # asdict recursively handles nesting automatically
    return json.dumps(asdict(self.getYhsPropertyBag()), indent=4)

  # Deserialization (JSON -> Object) 
  def deserialize_and_set_bag(self, json_str: str) -> YhsPropertyBag:
    data_dict = json.loads(json_str)
    
    propertyBag = dacite.from_dict(data_class=YhsPropertyBag, data=data_dict)
    self.NaPlateDataName = propertyBag.NaPlateDataName
    self.PnwVPAvel = propertyBag.PnwVPAvel
    self.PnwRotPole = propertyBag.PnwRotPole  
  
  def getYhsPropertyBag(self):
    propertyBag = YhsPropertyBag(self.NaPlateDataName, self.PnwVPAvel, self.PnwRotPole)
    return propertyBag
  
  def setYhsPropertyBag(self, propertyBag):
    self.NaPlateDataName = propertyBag.NaPlateDataName
    self.PnwVPAvel = propertyBag.PnwVPAvel
    self.PnwRotPole = propertyBag.PnwRotPole

  def setupNAPLateData(self, fileName):
    script_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "na_plate_gplates")
    na_file_path = os.path.join(script_dir, fileName)
    self.yhs_path_data = ry.load_data(na_file_path)
    self.NaPlateDataName = fileName

  def getPnwGpsRotPoleAndVelocity(self, sample_center, sample_radius):
      self.PnwRotPole, self.PnwVPAvel = epr.getPnwGpsRotPoleAndVelocity(sample_center, sample_radius)

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

  # Plot and label the NA Euler rotation pole
  def display_NA_pole_info(self):
    label_text1 = f"0 Ma YHS ({self.yhs_loc.long:.3f}, {self.yhs_loc.lat:.3f}), "
    # label_text2 = f"az: {self.NAPAvel.azimuth:.1f} deg, v: {(self.NAPAvel.vel):.1f} km/Ma"
    self.parent.geoWhiteboard.draw_target(self.yhs_loc.long, self.yhs_loc.lat, label_text1)
    self.parent.geoWhiteboard.draw_target(self.PnwRotPole.long, self.PnwRotPole.lat, 
                                 f"0 Ma pole ({self.PnwRotPole.long:0.3f}, {self.PnwRotPole.lat:0.3f})")
    
  def get_yhs_loc(self, currentMa, deltaMa, steps): # yrs is years
    self.checkLayersCreated()

    # get a local rot pole we can move as needed 
    runPnwRotPole = replace(self.PnwRotPole)
    #runPnwRotPole.print("runPnwRotPole: ") # these printouts useful for setting up GPlates
    # get the PnwVPole which is time invariant 
    self.PnwVPole = ek.getEulerPoleFromPlocAndPavel(runPnwRotPole.ploc(), self.PnwVPAvel)
    #self.PnwVPole.print("PnwVPole: ")

    # 1: Move yhs loc by NA speed scaled by ma from 0 Ma location (red line)
    # self.NAPole = ek.getEulerPoleFromPlocAndPavel(self.yhs_loc, self.NAPAvel)
    # loc_1 = self.NaPoleLayer.RenderPoleMotionForMa(self.yhs_loc, self.NAPole, -currentMa)
    loc_1 = self.NaPoleLayer.RenderYHSPoleMotionForMa(self.yhs_path_data, -currentMa)
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

  def displayGPSDataAndPoles(self, test, center, radius):
    setGeod(realWorld = True)    

    if test: # show synthetic data
      #pole for generating samples
      v_in = [0.767, 3.545] # v pavel from typical calibration
      v_pavel = PAvel.from_V(v_in) 
      pnwVPole = ek.getEulerPoleFromPlocAndPavel(center, v_pavel)

      sample_count = 400
      crop = 1.0 # no crop
      lat_list, long_list, ve_list, vn_list =\
        tu.create_random_sample_ring(pnwVPole, center, sample_count, radius, pnwVPole.omega, crop)
      s_e = None
      s_n = None
    else: # show GPS data
      lat_list, long_list, ve_list, vn_list, s_e, s_n =\
        tu.get_GPS_rotation_data(center.long, center.lat, radius * 1000.0)

    if len(lat_list) < 3:
      return False
    
    mod_ve_list = np.array(ve_list) - self.delta_ve
    mod_vn_list = np.array(vn_list) - self.delta_vn

    self.parent.rotDestLayer.dataProvider().truncate()
    self.parent.yhsRotFeatureList = []
    for i in range(len(lat_list)):
      feature = self.parent.rotData.createRotFeature(
          PLoc(long_list[i], lat_list[i]), PVel(mod_ve_list[i], mod_vn_list[i]))
      self.parent.yhsRotFeatureList.append(feature)        

    # get euler pole and gauss newton results and display
    self.PnwRotPole = epr.fit_euler_pole_linear(lat_list, long_list, mod_ve_list, mod_vn_list, s_e, s_n)
    self.PnwRotPole.print("self.PnwRotPole: ")

    offsets = gn.solve_gauss_newton_2D_transform_geo_wtd(long_list, lat_list, mod_ve_list, mod_vn_list, s_e, s_n, self.PnwRotPole )
    print(f"offsets: {offsets}")

    label_text1 = f"{self.PnwRotPole.long:.4f}, {self.PnwRotPole.lat:.4f}, {self.PnwRotPole.omega:.3f} deg, "
    label_text2 = f"e: {offsets[0]:.3f} km, n: {offsets[1]:.3f} km, {radius} km"
    self.parent.geoWhiteboard.draw_target(self.PnwRotPole.long, self.PnwRotPole.lat, label_text1 + label_text2)
    #print(label_text1 + label_text2)

    # if translation correction added, add a delta_V vector to the target to show that
    if self.delta_ve != 0.0 or self.delta_vn != 0.0:    
      feature = self.parent.rotData.createRotFeature(
        PLoc(self.PnwRotPole.long, self.PnwRotPole.lat), PVel(self.delta_ve + offsets[0], self.delta_vn + offsets[1])) 
      self.parent.yhsRotFeatureList.append(feature)
      self.PnwVPAvel  = getPAvel(self.delta_ve, self.delta_vn)
    else:
      self.PnwVPAvel = PAvel(0, 0)
  
    self.delta_ve = offsets[0]
    self.delta_vn = offsets[1]
    return True

  def clearGPSparams(self):
    self.delta_ve = 0.0
    self.delta_vn = 0.0  
