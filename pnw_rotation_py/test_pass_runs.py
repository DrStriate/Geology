import test_utils as tu
from .plate_motion import PLoc, PDist
import euler_pole_regression as epr
import gauss_newton as gn

def run_GPS_test_pass(self):
  # get rot data
  diam = 550 # km
  lat_list, long_list, ve_list, vn_list =\
        tu.get_GPS_rotation_data(tu.OC_NA_Pole['lat'], tu.OC_NA_Pole['long'], diam * 1000)
  finish_test_setup(self, lat_list, long_list, ve_list, vn_list, diam)

def run_quad_test_pass(self):
  # get rot data
  euler_pole = tu.OC_NA_Pole #{"lat" : 45.0,  "long" : -90, "omega" : 1.23 }
  bearings  = [45.0, 135.0, 225.0, 315.0]
  diam = 50 # km
  long_list, lat_list, ve_list, vn_list = tu.create_simple_sample_quad(euler_pole, bearings, diam * 1000)
  finish_test_setup(self, lat_list, long_list, ve_list, vn_list, diam)

def run_rand_disk_test_pass(self):
  # get rot data
  euler_pole = tu.OC_NA_Pole #{"lat" : 45.0,  "long" : -90, "omega" : 1.23 }
  sample_count = 400
  diam = 400 # km
  crop = 1.0 # no crop
  test_omega = 1.23

  lat_list, long_list, ve_list, vn_list = \
    tu.create_random_sample_ring(euler_pole, sample_count, diam * 1000, test_omega, crop)
  finish_test_setup(self, lat_list, long_list, ve_list, vn_list, diam)

def run_cropped_disk_test_test_pass(self):
    #test setup
  euler_pole = tu.OC_NA_Pole #{"lat" : 45.0,  "long" : -90, "omega" : 1.23 }
  sample_count = 400
  diam = 400 # m
  test_omega = 1.23
  crop = 0.5 # 50% cropped out

  lat_list, long_list, ve_list, vn_list = \
    tu.create_random_sample_ring(
      euler_pole, 
      sample_count, 
      diam * 1000, 
      test_omega, 
      crop,
      None)
  finish_test_setup(self, lat_list, long_list, ve_list, vn_list, diam)
  
def finish_test_setup(self, lat_list, long_list, ve_list, vn_list, diam):
     # set rot data in Qgis
  for i in range(len(lat_list)):
    feature = self.rotData.createRotFeature(
        PLoc(long_list[i], lat_list[i]), PDist(ve_list[i], vn_list[i]), 0.001)
    self.yhsRotFeatureList.append(feature)        

  # get euler pole and gauss newton results and display
  pole = epr.fit_euler_pole_linear(lat_list, long_list, ve_list, vn_list)
  gn_out = gn.solve_gauss_newton_2D_transform_geo(long_list, lat_list, ve_list, vn_list, tu.OC_NA_Pole)
  e_km = gn_out['t_x'] / 1000
  n_km = gn_out['t_y'] / 1000
  label_text1 = f"{pole['long']:.4f}, {pole['lat']:.4f}, {pole['omega']:.2f} deg, "
  label_text2 = f"e: {e_km:.2f} km, n: {n_km:.2f} km, {diam} km"
  
  #show results in Qgis
  self.geoWhiteboard.draw_target(pole['long'], pole['lat'], label_text1 + label_text2)
