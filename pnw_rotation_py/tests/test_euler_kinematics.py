import numpy as np
import euler_pole.euler_pole_regression as epr
import euler_pole.euler_kinematics as ek

def test_euler_kinematics():
  #test setup
  euler_point = {"lat" : 45.0,  "long" : -90}
  omega = 1.23 # degrees / Ma (or mm/yr)
  sample_bearings  = [45.0, 135.0, 225.0, 315.0]
  sample_dist = 50000 # m

  sample_lats = []
  sample_lons = []
  sample_v_east = []
  sample_v_north = [] # mm/ yr

  Omega = {"omega": np.radians(omega), "phi": np.radians(euler_point['lat']), "lamb": np.radians(euler_point['long'])}
  for i in range(len(sample_bearings)):
    sample = ek.create_sample(euler_point['long'], euler_point['lat'], sample_bearings[i], sample_dist)
    print(f"{i}: sample: {sample}")
    sample_lats.append(sample['lat'])
    sample_lons.append(sample['lon'])
    p = {"phi": np.radians(sample_lats[i]), "lamb": np.radians(sample_lons[i])}
    v = ek.calculate_v_from_Eigen_pole(Omega, p); 
    sample_v_east.append(v['v_e'])
    sample_v_north.append(v['v_n'])
    print(f"{i}: v: {v}")

  pole_result5 = epr.fit_euler_pole_linear(sample_lats, sample_lons, sample_v_east, sample_v_north)
  ek.print_result ("fit_euler_pole_linear", pole_result5)

# test_euler_kinematics()