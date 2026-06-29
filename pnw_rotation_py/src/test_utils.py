def create_simple_sample_quad(euler_pole, sample_bearings, sample_dist):
  sample_e = []
  sample_n = []
  sample_v_east = []
  sample_v_north = [] # mm/ yr

  Omega = {"omega": euler_pole['omega'], "phi": np.radians(euler_pole['lat']), "lamb": np.radians(euler_pole['long'])}
  # print("")
  for i in range(len(sample_bearings)):
    sample = ek.create_sample(euler_pole['long'], euler_pole['lat'], sample_bearings[i], sample_dist)
    p = {"phi": np.radians(sample['lat']), "lamb": np.radians(sample['lon'])}
    v = ek.calculate_v_from_Eigen_pole(Omega, p, Omega['omega']);
    sample_e.append(sample['lon'])
    sample_n.append(sample['lat'])
    sample_v_east.append(v['v_e'])
    sample_v_north.append(v['v_n'])
    #print(f"{i}: sample['lat']: {sample['lat']:.3f}, sample['lon']: {sample['lon']:.3f}, v_e: {v['v_e']:.2f}  v_n: {v['v_n']:.2f}")
  return sample_e, sample_n, sample_v_east, sample_v_north
  
def create_random_sample_ring(euler_pole, count, 
                              max_dist, 
                              test_omega, 
                              crop = 1.0, 
                              source_poll = None):
  rng = np.random.default_rng(seed=42)
  rands = rng.random(size=(count, 2))

  sample_n = []
  sample_e = []
  sample_v_east = []
  sample_v_north = [] # mm/ yr

  Omega = {"omega": euler_pole['omega'], "phi": np.radians(euler_pole['lat']), "lamb": np.radians(euler_pole['long'])}
  if source_poll:
    Omega_source = {"omega": source_poll['omega'], "phi": np.radians(source_poll['lat']), "lamb": np.radians(source_poll['long'])}
  else:
    Omega_source = Omega

  max_long =  ek.create_sample(euler_pole['long'], euler_pole['lat'], 90.0, max_dist)['lon']
  min_long =  ek.create_sample(euler_pole['long'], euler_pole['lat'], 270.0, max_dist)['lon']
  crop_long = min_long + (max_long - min_long) * crop
  cropped_samples = 0;
  for i in range(len(rands)):
    sample = ek.create_sample(euler_pole['long'], euler_pole['lat'], 360.0 * rands[i][0], max_dist * rands[i][1])
    p = {"phi": np.radians(sample['lat']), "lamb": np.radians(sample['lon'])}
    v = ek.calculate_v_from_Eigen_pole(Omega_source, p, test_omega); 

    if sample['lon'] < crop_long:
      sample_e.append(sample['lon'])
      sample_n.append(sample['lat'])
      sample_v_east.append(v['v_e'])
      sample_v_north.append(v['v_n'])

    # print(f"{i}: sample['lat']: {sample['lat']:.3f}, sample['lon']: {sample['lon']:.3f}, v_e: {v['v_e']:.2f}  v_n: {v['v_n']:.2f}")
    cropped_samples += 1
  #print(f"samples = {cropped_samples} out of {count}")
  
  return sample_n, sample_e, sample_v_east, sample_v_north
