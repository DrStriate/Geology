import geo_helper as gh

# Gauss-Newton 2d solver for translation, rotation and scale in 2D 
import numpy as np
def solve_gauss_newton_2D_transform_geo(sample_long, sample_lat, v_e, v_n, euler_pole, normalize = True): # meters and mm/Y units 
  sample_e, sample_n = gh.getSamplePoints(sample_long, sample_lat, euler_pole)
  return solve_gauss_newton_2D_transform(sample_e, sample_n, v_e, v_n, normalize)

# lats and longs should be normalized relative to "center" of rotation for best results
def solve_gauss_newton_2D_transform(sample_e, sample_n, v_e, v_n, normalize = True, ): # meters and mm/Y units 
  N = len(v_e)
  x = {'t_x' : 0, 't_y': 0, 's' : 0, 'r' : 0}
  if N < 4: # need at least 4 points to solve
    return x
  
  c = [0.0, 0.0]
  if normalize:
    c[0] = np.mean(sample_e)
    c[1] = np.mean(sample_n)
  
  J = np.zeros((2 * N, 4))
  R = np.zeros(2 * N)

  j_idx = 0
  for i in range(N):
    u = sample_e[i] - c[0]
    v = sample_n[i] - c[1]

    # Calculate Jacobian elements for Dx
    J[j_idx, 0] = 1.0
    J[j_idx, 1] = 0.0
    J[j_idx, 2] = u
    J[j_idx, 3] = v

    R[j_idx] = v_e[i]
    j_idx += 1

    J[j_idx, 0] = 0.0
    J[j_idx, 1] = 1.0
    J[j_idx, 2] = v
    J[j_idx, 3] = -u

    R[j_idx] = v_n[i]
    j_idx += 1

    # JT = J.T
    # JTJ = JT.dot(J)
    # JTR = JT.dot(R)

  x, residuals, rank, s = np.linalg.lstsq(J, R, rcond=None)
  return  {'t_x' : x[0], 't_y': x[1], 's' : x[2], 'r' : x[3]}

# Gauss-Newton 2d *weighted* solver for translation, rotation and scale in 2D 
import numpy as np
def solve_gauss_newton_2D_transform_geo_wtd(sample_long, sample_lat, v_e, v_n, w_e, w_n, center_ploc, normalize = True): # meters and mm/Y units 
  sample_e, sample_n = gh.getSamplePoints(sample_long, sample_lat, center_ploc)
  return solve_gauss_newton_2D_transform_wtd(sample_e, sample_n, v_e, v_n, w_e, w_n, normalize)

# lats and longs should be normalized relative to "center" of rotation for best results
def solve_gauss_newton_2D_transform_wtd(sample_e, sample_n, v_e, v_n, w_e, w_n, normalize = True, ): # meters and mm/Y units 
  N = len(v_e)
  x = {'t_x' : 0, 't_y': 0, 's' : 0, 'r' : 0}
  if N < 4: # need at least 4 points to solve
    return x
  
  c = [0.0, 0.0]
  if normalize:
    c[0] = np.mean(sample_e)
    c[1] = np.mean(sample_n)
  
  J = np.zeros((2 * N, 4))
  R = np.zeros(2 * N)

  j_idx = 0
  for i in range(N):
    u = sample_e[i] - c[0]
    v = sample_n[i] - c[1]

    # Current root weights for this point
    sw_e = 1.0 / w_e[i]
    sw_n = 1.0 / w_n[i]

    # Weighted Jacobian and Residual elements for Easting (Dx)
    J[j_idx, 0] = 1.0 * sw_e
    J[j_idx, 1] = 0.0 * sw_e
    J[j_idx, 2] = u   * sw_e
    J[j_idx, 3] = v   * sw_e
    
    R[j_idx] = v_e[i] * sw_e
    j_idx += 1
    
    # Weighted Jacobian and Residual elements for Northing (Dy)
    J[j_idx, 0] = 0.0 * sw_n
    J[j_idx, 1] = 1.0 * sw_n
    J[j_idx, 2] = v   * sw_n
    J[j_idx, 3] = -u  * sw_n
    
    R[j_idx] = v_n[i] * sw_n
    j_idx += 1
        
  # np.linalg.lstsq solves the system: (sqrt(W)*J)^T * (sqrt(W)*J) * x = (sqrt(W)*J)^T * (sqrt(W)*r)
  # Which simplifies exactly to: J^T * W * J * x = J^T * W * r
  x, residuals, rank, s = np.linalg.lstsq(J, R, rcond=None)
  
  return {'t_x' : x[0], 't_y': x[1], 's' : x[2], 'r' : np.degrees(x[3])}

def print_x(x):
  print(f"\nt_x:\t {x['t_x']:.5f}")
  print(f"t_y:\t {x['t_y']:.5f}")
  #print(f"s:  \t {x['s']:.5f}")
  print(f"r:  \t {x['r']:.5f}°")