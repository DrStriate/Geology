
# Gauss-Newton 2d solver for translation, rotation and scale in 2D 
import numpy as np

# lats and longs should be normalized relatice to center of rotation for best results
def solve_gauss_newton_2D_transform(sample_e, sanmple_n, v_e, v_n, normalize = True): # meters and mm/Y units 
  N = len(sanmple_n)
  x = {'t_x' : 0, 't_y': 0, 's' : 0, 'r' : 0}
  if N < 4: # need at least 4 points to solve
    return x
  
  c = [0.0, 0.0]
  if normalize:
    c[0] = np.mean(sample_e)
    c[1] = np.mean(sanmple_n)
  
  j = np.zeros((2 * N, 4))
  r = np.zeros(2 * N)

  j_idx = 0
  for i in range(N):
    u = sample_e[i] - c[0]
    v = sanmple_n[i] - c[1]

    # Calculate Jacobian elements for Dx
    j[j_idx, 0] = 1.0
    j[j_idx, 1] = 0.0
    j[j_idx, 2] = u
    j[j_idx, 3] = v

    r[j_idx] = v_e[i]
    j_idx += 1

    j[j_idx, 0] = 0.0
    j[j_idx, 1] = 1.0
    j[j_idx, 2] = v
    j[j_idx, 3] = -u

    r[j_idx] = v_n[i]
    j_idx += 1

  x, residuals, rank, s = np.linalg.lstsq(j, r, rcond=None)
  return  {'t_x' : x[0], 't_y': x[1], 's' : x[2], 'r' : x[3]}