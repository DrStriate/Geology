
# Gauss-Newton 2d solver for translation, rotation and scale in 2D 
# Below is c++ code from 2D LSDSLAM project using Displacements (extracted from Displacements\Core\GaussNewton2D.cpp)
'''
  //printf("GN Regress\n");
  float4 X { 0.0f, 0.0f, 0.0f, 0.0f };
  const int N = dArray.size();
  if (N >= 4) // need at least 4 samples to regress
  {
    const float2 c{ w / 2.0f, h / 2.0f };

    // pass2 - process data
    Eigen::MatrixXf J(2 * N, 4);
    Eigen::VectorXf R(2 * N);
    Eigen::VectorXf wt(2 * N);

    int jIdx = 0;
    for (int n = 0; n < N; n++)
    {
      const float iu = pArray[n].x;
      const float iv = pArray[n].y;

      const float u = iu - c.x;
      const float v = iv - c.y;

      const float uHat = dArray[n].x;
      const float vHat = dArray[n].y;
      const float dMag = dArray[n].z;
      const float s0 = dArray[n].w;

      float wu = wArray[n].x;
      float wv = wArray[n].y;
      float w = wArray[n].z;

      float ru = uHat * dMag; // Dx
      float rv = vHat * dMag; // Dy

      // Calculate Jacobian elements for Dx
      J(jIdx, 0) = -1.0f; //  dru/dtx
      J(jIdx, 1) = 0.0f;  //  dru/dty
      J(jIdx, 2) = u;     //  dru/ds
      J(jIdx, 3) = v;     //  dru/dtheta

      R(jIdx) = ru;  // Dx
      wt(jIdx) = wu; // wu;
      jIdx++;

      // Calculate Jacobian elements for Dy
      J(jIdx, 0) = 0.0f;  //  drv/dtx
      J(jIdx, 1) = -1.0f; //  drv/dty
      J(jIdx, 2) = v;     //  drv/ds
      J(jIdx, 3) = -u;    //  drv/dtheta

      R(jIdx) = rv;  // Dy
      wt(jIdx) = wv; // wv;
      jIdx++;
    }

    // std::cout << "J: \n" << J << std::endl;
    // std::cout << "R: \n" << R << std::endl;
    Eigen::DiagonalMatrix<float, Eigen::Dynamic> W(2 * N);
    W.diagonal() = wt;
    Eigen::MatrixXf A = J.transpose() * W * J;
    // std::cout << "A: \n" << A << std::endl;
    Eigen::MatrixXf b = -J.transpose() * W * R;
    // std::cout << "b: \n" << b << std::endl;

    // Residual rms
    if (R2)
      *R2 = sqrt(R.transpose() * W * R);
        
    // Eigen::Vector4f xVector = A.inverse() * b;
    // Eigen::Vector4f xVector = A.colPivHouseholderQr().solve(b).col(0);
    // Eigen::Vector4f xVector = (A.transpose() * A).ldlt().solve(A.transpose() * b);
    Eigen::Vector4f xVector = A.colPivHouseholderQr().solve(b);
    // std::cout << "X: " << xVector.transpose();

    X = {xVector(0),
         xVector(1),
         xVector(2),
         xVector(3)};
  }
  return X;
};
'''
import numpy as np

# lats and longs should be normalized relatice to center of rotation for best results
def solve_gauss_newton_2D_transform(lats, longs, v_n, v_e, normalize = True): # degrees and mm/Y units 
  N = len(lats)
  x = {'t_x' : 0, 't_y': 0, 's' : 0, 'r' : 0}
  if N < 4: # need at least 4 points to solve
    return x
  
  c = [0.0, 0.0]
  if normalize:
    c[0] = np.mean(longs)
    c[1] = np.mean(lats)
  
  j = np.zeros((2 * N, 4))
  r = np.zeros(2 * N)

  j_idx = 0
  for i in range(N):
    u = longs[i] - c[0]
    v = lats[i] - c[1]

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

  # Solve the linear system A * omega = B using standard least squares
  # This solves the normal equation: omega = (A^T * A)^(-1) * A^T * B
  x, residuals, rank, s = np.linalg.lstsq(j, r, rcond=None)
  return  {'t_x' : x[0], 't_y': x[1], 's' : x[2], 'r' : x[3]}