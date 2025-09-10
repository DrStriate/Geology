#include "gaussNewton2D.h"
#include <iostream>
#include "disparityFn.h"
#include "gradientFn.h"
#include "displacementFn.h"
#include <Eigen/Dense>

#define sqr(x) ((x)*(x))

float4 GaussNewton2D::getTransform12(
    Image<float4> &displacementImage1,
    Image<float4> &displacementImage2,
    float* R2)
{
  // pass 1 filter and get valid data count
  const int w = displacementImage2.Width();
  const int h = displacementImage2.Height();

  std::vector<float4> dArray;
  std::vector<float2> pArray;
  std::vector<float4> wArray;

  const float4 *dImage1 = displacementImage1.HData();
  const float4 *dImage2 = displacementImage2.HData();
  int idx = 0;
  for (int j = 0; j < displacementImage1.Height(); j++)
  {
    for (int i = 0; i < displacementImage1.Width(); i++)
    {
      const float4 d1 = dImage1[idx];
      const float4 d2 = dImage2[idx++];
      float4 D; // ux, uy, D
      float4 W; // wu wv
      float2 P; // pu pv
      if (DisparityFn::getDisparityCoeffs(d1, d2, i, j, D, W, P))
      {
        dArray.push_back(D);
        wArray.push_back(W);
        pArray.push_back(P);
      }
    }
  }

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

// Note: no rotation applied to gradients so pre-apply as needed for accuracy if image pre-transformed
float4 GaussNewton2D::getTransform12(
  Image<float>& laplacianImage1,
  Image<float2>& gradientImage1,
  Image<float>& laplacianImage2,
  Image<float2>& gradientImage2,
  float* R2,
  float sigma)
{
  // pass 1 filter and get valid data count
  const int w = laplacianImage1.Width();
  const int h = laplacianImage1.Height();

  // 1 - ref, 2 - frame
  Image<float4> displacementImage1;
  DisplacementFn::getDisplacement(displacementImage1, laplacianImage1, gradientImage1, sigma);
  Image<float4> displacementImage2;
  DisplacementFn::getDisplacement(displacementImage2, laplacianImage2, gradientImage2, sigma);

  std::vector<float4> dArray;
  std::vector<float2> pArray;
  std::vector<float4> wArray;

  const float4* dImage1 = displacementImage1.HData();
  const float4* dImage2 = displacementImage2.HData();
  int idx = 0;
  for (int j = 0; j < h; j++)
  {
    for (int i = 0; i < w; i++)
    {
      const float4 d1 = dImage1[idx];
      const float4 d2 = dImage2[idx++];
      float4 D; // ux, uy, D
      float4 W; // wu wv
      float2 P; // pu pv
      if (DisparityFn::getDisparityCoeffs(d1, d2, i, j, D, W, P))
      {
        dArray.push_back(D);
        wArray.push_back(W);
        pArray.push_back(P);
      }
    }
  }

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

  /*
  for (int j = 0; j < h; j++)
  {
    for (int i = 0; i < w; i++)
    {
      int idx = j * h + i;

      // 1 - ref, 2 - frame
      float rLaplacian = laplacianData1[idx];
      float rgu = gradientData1[idx].x;
      float rgv = gradientData1[idx].y;
      float4 d1 = DisplacementFn::getDisplacement(rLaplacian, rgu, rgv, sigma);


      float fLaplacian = laplacianData2[idx];
      float fgu = gradientData2[idx].x;
      float fgv = gradientData2[idx].y;
      float4 d2 = DisplacementFn::getDisplacement(fLaplacian, fgu, fgv, sigma);

      float4 D; // uu, uv, D, wt
      float4 W; // wu, wv
      float2 P; // pu, pv
      if (DisparityFn::getDisparityCoeffs(d1, d2, i, j, D, W, P))
      {
        dArray.push_back(D);
        wArray.push_back(W);
        pArray.push_back(P);
      }
    }
  }

  // printf("GN Regress\n");
  float4 X{ 0.0f, 0.0f, 0.0f, 0.0f};
  const int N = dArray.size();
  if (N >= 6) // need at least 6 samples to regress
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
      J(jIdx, 0) = -1.0f; //  dru/dTx
      J(jIdx, 1) = 0.0f;  //  dru/dTy
      J(jIdx, 2) = u;     //  dru/dTz
      J(jIdx, 3) = v;     //  dru/dRz

      R(jIdx) = ru;  // Dx
      wt(jIdx) = wu; // wu;
      jIdx++;

      // Calculate Jacobian elements for Dy
      J(jIdx, 0) = 0.0f;  //  drv/dTx
      J(jIdx, 1) = -1.0f; //  drv/dTy
      J(jIdx, 2) = v;     //  drv/dTz
      J(jIdx, 3) = -u;    //  drv/dRz

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
         xVector(3) };
  }
  return X;
  */
};

float4 GaussNewton2D::getTransformLSD(
  Image<float>& refImage,
  Image<float>& frameImage,
  float* R2,
  bool useGpu)
{
  // pass 1 filter and get valid data count
  const int w = refImage.Width();
  const int h = refImage.Height();
  std::shared_ptr<Image<float2>> refGradientImage = std::make_shared<Image<float2>>(w, h);
  GradientFn::Convolve(*refGradientImage, refImage, useGpu);

  std::shared_ptr<Image<float2>> frameGradientImage = std::make_shared<Image<float2>>(w, h);
  GradientFn::Convolve(*frameGradientImage, frameImage, useGpu);

  std::vector<float> dArray;
  std::vector<float2> pArray;
  std::vector<float4> wArray;

  const float* dRefImage = refImage.HData();
  const float* dFrameImage = frameImage.HData();
  const float2* dRefGradientImage = refGradientImage->HData();
  const float2* dFrameGradientImage = frameGradientImage->HData();
  int idx = 0;
  for (int j = 0; j < h; j++)
  {
    for (int i = 0; i < w; i++)
    {
      const float iRef = dRefImage[idx];
      const float iFrame = dFrameImage[idx];
      const float2 gRef = dRefGradientImage[idx];
      const float2 gFrame = dFrameGradientImage[idx++];
      if (sqrt(sqr(gRef.x) + sqr(gRef.y)) > GRAD_THRESH)
      {
        dArray.push_back(iFrame - iRef);
        pArray.push_back(float2{ (float)i, (float)j });
        wArray.push_back(float4{ gRef.x, gRef.y, gFrame.x, gFrame.y });
      }
    }
  }
  //printf("GN Regress\n");
  float4 X { 0.0f, 0.0f, 0.0f, 0.0f };
  const int N = dArray.size();
  if (N > 4) // need at least 4 samples to regress
  {
    const float2 c{ w / 2.0f, h / 2.0f };

    // pass2 - process data
    Eigen::MatrixXf J(N, 4);
    Eigen::VectorXf R(N);
    Eigen::VectorXf wt(N);

    for (int n = 0; n < N; n++)
    {
      const float iu = pArray[n].x;
      const float iv = pArray[n].y;

      const float u = iu - c.x;
      const float v = iv - c.y;

      float grx = wArray[n].x;
      float gry = wArray[n].y;
      float gfx = wArray[n].z;
      float gfy = wArray[n].w;

      float ivr = sqr(grx) + sqr(gry); // inverse variance est. of ref and frame
      float ivf = sqr(gfx) + sqr(gfy);
      float w = (ivf > 0? 1.0f / (1.0f / ivr + 1.0f / ivf) : 0.0f);

      float r = dArray[n];

      // Calculate Jacobian elements 
      J(n, 0) = -gfx; //  dru/dtx
      J(n, 1) = -gfy;  //  dru/dty
      J(n, 2) = u * gfx + v * gfy;     //  dru/ds
      J(n, 3) = v * gfx - u * gfy;     //  dru/dtheta

      R(n) = r;  // Dx
      wt(n) = w; // wu
    }
    //std::cout << "J: \n" << J << std::endl;  
    //std::cout << "R: \n" << R << std::endl;
    Eigen::DiagonalMatrix<float, Eigen::Dynamic> W(2 * N);
    W.diagonal() = wt;
    Eigen::MatrixXf A = J.transpose() * W * J;
    //std::cout << "A: \n" << A << std::endl;
    Eigen::MatrixXf b = -J.transpose() * W * R;
    //std::cout << "b: \n" << b << std::endl;

    // Residual rms
    if (R2)
      *R2 = sqrt(R.transpose() * W * R);

    // Eigen::Vector4f xVector = A.inverse() * b;
    // Eigen::Vector4f xVector = A.colPivHouseholderQr().solve(b).col(0);
    // Eigen::Vector4f xVector = (A.transpose() * A).ldlt().solve(A.transpose() * b);
    Eigen::Vector4f xVector = A.colPivHouseholderQr().solve(b);
    // std::cout << "X: " << xVector.transpose();

    X = { xVector(0),
         xVector(1),
         xVector(2),
         xVector(3) };
  }
  return X;
};
