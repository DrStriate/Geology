#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <filesystem>
#include "../../eigen-3.4.0/Eigen/Dense"

#define sqr(x) ((x) * (x))

const struct mapBounds
{
  float minLat = 41.0;
  float maxLat = 50.0;
  float minLon = 236.0;
  float maxLon = 250.0;
} gpsBounds;

struct GPS_VData_Point
{
  float lon;
  float lat;
  float Ve;
  float Vn;
  float Se;
  float Sn;
  float Ren;
} dataPoint;

bool readDataFile(const std::string &filename, std::vector<GPS_VData_Point> &gpsData)
{
  std::ifstream file(filename);

  if (!file.is_open())
  {
    std::cerr << "Path: " << std::filesystem::current_path() << std::endl;
    std::cerr << "Error: Could not open file " << filename << std::endl;
    return false;
  }

  gpsData.resize(0);
  std::string line;
  // Read the file line by line
  while (std::getline(file, line))
  {
    if (line[0] != '/' && line[1] != '/')
    {
      std::istringstream iss(line);
      std::vector<float> entries;
      float entry;

      // Loop to extract numbers until the end of the stream
      while (iss >> entry)
        entries.push_back(entry);

      if (entries.size() == 7)
      {
        GPS_VData_Point dataPoint{entries[0], entries[1], entries[2], entries[3], entries[4], entries[5], entries[6]};
        if (dataPoint.lat > gpsBounds.minLat &&
            dataPoint.lat < gpsBounds.maxLat &&
            dataPoint.lon > gpsBounds.minLon &&
            dataPoint.lon < gpsBounds.maxLon)
          gpsData.push_back(dataPoint);
      }
    }
  }
  std::cout << "Loaded " << gpsData.size() << " points\n";
  file.close();
  return true;
}

bool getTransform12(
    std::vector<GPS_VData_Point> &pArray,
    Eigen::Vector4f &xVector,
    float *R2)
{
  xVector.setZero();
  const int N = pArray.size();
  if (N < 4) // need at least 4 samples to regress
    return false;

  // Starting center point estimation
  float cx = (gpsBounds.maxLon + gpsBounds.minLon) / 2.0f;
  float cy = (gpsBounds.maxLat + gpsBounds.minLat) / 2.0f;

  // pass2 - process data
  Eigen::MatrixXf J(2 * N, 4);
  Eigen::VectorXf R(2 * N);
  Eigen::VectorXf wt(2 * N);

  int jIdx = 0;
  for (int n = 0; n < N; n++)
  {
    const float u = pArray[n].lon - cx;
    const float v = pArray[n].lat - cy;

    float wu = 1.0f / sqr(pArray[n].Se);
    float wv = 1.0f / sqr(pArray[n].Sn);
    float ru = pArray[n].Ve;
    float rv = pArray[n].Vn;

    // Calculate Jacobian elements for Dx
    J(jIdx, 0) = 1.0f; //  dru/dtx
    J(jIdx, 1) = 0.0f; //  dru/dty
    J(jIdx, 2) = u;    //  dru/ds
    J(jIdx, 3) = v;    //  dru/dtheta

    R(jIdx) = ru;  // Dx
    wt(jIdx) = wu; // wu;
    jIdx++;

    // Calculate Jacobian elements for Dy
    J(jIdx, 0) = 0.0f; //  drv/dtx
    J(jIdx, 1) = 1.0f; //  drv/dty
    J(jIdx, 2) = v;    //  drv/ds
    J(jIdx, 3) = u;   //  drv/dtheta

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

  // xVector = A.inverse() * b;
  // xVector = A.colPivHouseholderQr().solve(b).col(0);
  // xVector = (A.transpose() * A).ldlt().solve(A.transpose() * b);
  xVector = A.colPivHouseholderQr().solve(b);

  xVector[0] += cx;
  xVector[1] += cy;

  return true;
};

int main()
{
  std::string gpsDataFileName = "./data/nshm2023_wus_v1.txt";
  std::vector<GPS_VData_Point> gpsData;
  if (readDataFile(gpsDataFileName, gpsData))
  {
    Eigen::Vector4f xVector;
    float R2;
    if (getTransform12(gpsData, xVector, &R2))
      std::cout << xVector.transpose() << std::endl;
  }
  return 0;
};