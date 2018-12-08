/* A tiny test */

#include <iostream>
using namespace std;

#include "rkf.h"

int main()
{
  Eigen::Matrix2f F, P0, Q0;
  Eigen::Matrix<float,1,2> H;
  Eigen::Vector2f x0, x0_kalman;
  Eigen::Matrix<float,1,1> R0;
  F << 1, 0.01,
    0, 1;
  H << 1, 0;
  x0 << 0.01, 0.01;
  P0 << 0.001, 0.001, 0.001, 0.001;
  Q0 << 2.4999997e-07, 4.9999999e-05, 4.9999999e-05, 9.9999998e-03;
  R0 << 0.01;
  x0_kalman = Eigen::Vector2f::Zero(2);

  RKF<2,1> rkf(F, Eigen::MatrixXf::Zero(2,2), H, x0_kalman, P0, Q0, R0);
  Eigen::Matrix<float,1,1> meas;
  rkf.time_update (Eigen::Vector2f::Zero(2));
  meas << 0.1098551;
  rkf.measurement_update(meas);
  cout << "estimate:" <<rkf.current_estimate() << endl;
  rkf.time_update (Eigen::Vector2f::Zero(2));
  meas << 0.3226539;
  rkf.measurement_update(meas);
  cout << "estimate:" << rkf.current_estimate() << endl;
}
