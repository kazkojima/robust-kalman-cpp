/*
  Robust Kalman filter eigen/cppoptlib/c++ implementation.
  Based on the python implementation by Milos Stojanovic (github: milsto).
  Also uses CppNumericalSolver by Parick Wieschollek (github: PatWie).
*/

#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <cppoptlib/meta.h>
#include <cppoptlib/problem.h>
#include <cppoptlib/solver/neldermeadsolver.h>
//#include <cppoptlib/solver/bfgssolver.h>
//#include <cppoptlib/solver/gradientdescentsolver.h>

template <int N, int M>
class RKF
{
public:
  RKF (Eigen::Matrix<float,N,N> Fm, // State transition matrix
       Eigen::Matrix<float,N,Eigen::Dynamic> Bm, // Input transition matrix
       Eigen::Matrix<float,M,N> Hm, // Observation matrix
       Eigen::Matrix<float,N,1> x0, // Initial state vector
       Eigen::Matrix<float,N,N> P0, // Initial state covariance matrix
       Eigen::Matrix<float,N,N> Q0, // (Initial) state noise covariance
       Eigen::Matrix<float,M,M> R0  // (Initial) observation noise covariance
       );
  void time_update (Eigen::VectorXf inputs);
  void measurement_update(Eigen::Matrix<float,M,1> meas);
  Eigen::Matrix<float,N,1> current_estimate (void);

private:
  Eigen::Matrix<float,N,N> F; // State transition matrix
  Eigen::Matrix<float,N,Eigen::Dynamic> B; // Input transition matrix
  Eigen::Matrix<float,M,N> H; // Observation matrixEigen::VectorXf x; // state vector
  Eigen::VectorXf x; // state vector
  Eigen::Matrix<float,N,N> P; // state covariance matrix
  Eigen::Matrix<float,N,N> Q; // state noise covariance
  Eigen::Matrix<float,M,M> R; // observation noise covariance
  Eigen::VectorXf inovation; // Residual or inovation
  Eigen::Matrix<float,M,M> Pinov; // inovation covariance matrix
  Eigen::Matrix<float,N,M> K; // Kalman gain

  bool noinput = false;

  static float robust_score (float z, float delta = 1.5f);

  class estimate_criterion : public cppoptlib::Problem<float>
  {
  public:
    estimate_criterion (Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> X,
			Eigen::Matrix<float,Eigen::Dynamic,1> Y) :
      EX(X), EY(Y)
    {
      //cout << "X:" << X << endl;
      //cout << "Y:" << Y << endl;
    }
    float value(const Eigen::Matrix<float,Eigen::Dynamic,1> &x)
    {
      float crit = 0.0f;

      //cout << "&x:" << x << endl;
      Eigen::Matrix<float,Eigen::Dynamic,1> nm = EY - EX*x;
      //cout << "nm:" << nm << endl;
      for (int i = 0; i < nm.rows(); i++)
	crit += robust_score (nm(i));
      //cout << "->" << crit << endl;
      return crit;
    }
  private:
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> EX;
    Eigen::Matrix<float,Eigen::Dynamic,1> EY;
  };
};

template <int N, int M>
float
RKF<N, M>::robust_score (float z, float delta)
{
  if (fabsf (z) >= delta)
    return delta*fabsf (z) - powf(delta,2)/2.0f;
  else
    return powf (z, 2)/2.0f;
}

template <int N, int M>
RKF<N, M>::RKF (Eigen::Matrix<float,N,N> Fm,
		Eigen::Matrix<float,N,Eigen::Dynamic> Bm,
		Eigen::Matrix<float,M,N> Hm,
		Eigen::Matrix<float,N,1> x0,
		Eigen::Matrix<float,N,N> P0,
		Eigen::Matrix<float,N,N> Q0,
		Eigen::Matrix<float,M,M> R0)
{
  F = Fm;
  B = Bm;
  H = Hm;
  x = x0;
  P = P0;
  Q = Q0;
  R = R0;
  if (B.squaredNorm() == 0)
    noinput = true;
}

template <int N, int M>
void
RKF<N, M>::time_update (Eigen::VectorXf inputs)
{
  if (noinput)
    x = F*x;
  else
    x = F*x + B*inputs;
  P = F*P*(F.transpose()) + Q;
  //cout << "x:" << x << endl;
  //cout << "P:" << P << endl;
}

template <int N, int M>
void
RKF<N, M>::measurement_update (Eigen::Matrix<float,M,1> meas)
{
  inovation = meas - H*x;
  Pinov = H*P*(H.transpose()) + R;
  K = (P*(H.transpose()))*Pinov.inverse();
  //cout << "K:" << K << endl;
  Eigen::MatrixXf epsilon_covariance(N+M, N+M);
  epsilon_covariance <<
    P, Eigen::MatrixXf::Zero(N, M),
    Eigen::MatrixXf::Zero(M, N), R;
  // LL^T Cholesky decomposition
  Eigen::LLT<Eigen::MatrixXf> llt;
  llt.compute(epsilon_covariance);
  Eigen::MatrixXf L = llt.matrixL();
  Eigen::MatrixXf IH(N+M, N);
  IH <<
    Eigen::MatrixXf::Identity(N, N),
    H;
  //cout << "L:" << (L.inverse())*IH << endl;
  Eigen::VectorXf ex(N+M);
  ex <<
    x,
    meas;

  // Nelder Mead solver
  estimate_criterion ecf((L.inverse())*IH, (L.inverse())*ex);
  //GradientDescentSolver<estimate_criterion> solver;
  //BfgsSolver<estimate_criterion> solver;
  cppoptlib::NelderMeadSolver<estimate_criterion> solver;
  solver.minimize(ecf, x);
  //cout << "x:" << x << endl;

  P = P - K*H*P;
  //cout << "P:" << P << endl;
}

template <int N, int M>
Eigen::Matrix<float,N,1>
RKF<N, M>::current_estimate (void)
{
  return x;
}
