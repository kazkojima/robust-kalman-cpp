# robust-kalman-cpp

Robust Kalman filter eigen/cppoptlib/c++ implementation based on the [python implementation](https://github.com/milsto/robust-kalman) by Milos Stojanovic (github: milsto).  Also uses CppNumericalSolver by Parick Wieschollek (github: PatWie).  See milsto's nice explanation in his github for detail.

Highly experimental. This header file only implementation gives
```
template <int N, int M> class RKF
```
where N is number of states and M is number of outputs.

When compiling, Eigen and [CppNumericalSolvers](https://github.com/PatWie/CppNumericalSolvers) headers have to be seen by your compiler.  For example, if Eigen is installed as /usr/include/eigen3/Eigen and CppNumericalSolvers is installed as /usr/local/include/cppoptlib, the optional compiler option -I/usr/include/eigen3 -I/usr/local/include/ will be needed.
