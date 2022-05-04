// ------------------
// Misc
// ------------------
#include "common.h"

// ------------------
// Sampling
// ------------------

Eigen::VectorXi sample_int(
    unsigned int n,
    const unsigned int& seed
){
  Eigen::VectorXi x(n);

  std::mt19937_64 engine(seed);  // Mersenne twister random number engine
  std::uniform_int_distribution<int> distribution(1,100000000);

  for(unsigned int i(0); i<n; ++i){
    x(i) = distribution(engine);
  }

  return x;
}

// ------------------
// Covariance estimates
// ------------------
Eigen::MatrixXd covariance(
    Eigen::MatrixXd& boot
){
  unsigned int n = boot.rows() - 1;
  Eigen::MatrixXd centered = boot.rowwise() - boot.colwise().mean();
  Eigen::MatrixXd cov = (centered.adjoint() * centered) / n;
  return cov;
}

double median(Eigen::VectorXd x){
  std::sort(x.data(), x.data()+x.size());
  return x.size() % 2 == 0 ?
  x.segment((x.size()-2)/2, 2).mean() :
    x(x.size()/2);
}

double mad(Eigen::VectorXd x){
  double center = median(x);
  Eigen::VectorXd y(x.size());
  for(unsigned int i(0); i<x.size(); ++i){
    y(i) = std::abs(x(i) - center);
  }
  return 1.4826 * median(y);
}
