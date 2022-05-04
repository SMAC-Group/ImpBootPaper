// ------------------
// Misc
// ------------------
#include "common.h"

// ------------------
// Sampling
// ------------------
Eigen::VectorXi sample_int(unsigned int n, const unsigned int& seed);

// ------------------
// Covariance estimates
// ------------------
Eigen::MatrixXd covariance(Eigen::MatrixXd& boot);
double median(Eigen::VectorXd x);
double mad(Eigen::VectorXd x);

// ------------------
// Link function
// ------------------
inline double logit(double x){return std::log(x / (0.1e1 - x));}
inline double logistic(double x){return 0.1e1 / (0.1e1 + std::exp(-x));}
inline double loglink(double x){return log(x);}
inline double d_loglink(double x){return 0.1e1 / x;}
inline double invloglink(double x){return exp(x);}
inline double d_invloglink(double x){return exp(x);}
