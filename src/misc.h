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
);

// ------------------
// Covariance matrix
// ------------------
Eigen::MatrixXd covariance(
    Eigen::MatrixXd& boot
);

// ------------------
// Link function
// ------------------
inline double logit(double x){return std::log(x / (0.1e1 - x));}
inline double logistic(double x){return 0.1e1 / (0.1e1 + std::exp(-x));}



