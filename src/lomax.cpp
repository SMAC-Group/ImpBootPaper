// ------------------
// Description
// ------------------
// C++ implementation of the maximum likelihood estimator
// and the Switched Z-estimator for the Lomax distribution, also
// known as Pareto of type II. The Lomax distribution is parametrized
// by two positive parameters. We use log-transformation for the
// optimization procedures.

// ------------------
// Header
// ------------------
#include "common.h"

// ------------------
// Implementation
// ------------------
// ------------------
// Data Generating Processes
// ------------------
Eigen::VectorXd y_lomax(
    unsigned int n,
    double alpha,
    double eta,
    unsigned int seed
){
  Eigen::VectorXd y(n);
  double b,q,u;

  b = invloglink(alpha);
  q = invloglink(eta);

  std::mt19937_64 engine(seed);  // Mersenne twister random number engine
  std::uniform_real_distribution<double> unif(0.0, 1.0);

  for(unsigned int i(0);i<n;++i){
    u = unif(engine);
    y(i) = b * std::pow(u, -1.0 / q) - b;
  }

  return y;
}

Eigen::MatrixXd d_y_lomax(
    unsigned int n,
    double alpha,
    double eta,
    unsigned int seed
){
  Eigen::MatrixXd dy(n,2);
  double b,db,q,dq,u,t1,t2;

  b = invloglink(alpha);
  db = d_invloglink(alpha);
  q = invloglink(eta);
  dq = d_invloglink(eta);

  std::mt19937_64 engine(seed);  // Mersenne twister random number engine
  std::uniform_real_distribution<double> unif(0.0, 1.0);

  t2 = dq * b / q / q;

  for(unsigned int i(0);i<n;++i){
    u = unif(engine);
    t1 = std::pow(u, -1.0 / q);
    dy(i,0) = db * t1 - db;
    dy(i,1) = t2 * log(u) * t1;
  }

  return dy;
}

//' @title Random generation of Lomax distribution
//'
//' @param n sample size
//' @param b positive parameter
//' @param q positive parameter
//' @param seed integer representing the state fir random number generation
//' @export
// [[Rcpp::export]]
Eigen::VectorXd rlomax(
    unsigned int n,
    double b,
    double q,
    unsigned int seed
){
  return y_lomax(n,loglink(b),loglink(q),seed);
}

// ------------------
// Log-likelihood
// ------------------
double ll_lomax(
    const Eigen::Vector2d& theta,
    const Eigen::VectorXd& y
){
  unsigned int n = y.size();
  Eigen::VectorXd xi(n);
  double of,b,q,t0;

  b = invloglink(theta(0));
  q = invloglink(theta(1));
  t0 = q + 1.0;
  xi = (1.0 + y.array() / b).log();
  of = std::log(b) - std::log(q) + t0 * xi.mean();

  return of;
}

Eigen::Vector2d d_ll_lomax(
    const Eigen::Vector2d& theta,
    const Eigen::VectorXd& y
){
  unsigned int n = y.size();
  Eigen::VectorXd xi(n),xi1(n);
  Eigen::ArrayXd v(n);
  double b,db,q,dq,t0;
  Eigen::Vector2d grad;

  b = invloglink(theta(0));
  db = d_invloglink(theta(0));
  q = invloglink(theta(1));
  dq = d_invloglink(theta(1));
  t0 = q + 1.0;
  v = 1.0 + y.array() / b;
  xi = v.log();
  xi1 = 1.0 / v;

  grad(0) = db / b - db * t0 * y.dot(xi1) / b / b / n;
  grad(1) = -dq / q + dq * xi.mean();

  return grad;
}

Eigen::Vector2d psi_lomax(
    const Eigen::Vector2d& theta,
    const Eigen::Vector2d& pi,
    unsigned int n,
    unsigned int seed
){
  Eigen::VectorXd xi(n),xi1(n),y(n);
  Eigen::ArrayXd v(n);
  double b,db,q,dq,t0;
  Eigen::Vector2d grad;

  y = y_lomax(n,theta(0),theta(1),seed);

  b = invloglink(pi(0));
  db = d_invloglink(pi(0));
  q = invloglink(pi(1));
  dq = d_invloglink(pi(1));
  t0 = q + 1.0;
  v = 1.0 + y.array() / b;
  xi = v.log();
  xi1 = 1.0 / v;

  grad(0) = db / b - db * t0 * y.dot(xi1) / b / b / n;
  grad(1) = -dq / q + dq * xi.mean();

  return grad;
}

Eigen::Matrix2d d_psi_lomax(
    const Eigen::Vector2d& theta,
    const Eigen::Vector2d& pi,
    unsigned int n,
    unsigned int seed
){
  Eigen::VectorXd xi(n),xi1(n),y(n);
  Eigen::MatrixXd d_y(n,2);
  Eigen::ArrayXd v(n);
  double b,db,q,dq,t0;
  Eigen::Matrix2d jac;

  y = y_lomax(n,theta(0),theta(1),seed);
  d_y = d_y_lomax(n,theta(0),theta(1),seed);

  b = invloglink(pi(0));
  db = d_invloglink(pi(0));
  q = invloglink(pi(1));
  dq = d_invloglink(pi(1));
  t0 = q + 1.0;
  v = 1.0 + y.array() / b;
  xi1 = 1.0 / v;
  xi = y.array() / v / v / b;

  jac.row(0) = - db * t0 * (xi1 - xi).transpose() * d_y / b / b / n;
  jac.row(1) = dq * xi1.transpose() * d_y / b / n;

  return jac;
}

// ------------------
// Maximum likelihood
// ------------------

class mle_lomax: public Numer::MFuncGrad
{
private:
  const Eigen::VectorXd y;

public:
  mle_lomax(const Eigen::VectorXd y_) : y(y_) {}
  double f_grad(Numer::Constvec& theta, Numer::Refvec gr);
};

double mle_lomax::f_grad(
    Numer::Constvec& theta,
    Numer::Refvec gr
){
  const double of = ll_lomax(theta,y);
  gr = d_ll_lomax(theta,y);
  return of;
}

//' @title Maximum likelihood estimation of beta regression
//'
//' @param start initial values
//' @param y observations
//' @param maxit maximum number of iteration
//' @param eps_f tolerance
//' @param eps_g tolerance
//' @export
// [[Rcpp::export]]
Rcpp::List optim_mle_lomax(
    Eigen::VectorXd& start,
    Eigen::VectorXd& y,
    int maxit = 500,
    double eps_f = 1e-10,
    double eps_g = 1e-10
){
  double fopt;
  mle_lomax f(y);
  Eigen::Vector2d theta;
  theta(0) = std::log(start(0));
  theta(1) = std::log(start(1));
  int status = Numer::optim_lbfgs(f, theta, fopt, maxit, eps_f, eps_g);
  theta(0) = std::exp(theta(0));
  theta(1) = std::exp(theta(1));
  return Rcpp::List::create(
    Rcpp::Named("par") = theta,
    Rcpp::Named("value") = fopt,
    Rcpp::Named("conv") = status
  );
}

Rcpp::List of_mle_lomax(
    Eigen::VectorXd& theta,
    Eigen::VectorXd& y
){
  double fopt;
  Eigen::VectorXd gr(2);
  mle_lomax f(y);
  fopt = f.f_grad(theta,gr);
  return Rcpp::List::create(
    Rcpp::Named("f") = fopt,
    Rcpp::Named("grad") = gr
  );
}

// ------------------
// Switched Z-estimator
// ------------------

class swiz_lomax: public Numer::MFuncGrad
{
private:
  const Eigen::Vector2d pi;
  const unsigned int n;
  const unsigned int seed;

public:
  swiz_lomax(const Eigen::Vector2d pi_, const unsigned int n_, const unsigned int seed_) :
  pi(pi_), n(n_), seed(seed_) {}
  double f_grad(Numer::Constvec& theta, Numer::Refvec gr);
};

double swiz_lomax::f_grad(
    Numer::Constvec& theta,
    Numer::Refvec gr
){
  Eigen::Vector2d p;
  Eigen::Matrix2d dp;
  p = psi_lomax(theta,pi,n,seed);
  dp = d_psi_lomax(theta,pi,n,seed);
  const double of = p.squaredNorm() / 0.2e1;
  gr = dp.transpose() * p;
  return of;
}

Rcpp::List optim_swiz_lomax(
    Eigen::VectorXd& theta,
    Eigen::VectorXd& pi,
    unsigned int n,
    unsigned int seed
){
  double fopt;
  // Eigen::Vector4d gr;
  swiz_lomax f(pi,n,seed);
  // fopt = f.f_grad(theta,gr);

  int status = Numer::optim_lbfgs(f, theta, fopt);
  return Rcpp::List::create(
    // Rcpp::Named("f") = fopt,
    // Rcpp::Named("grad") = gr
    Rcpp::Named("par") = theta,
    Rcpp::Named("value") = fopt,
    Rcpp::Named("conv") = status
  );
}

Rcpp::List of_swiz2_lomax(
    Eigen::VectorXd& theta,
    Eigen::VectorXd& pi,
    unsigned int n,
    unsigned int seed
){
  Eigen::Vector2d p;
  Eigen::Matrix2d dp;
  p = psi_lomax(theta,pi,n,seed);
  dp = d_psi_lomax(theta,pi,n,seed);

  return Rcpp::List::create(
    Rcpp::Named("psi") = p,
    Rcpp::Named("jac") = dp
  );
}


Rcpp::List of_swiz_lomax(
    Eigen::VectorXd& theta,
    Eigen::VectorXd& pi,
    unsigned int n,
    unsigned int seed
){
  double fopt;
  Eigen::Vector2d gr;
  swiz_lomax f(pi,n,seed);
  fopt = f.f_grad(theta,gr);

  return Rcpp::List::create(
    Rcpp::Named("f") = fopt,
    Rcpp::Named("grad") = gr
  );
}

//' @title SwiZ distribution for Lomax distribution
//'
//' @param pi initial estimator
//' @param n sample size
//' @param B number of SwiZ estimates
//' @param seed integer representing the state fir random number generation
//' @param ncores number of cores (OpenMP parallelisation)
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd swiz_dist_lomax(
    Eigen::VectorXd& pi,
    unsigned int n,
    unsigned int B,
    unsigned int seed,
    unsigned int ncores
){
  Eigen::MatrixXd boot(B,2);

  #pragma omp parallel for num_threads(ncores)
  for(unsigned int i=0; i<B; ++i){
    Eigen::Vector2d theta;
    theta(0) = std::log(pi(0));
    theta(1) = std::log(pi(1));
    unsigned int se = seed + i;
    double fopt;
    swiz_lomax f(theta,n,se);
    int maxit = 500;
    double eps_f = 1e-10;
    double eps_g = 1e-10;
    Numer::optim_lbfgs(f, theta, fopt, maxit, eps_f, eps_g);
    boot(i,0) = std::exp(theta(0));
    boot(i,1) = std::exp(theta(1));
  }

  return boot;
}

// ------------------
// Parametric Bootstrap
// ------------------

//' @title Parametric bootstrap distribution for Lomax distribution
//'
//' @param start initial estimator
//' @param n sample size
//' @param B number of SwiZ estimates
//' @param seed integer representing the state fir random number generation
//' @param ncores number of cores (OpenMP parallelisation)
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd par_bootstrap_mle_lomax(
    Eigen::VectorXd& start,
    unsigned int n,
    unsigned int B,
    unsigned int seed,
    unsigned int ncores
){
  Eigen::MatrixXd boot(B,2);

  #pragma omp parallel for num_threads(ncores)
  for(unsigned int i=0; i<B; ++i){
    Eigen::Vector2d theta;
    theta(0) = std::log(start(0));
    theta(1) = std::log(start(1));
    unsigned int se = seed + i;
    // Eigen::VectorXd yt = y(n,start(0),start(1),se);
    Eigen::VectorXd y = rlomax(n,start(0),start(1),se);
    double fopt;
    mle_lomax f(y);
    int maxit = 300;
    double eps_f = 1e-8;
    double eps_g = 1e-7;
    Numer::optim_lbfgs(f, theta, fopt, maxit, eps_f, eps_g);
    boot(i,0) = std::exp(theta(0));
    boot(i,1) = std::exp(theta(1));
  }

  return boot;
}

//' @title Studentized parametric bootstrap distribution for beta regression
//'
//' @param theta initial estimator
//' @param boot matrix, parametric bootstrap
//' @param n sample size
//' @param B number of SwiZ estimates
//' @param seed integer representing the state fir random number generation
//' @param ncores number of cores (OpenMP parallelisation)
//' @param robust if true uses robust estimation of covariance
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd par_boott_lomax(
    Eigen::VectorXd& theta,
    Eigen::MatrixXd& boot,
    unsigned int n,
    unsigned int B,
    unsigned int seed,
    unsigned int ncores,
    bool robust=false
){
  unsigned int S = boot.rows();
  Eigen::MatrixXd boott(S,2);
  Eigen::VectorXi se(S);
  se = sample_int(S,seed);

  #pragma omp parallel for num_threads(ncores)
  for(unsigned int i=0; i<S; ++i){
    Eigen::MatrixXd new_boot(B,2);
    Eigen::VectorXd start = boot.row(i);
    new_boot = par_bootstrap_mle_lomax(start,n,B,se(i),1);
    Eigen::ArrayXd sd(2);
    if(robust){
      sd(0) = mad(new_boot.col(0));
      sd(1) = mad(new_boot.col(1));
    } else {
      Eigen::MatrixXd cov = covariance(new_boot);
      sd = cov.diagonal().array().sqrt();
    }
    boott.row(i) = (start - theta).array() / sd;
  }

  return boott;
}


// ------------------
// For BCa
// ------------------
//' @title BCa acceleration parameter for Lomax distribution
//'
//' @param start MLE
//' @param y observations
//' @param which binary, 0=alpha, 1=lambda
//' @export
// [[Rcpp::export]]
double acceleration_lomax(
    Eigen::VectorXd& start,
    Eigen::VectorXd& y,
    unsigned int which // 0 : alpha, 1 : lambda
){
  unsigned int n = y.size();
  int maxit = 300;
  double eps_f = 1e-8;
  double eps_g = 1e-7;
  Eigen::ArrayXd boot(n);
  Eigen::VectorXd yy(n-1);
  double t1(0);

  for(unsigned int i(0);i<n;++i){
    // I exploit the fact data has no order (iid)
    yy = y.tail(n-1);
    if(i != 0){
      yy(i-1) = y(0);
    }
    Eigen::Vector2d theta;
    theta(0) = std::log(start(0));
    theta(1) = std::log(start(1));
    double fopt;
    mle_lomax f(yy);
    Numer::optim_lbfgs(f, theta, fopt, maxit, eps_f, eps_g);
    boot(i) = std::exp(theta(which));
    t1 += boot(i) / n;
  }
  boot -= t1;
  double t2 = (boot * boot).sum();
  double t3 = (boot * boot * boot).sum();
  return t3 / 6.0 / std::pow(t2, 3/2);
}

// ------------------
// Old attempts (Fisher information, ...)
// ------------------
// double sigma2_alpha(double alpha, double lambda){
//   double t1,t2,t3;
//   t2 = lambda * lambda;
//   t1 = alpha * t2 + 2.0 * t2;
//   t3 = lambda + lambda * alpha;
//   return alpha / (t1 * (1.0 / alpha / t1 - 1.0 / t3 / t3));
// }
//
// double sigma2_lambda(double alpha, double lambda){
//   double t1,t2,t3;
//   t2 = lambda * lambda;
//   t1 = alpha * t2 + 2.0 * t2;
//   t3 = lambda + lambda * alpha;
//   return 1.0 / (alpha * alpha * (1.0 / (alpha * t1) - 1.0 / t3 / t3));
// }
//
// Eigen::MatrixXd par_boott_fisher(
//     Eigen::VectorXd& theta,
//     Eigen::MatrixXd& boot,
//     unsigned int n,
//     unsigned int B,
//     unsigned int seed,
//     unsigned int ncores
// ){
//   unsigned int S = boot.rows();
//   double sd;
//   Eigen::MatrixXd boott(S,2);
//   Eigen::VectorXi se(S);
//   se = sample_int(S,seed);
//
// #pragma omp parallel for num_threads(ncores)
//   for(unsigned int i=0; i<S; ++i){
//     Eigen::MatrixXd new_boot(B,2);
//     Eigen::VectorXd start = boot.row(i);
//     sd = sqrt(sigma2_alpha(start(0),start(1)) / n);
//     boott(i,0) = (start(0) - theta(0)) / sd;
//     sd = sqrt(sigma2_lambda(start(0),start(1)) / n);
//     boott(i,1) = (start(1) - theta(1)) / sd;
//   }
//
//   return boott;
// }
