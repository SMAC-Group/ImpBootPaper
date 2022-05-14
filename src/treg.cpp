// ------------------
// Description
// ------------------
// C++ implementation of the maximum likelihood estimator
// and the Switched Z-estimator for the Student's t regression.
// The t regression is parametrized regression coefficients and
// two positive parameters: the variance and the degrees of freedom.
// We use std::log-transformation for the positive parameters for
// the optimization procedures.

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
Eigen::VectorXd y_treg(
    const Eigen::MatrixXd& x,
    const Eigen::VectorXd& beta,
    double alpha,
    double eta,
    unsigned int seed
){
  unsigned int n = x.rows();
  Eigen::VectorXd y(n), xb(n);
  double t0,nu;
  double u,v,w,c,r,e;

  t0 = sqrt(invloglink(alpha));
  nu = invloglink(eta);
  xb = x * beta;

  std::mt19937_64 engine(seed);  // Mersenne twister random number engine
  std::uniform_real_distribution<double> unif(-1.0, 1.0);

  for(unsigned int i(0);i<n;++i){
    // Bailey's polar algorithm
    do{
      u = unif(engine);
      v = unif(engine);
      w = u * u + v * v;
    } while (w > 1.0);
    c = u / std::sqrt(w);
    r = std::sqrt(nu * std::pow(w, -2.0 / nu) - nu);
    e = c * r;
    y(i) = xb(i) + t0 * e;
  }

  return y;
}

// derivative
Eigen::MatrixXd d_y_treg(
    const Eigen::MatrixXd& x,
    const Eigen::VectorXd& beta,
    double alpha,
    double eta,
    unsigned int seed
){
  unsigned int n = x.rows();
  unsigned int p = beta.size() + 2;
  Eigen::MatrixXd dy(n,p);
  double t0,t1,t2,t3,t4,t5,nu;
  double u,v,w,c,r,e;

  nu = invloglink(eta);
  t0 = sqrt(invloglink(alpha));
  t1 = sqrt(nu);
  t2 = d_invloglink(alpha) / t0 / 2.0;
  t3 = d_invloglink(eta) * t0 / t1 / 2.0;

  std::mt19937_64 engine(seed);  // Mersenne twister random number engine
  std::uniform_real_distribution<double> unif(-1.0, 1.0);

  dy.leftCols(p - 2) = x;

  for(unsigned int i(0);i<n;++i){
    // Bailey's polar algorithm
    do{
      u = unif(engine);
      v = unif(engine);
      w = u * u + v * v;
    } while (w > 1.0);
    c = u / std::sqrt(w);
    t4 = std::pow(w, -2.0 / nu);
    t5 = std::sqrt(t4 - 1.0);
    r = t1 * t5;
    e = c * r;
    dy(i, p-2) = t2 * e;
    dy(i, p-1) = t3 * c * t5 + t3 * c * 2.0 * std::log(w) * t4 / t5 / nu;
  }

  return dy;
}

Rcpp::List y_dy_treg(
    Eigen::MatrixXd& x,
    Eigen::VectorXd& beta,
    double alpha,
    double eta,
    unsigned int seed
){
  Eigen::VectorXd y = y_treg(x,beta,alpha,eta,seed);
  Eigen::MatrixXd d_y = d_y_treg(x,beta,alpha,eta,seed);

  return Rcpp::List::create(
    Rcpp::Named("y") = y,
    Rcpp::Named("dy") = d_y
  );
}

//' @title Random generation of t regression responses
//'
//' @param x a matrix of design (first column must be one to include an intercept)
//' @param beta a vector of regression coefficient
//' @param sig2 variance parameter
//' @param nu degrees of freedom
//' @param seed integer representing the state fir random number generation
//' @export
// [[Rcpp::export]]
Eigen::VectorXd rtreg(
    const Eigen::MatrixXd& x,
    const Eigen::VectorXd& beta,
    double sig2,
    double nu,
    unsigned int seed
){
  return y_treg(x,beta,loglink(sig2),loglink(nu),seed);
}

// ------------------
// std::log-likelihood
// ------------------
double ll_treg(
    const Eigen::VectorXd& theta,
    const Eigen::VectorXd& y,
    const Eigen::MatrixXd& x
){
  unsigned int n = y.size();
  unsigned int p = theta.size();
  Eigen::ArrayXd xi(n);
  double of,nu,sig2,t0;

  sig2 = invloglink(theta(p - 2));
  nu = invloglink(theta(p - 1));
  t0 = nu / 2.0 + 1.0 / 2.0;
  xi = y - x * theta.head(p - 2);
  of = std::lgamma(t0 - 0.5) - std::lgamma(t0)
    + std::log(nu * sig2) / 2.0 + t0 * (1.0 + xi * xi / nu / sig2).log().matrix().mean();

  return of;
}

Eigen::VectorXd d_ll_treg(
    const Eigen::VectorXd& theta,
    const Eigen::VectorXd& y,
    const Eigen::MatrixXd& x
){
  unsigned int n = y.size();
  unsigned int p = theta.size();
  Eigen::ArrayXd xi(n),xi1(n),xi2(n);
  double nu,sig2,t0,t1,t2,t3;
  Eigen::VectorXd grad(p);

  sig2 = invloglink(theta(p - 2));
  nu = invloglink(theta(p - 1));
  t0 = nu + 1.0;
  t1 = t0 / nu / sig2;
  t2 = d_invloglink(theta(p - 2)) / sig2 / 2.0;
  t3 = d_invloglink(theta(p - 1)) / 2.0;
  xi = y - x * theta.head(p - 2);
  xi2 = xi * xi;
  xi1 = 1.0 + xi2 / nu / sig2;

  grad.head(p - 2) = -t1 * x.transpose() * (xi / xi1).matrix() / n;
  grad(p - 2) = t2  - t2 * t1 * (xi2 / xi1).matrix().mean();
  grad(p - 1) = t3 * boost::math::digamma(nu / 2.0, my_policy()) - t3 * boost::math::digamma(t0 / 2.0, my_policy())
    + t3 / nu + t3 * xi1.log().matrix().mean() - t1 * t3 * (xi2 / xi1).matrix().mean() / nu;

  return grad;
}

Eigen::VectorXd psi_treg(
    const Eigen::VectorXd& theta,
    const Eigen::VectorXd& pi,
    const Eigen::MatrixXd& x,
    unsigned int seed
){
  unsigned int n = x.rows();
  unsigned int p = theta.size();
  Eigen::VectorXd y(n);
  y = y_treg(x,theta.head(p-2),theta(p-2),theta(p-1),seed);

  Eigen::ArrayXd xi(n),xi1(n),xi2(n);
  double nu,sig2,t0,t1,t2,t3;
  Eigen::VectorXd grad(p);

  sig2 = invloglink(pi(p - 2));
  nu = invloglink(pi(p - 1));
  t0 = nu + 1.0;
  t1 = t0 / nu / sig2;
  t2 = d_invloglink(pi(p - 2)) / sig2 / 2.0;
  t3 = d_invloglink(pi(p - 1)) / 2.0;
  xi = y - x * pi.head(p - 2);
  xi2 = xi * xi;
  xi1 = 1.0 + xi2 / nu / sig2;

  grad.head(p-2) = -t1 * x.transpose() * (xi / xi1).matrix() / n;
  grad(p-2) = t2  - t2 * t1 * (xi2 / xi1).matrix().mean();
  grad(p-1) = t3 * boost::math::digamma(nu / 2.0, my_policy()) - t3 * boost::math::digamma(t0 / 2.0, my_policy())
    + t3 / nu + t3 * xi1.log().matrix().mean() - t1 * t3 * (xi2 / xi1).matrix().mean() / nu;

  return grad;
}

Eigen::MatrixXd d_psi_treg(
    const Eigen::VectorXd& theta,
    const Eigen::VectorXd& pi,
    const Eigen::MatrixXd& x,
    unsigned int seed
){
  unsigned int n = x.rows();
  unsigned int k = x.cols();
  unsigned int p = theta.size();
  Eigen::VectorXd y(n);
  y = y_treg(x,theta.head(p-2),theta(p-2),theta(p-1),seed);
  Eigen::MatrixXd d_y(n,p);
  d_y = d_y_treg(x,theta.head(p-2),theta(p-2),theta(p-1),seed);

  Eigen::ArrayXd xi(n),xi1(n),xi2(n),v1(n),v2(n),v(n);
  double nu,sig2,t0,t1,t2,t4;
  Eigen::MatrixXd jac(p,p);

  sig2 = invloglink(pi(p - 2));
  nu = invloglink(pi(p - 1));
  t4 = nu * sig2;
  t0 = nu + 1.0;
  t1 = t0 / t4;
  t2 = d_invloglink(pi(p-2)) / sig2;
  xi = y - x * pi.head(p-2);
  xi2 = xi * xi;
  xi1 = 1.0 + xi2 / t4;
  v1 = xi / xi1;
  v2 = xi2 / xi1;

  for(unsigned int i(0);i<p;++i){
    v = d_y.col(i);
    jac.col(i).head(k) = -t1 * x.transpose() * (v / xi1 - 2.0 * v * v2 / xi1 / t4).matrix() / n;
    jac(p-2,i) = -t1 * t2 * (v * v1 - v * xi * v2 / xi1 / t4).matrix().mean();
    // jac(p-1,i) = ;
  }

  // jac.row(p-2) = -t1 * t2 * (v1 - xi * v2 / xi1 / t4).matrix().transpose() * d_yt / n / sig2;
  jac.row(p-1) = d_invloglink(pi(p - 1)) * (v1 - t0 * v1 / nu + t0 * v1 * v2 / t4 / nu).matrix().transpose() * d_y / n / t4;

  return jac;
}

// ------------------
// Maximum likelihood
// ------------------

class mle_treg: public Numer::MFuncGrad
{
private:
  const Eigen::VectorXd y;
  const Eigen::MatrixXd x;

public:
  mle_treg(const Eigen::VectorXd y_, const Eigen::MatrixXd x_) : y(y_), x(x_) {}
  double f_grad(Numer::Constvec& theta, Numer::Refvec gr);
};

double mle_treg::f_grad(
    Numer::Constvec& theta,
    Numer::Refvec gr
){
  const double of = ll_treg(theta,y,x);
  gr = d_ll_treg(theta,y,x);
  return of;
}

//' @title Maximum likelihood estimation of t regression
//'
//' @param start initial values (coef + log-variance + log-dof)
//' @param y responses
//' @param x matrix of design
//' @param maxit maximum number of iteration
//' @param eps_f tolerance
//' @param eps_g tolerance
// [[Rcpp::export]]
Rcpp::List optim_mle_treg(
    Eigen::VectorXd& start,
    Eigen::VectorXd& y,
    Eigen::MatrixXd& x,
    int maxit = 500,
    double eps_f = 1e-10,
    double eps_g = 1e-10
){
  unsigned int p = start.size();
  Eigen::VectorXd theta = start;
  theta(p-1) = std::log(start(p-1));
  theta(p-2) = std::log(start(p-2));
  double fopt;
  mle_treg f(y,x);
  int status = Numer::optim_lbfgs(f, theta, fopt, maxit, eps_f, eps_g);
  theta(p-1) = std::exp(theta(p-1));
  theta(p-2) = std::exp(theta(p-2));
  return Rcpp::List::create(
    Rcpp::Named("par") = theta,
    Rcpp::Named("value") = fopt,
    Rcpp::Named("conv") = status
  );
}

Rcpp::List of_mle_treg(
    Eigen::VectorXd& theta,
    Eigen::VectorXd& y,
    Eigen::MatrixXd& x
){
  double fopt;
  unsigned int p = theta.size();
  Eigen::VectorXd gr(p);
  mle_treg f(y,x);
  fopt = f.f_grad(theta,gr);
  return Rcpp::List::create(
    Rcpp::Named("f") = fopt,
    Rcpp::Named("grad") = gr
  );
}

// ------------------
// Switched Z-estimator
// ------------------

class swiz_treg: public Numer::MFuncGrad
{
private:
  const Eigen::VectorXd pi;
  const Eigen::MatrixXd x;
  const unsigned int seed;

public:
  swiz_treg(const Eigen::VectorXd pi_, const Eigen::MatrixXd x_, const unsigned int seed_) :
  pi(pi_), x(x_), seed(seed_) {}
  double f_grad(Numer::Constvec& theta, Numer::Refvec gr);
};

double swiz_treg::f_grad(
    Numer::Constvec& theta,
    Numer::Refvec gr
){
  unsigned int n = pi.size();
  Eigen::VectorXd p(n);
  Eigen::MatrixXd dp(n,n);
  p = psi_treg(theta,pi,x,seed);
  dp = d_psi_treg(theta,pi,x,seed);
  const double of = p.squaredNorm() / 0.2e1;
  gr = dp.transpose() * p;
  return of;
}

Rcpp::List optim_swiz_treg(
    Eigen::VectorXd& theta,
    Eigen::VectorXd& pi,
    Eigen::MatrixXd& x,
    unsigned int seed
){
  double fopt;
  // Eigen::Vector4d gr;
  swiz_treg f(pi,x,seed);
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

Rcpp::List of_swiz2_treg(
    Eigen::VectorXd& theta,
    Eigen::VectorXd& pi,
    Eigen::MatrixXd& x,
    unsigned int seed
){
  unsigned int n = theta.size();
  Eigen::VectorXd p(n);
  Eigen::MatrixXd dp(n,n);
  p = psi_treg(theta,pi,x,seed);
  dp = d_psi_treg(theta,pi,x,seed);

  return Rcpp::List::create(
    Rcpp::Named("psi") = p,
    Rcpp::Named("jac") = dp
  );
}


Rcpp::List of_swiz_treg(
    Eigen::VectorXd& theta,
    Eigen::VectorXd& pi,
    Eigen::MatrixXd& x,
    unsigned int seed
){
  double fopt;
  unsigned int p = theta.size();
  Eigen::VectorXd gr(p);
  swiz_treg f(pi,x,seed);
  fopt = f.f_grad(theta,gr);

  return Rcpp::List::create(
    Rcpp::Named("f") = fopt,
    Rcpp::Named("grad") = gr
  );
}

//' @title SwiZ distribution for t regression
//'
//' @param pi initial estimator
//' @param x matrix of design
//' @param B number of SwiZ estimates
//' @param seed integer representing the state fir random number generation
//' @param ncores number of cores (OpenMP parallelisation)
// [[Rcpp::export]]
Eigen::MatrixXd swiz_dist_treg(
    Eigen::VectorXd& pi,
    Eigen::MatrixXd& x,
    unsigned int B,
    unsigned int seed,
    unsigned int ncores
){
  unsigned int p = pi.size();
  Eigen::MatrixXd boot(B,p);

  #pragma omp parallel for num_threads(ncores)
  for(unsigned int i=0; i<B; ++i){
    Eigen::VectorXd theta = pi;
    theta(p-1) = std::log(pi(p-1));
    theta(p-2) = std::log(pi(p-2));
    unsigned int se = seed + i;
    double fopt;
    swiz_treg f(theta,x,se);
    int maxit = 500;
    double eps_f = 1e-10;
    double eps_g = 1e-10;
    Numer::optim_lbfgs(f, theta, fopt, maxit, eps_f, eps_g);
    boot.row(i) = theta;
    boot(i,p-1) = std::exp(theta(p-1));
    boot(i,p-2) = std::exp(theta(p-2));
  }

  return boot;
}

// ------------------
// Parametric Bootstrap
// ------------------

//' @title Parametric bootstrap distribution for t regression
//'
//' @param start initial estimator
//' @param x matrix of design
//' @param B number of SwiZ estimates
//' @param seed integer representing the state fir random number generation
//' @param ncores number of cores (OpenMP parallelisation)
// [[Rcpp::export]]
Eigen::MatrixXd par_bootstrap_mle_treg(
    Eigen::VectorXd& start,
    Eigen::MatrixXd& x,
    unsigned int B,
    unsigned int seed,
    unsigned int ncores
){
  unsigned int p = start.size();
  Eigen::MatrixXd boot(B,p);

  #pragma omp parallel for num_threads(ncores)
  for(unsigned int i=0; i<B; ++i){
    Eigen::VectorXd theta = start;
    theta(p-1) = std::log(start(p-1));
    theta(p-2) = std::log(start(p-2));
    unsigned int se = seed + i;
    // Eigen::VectorXd yt = y(x,start.head(p-2),start(p-2),start(p-1),se);
    Eigen::VectorXd y = rtreg(x,start.head(p-2),start(p-2),start(p-1),se);
    double fopt;
    mle_treg f(y,x);
    int maxit = 300;
    double eps_f = 1e-8;
    double eps_g = 1e-7;
    Numer::optim_lbfgs(f, theta, fopt, maxit, eps_f, eps_g);
    boot.row(i) = theta;
    boot(i,p-1) = std::exp(theta(p-1));
    boot(i,p-2) = std::exp(theta(p-2));
  }

  return boot;
}

//' @title Studentized parametric bootstrap distribution for t regression
//'
//' @param theta initial estimator
//' @param boot matrix, parametric bootstrap
//' @param x matrix of design
//' @param B number of SwiZ estimates
//' @param seed integer representing the state fir random number generation
//' @param ncores number of cores (OpenMP parallelisation)
//' @param robust if true uses robust estimation of covariance
// [[Rcpp::export]]
Eigen::MatrixXd par_boott_treg(
    Eigen::VectorXd& theta,
    Eigen::MatrixXd& boot,
    Eigen::MatrixXd& x,
    unsigned int B,
    unsigned int seed,
    unsigned int ncores,
    bool robust = false
){
  unsigned int S = boot.rows();
  unsigned int p = boot.cols();
  Eigen::MatrixXd boott(S,p);
  Eigen::VectorXi se(S);
  se = sample_int(S,seed);

  #pragma omp parallel for num_threads(ncores)
  for(unsigned int i=0; i<S; ++i){
    Eigen::MatrixXd new_boot(B,p);
    Eigen::VectorXd start = boot.row(i);
    new_boot = par_bootstrap_mle_treg(start,x,B,se(i),1);
    Eigen::ArrayXd sd(p);
    if(robust){
      for(unsigned int j=0;j<p;++j){
        sd(j) = mad(new_boot.col(j));
      }
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
// inspired from R bootstrap::bcanon
//' @title BCa acceleration parameter for t regression
//'
//' @param start MLE
//' @param y observations
//' @param x matrix of design
// [[Rcpp::export]]
Eigen::VectorXd acceleration_treg(
    Eigen::VectorXd& start,
    Eigen::VectorXd& y,
    Eigen::MatrixXd& x
){
  unsigned int n = y.size();
  unsigned int p = x.cols();
  int maxit = 300;
  double eps_f = 1e-8;
  double eps_g = 1e-7;
  Eigen::ArrayXXd u(n,p);
  Eigen::ArrayXXd uu(n,p);
  Eigen::ArrayXd t2(p), t3(p);
  Eigen::VectorXd yy(n-1);
  Eigen::MatrixXd xx(n-1,p);

  for(unsigned int i(0);i<n;++i){
    // I exploit the fact data has no order (iid)
    yy = y.tail(n-1);
    xx = x.bottomRows(n-1);
    if(i != 0){
      yy(i-1) = y(0);
      xx.row(i-1) = xx.row(0);
    }

    Eigen::VectorXd theta = start;
    theta(p-1) = std::log(start(p-1));
    theta(p-2) = std::log(start(p-2));
    double fopt;
    mle_treg f(yy,xx);
    Numer::optim_lbfgs(f, theta, fopt, maxit, eps_f, eps_g);
    u.row(i) = theta;
    u(i,p-1) = std::exp(theta(p-1));
    u(i,p-2) = std::exp(theta(p-2));
  }
  for(unsigned int i(0);i<p;++i){
    uu.col(i) = u.col(i).sum() / n - u.col(i);
  }
  uu *= uu;
  t2 = uu.colwise().sum();
  uu *= uu;
  t3 = uu.colwise().sum();
  return t3 / 6.0 / t2.pow(3.0 / 2.0);
}
