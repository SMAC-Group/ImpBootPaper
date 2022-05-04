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
double r_gamma(
    double alpha,
    unsigned int seed
){
  if(alpha < 0.0){
    return NAN;
  }

  boost::math::gamma_distribution<>  gamma_distr(alpha);
  std::mt19937_64 engine(seed);  // Mersenne twister random number engine
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  double u = unif(engine);

  return quantile(gamma_distr,u);
}

double d_r_gamma_a(
    double alpha,
    unsigned int seed
){
  if(alpha < 0.0){
    return NAN;
  }

  double eps = 1e-8;

  boost::math::gamma_distribution<>  gamma_distr_u(alpha+eps);
  boost::math::gamma_distribution<>  gamma_distr_l(alpha-eps);
  std::mt19937_64 engine(seed);  // Mersenne twister random number engine
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  double u = unif(engine);

  return (quantile(gamma_distr_u,u) - quantile(gamma_distr_l,u)) / 2.0 / eps;
}

Eigen::VectorXd y(
    const Eigen::MatrixXd& x,
    const Eigen::VectorXd& beta,
    double varphi,
    unsigned int seed
){
  unsigned int n = x.rows();
  Eigen::VectorXd y(n), eta(n);
  double phi,a,b,mu,p,q;

  phi = std::exp(varphi);
  eta = x * beta;
  Eigen::VectorXi seeds = sample_int(n, seed);

  for(unsigned int i(0);i<n;++i){
    mu = logistic(eta(i));
    p = mu * phi;
    q = phi - p;
    a = r_gamma(p,seeds(i));
    b = r_gamma(q,seeds(i)+1);
    y(i) = a / (b + a);
  }

  return y;
}

Eigen::MatrixXd d_y(
    const Eigen::MatrixXd& x,
    const Eigen::VectorXd& beta,
    double varphi,
    unsigned int seed
){
  unsigned int n = x.rows();
  unsigned int i1 = beta.size() + 1;
  Eigen::VectorXd eta(n);
  Eigen::MatrixXd dy(n,i1);
  double phi,mu,a,b,p,q,t2,t3,t4,t5,t6;
  Eigen::RowVectorXd v(i1),d_a(i1),d_b(i1);

  phi = std::exp(varphi);
  eta = x * beta;
  Eigen::VectorXi seeds = sample_int(n, seed);

  for(unsigned int i(0);i<n;++i){
    mu = logistic(eta(i));
    p = mu * phi;
    q = phi - p;
    a = r_gamma(p,seeds(i));
    b = r_gamma(q,seeds(i)+1);
    t3 = d_r_gamma_a(p,seeds(i));
    t4 = d_r_gamma_a(q,seeds(i)+1);
    t2 = a + b;
    t6 = 0.1e1 - mu;
    t5 = mu * t6;
    d_a.head(i1-1) = t3 * phi * t5 * x.row(i);
    d_a(i1-1) = mu * t3 * phi;
    d_b.head(i1-1) = -t4 * phi * t5 * x.row(i);
    d_b(i1-1) = phi * t4 * t6;
    v = b * d_a - a * d_b;
    dy.row(i) = v / t2 / t2;
  }

  return dy;
}


// [[Rcpp::export]]
Eigen::VectorXd rbetareg(
    const Eigen::MatrixXd& x,
    const Eigen::VectorXd& beta,
    double phi,
    unsigned int seed
){
  return y(x,beta,log(phi),seed);
}

// ------------------
// Log-likelihood
// ------------------
double ll(
    const Eigen::VectorXd& theta,
    const Eigen::VectorXd& y,
    const Eigen::MatrixXd& x
){
  unsigned int n = y.size();
  unsigned int p = theta.size();
  Eigen::VectorXd eta(n);
  double of(0.0),phi,mu,t1,t2;

  phi = std::exp(theta(p-1));
  eta = x * theta.head(p-1);

  for(unsigned int i(0);i<n;++i){
    mu = logistic(eta(i));
    t1 = mu * phi;
    t2 = phi - t1;
    of += std::lgamma(t1) - std::lgamma(phi) + std::lgamma(t2) - (t1 - 0.1e1) * std::log(y(i)) - (t2 - 0.1e1) * std::log(0.1e1 - y(i));
  }

  return of / n;
}

Eigen::VectorXd d_ll(
    const Eigen::VectorXd& theta,
    const Eigen::VectorXd& y,
    const Eigen::MatrixXd& x
){
  unsigned int n = y.size();
  unsigned int p = theta.size();
  Eigen::VectorXd eta(n);
  double phi,mu,t1,t2,t3,t4,delta;
  Eigen::VectorXd grad = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd d_mu(p-1);

  phi = std::exp(theta(p-1));
  eta = x * theta.head(p-1);
  t4 = boost::math::digamma(phi, my_policy());

  for(unsigned int i(0);i<n;++i){
    mu = logistic(eta(i));
    t1 = mu * phi;
    t2 = phi - t1;
    t3 = mu - mu * mu;
    delta = boost::math::digamma(t1, my_policy()) - boost::math::digamma(t2, my_policy()) - logit(y(i));
    d_mu = t3 * x.row(i);
    grad.head(p-1) += phi * delta * d_mu;
    grad(p-1) += mu * delta + boost::math::digamma(t2, my_policy()) - std::log(0.1e1 - y(i)) - t4;
  }

  grad(p-1) *= phi;

  return grad / n;
}

Eigen::VectorXd psi(
    const Eigen::VectorXd& theta,
    const Eigen::VectorXd& pi,
    const Eigen::MatrixXd& x,
    unsigned int seed
){
  unsigned int n = x.rows();
  unsigned int p = theta.size();
  Eigen::VectorXd eta(n);
  double phi,mu,t1,t2,t3,t4,t5,delta;
  Eigen::VectorXd grad = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd d_mu(p-1);
  Eigen::Vector2d tmp;

  // return NAs if theta has NAs
  if(theta.array().isFinite().count() != p){
    grad = Eigen::VectorXd::Constant(p, NAN);
    return grad;
  }

  phi = std::exp(pi(p-1));
  eta = x * pi.head(p-1);

  // Rcpp::Rcout << "Parameters: " << theta.transpose() << "\n";
  // Rcpp::Rcout << "Phi: " << exp(theta(p-1)) << "\n";

  double _tmp_log_min = logistic((x * theta.head(p-1)).minCoeff());
  double _tmp_log_max = logistic((x * theta.head(p-1)).maxCoeff());
  double _phi = exp(theta(p-1));

  // return NAs if eta is constant
  if(_tmp_log_min == _tmp_log_max){
    grad = Eigen::VectorXd::Constant(p, NAN);
    return grad;
  }

  double eps = 1e-8;
  double _tmp1 = _tmp_log_min * _phi;
  double _tmp2 = _phi - _tmp_log_max * _phi;
  tmp(0) = _tmp1 - eps; tmp(1) = _tmp2 - eps;
  double _min_observed = tmp.minCoeff();

  // Rcpp::Rcout << "min observed: " << tmp.transpose() << "\n";

  _tmp1 = _tmp_log_max * _phi;
  _tmp2 = _phi - _tmp_log_min * _phi;
  tmp(0) = _tmp1 + eps; tmp(1) = _tmp2 + eps;
  double _max_observed = tmp.maxCoeff();

  // Rcpp::Rcout << "max observed: " << tmp.transpose() << "\n";


  // Rcpp::Rcout << "min observed: " << _min_observed << "\n";
  // Rcpp::Rcout << "max observed: " << _max_observed << "\n";
  // Rcpp::Rcout << "min observed: " << std::numeric_limits<double>::min() << "\n";
  // Rcpp::Rcout << "Error check " << ((_min_observed < 1e-7) || (!std::isfinite(_min_observed))) << "\n";

  // return NAs if parameter for gamma is smaller than numerical zero
  if(_min_observed < 1e-7 || !std::isfinite(_min_observed) ||
     _max_observed > 1e7 || !std::isfinite(_max_observed)){//1.175494e-38){
    grad = Eigen::VectorXd::Constant(p, NAN);
    return grad;
  }

  Eigen::VectorXd yt = y(x,theta.head(p-1),theta(p-1),seed);

  t4 = boost::math::digamma(phi, my_policy());

  for(unsigned int i(0);i<n;++i){
    mu = logistic(eta(i));
    t1 = mu * phi;
    t2 = phi - t1;
    t3 = mu - mu * mu;
    t5 = boost::math::digamma(t2, my_policy());
    delta = boost::math::digamma(t1, my_policy()) - t5 - logit(yt(i));
    d_mu = t3 * x.row(i);
    grad.head(p-1) += phi * delta * d_mu;
    grad(p-1) += mu * delta + t5 - std::log(0.1e1 - yt(i)) - t4; // check the log computation
  }

  grad(p-1) *= phi;

  return grad / n;
}

Eigen::MatrixXd d_psi(
    const Eigen::VectorXd& theta,
    const Eigen::VectorXd& pi,
    const Eigen::MatrixXd& x,
    unsigned int seed
){
  unsigned int n = x.rows();
  unsigned int p = theta.size();
  Eigen::VectorXd eta(n);
  double phi,mu,t3,t4;
  Eigen::MatrixXd jac = Eigen::MatrixXd::Zero(p,p);
  Eigen::Vector2d tmp;

  // return NAs if theta has NAs
  if(theta.array().isFinite().count() != p){
    jac = Eigen::MatrixXd::Constant(p, p, NAN);
    return jac;
  }

  phi = std::exp(pi(p-1));
  eta = x * pi.head(p-1);

  double _tmp_log_min = logistic((x * theta.head(p-1)).minCoeff());
  double _tmp_log_max = logistic((x * theta.head(p-1)).maxCoeff());
  double _phi = std::exp(theta(p-1));

  // return NAs if eta is constant
  if(_tmp_log_min == _tmp_log_max){
    jac = Eigen::MatrixXd::Constant(p, p, NAN);
    // return jac;
  }

  double eps = 1e-8;
  double _tmp1 = _tmp_log_min * _phi;
  double _tmp2 = _phi - _tmp_log_max * _phi;
  tmp(0) = _tmp1 - eps; tmp(1) = _tmp2 - eps;
  double _min_observed = tmp.minCoeff();

  _tmp1 = _tmp_log_max * _phi;
  _tmp2 = _phi - _tmp_log_min * _phi;
  tmp(0) = _tmp1 + eps; tmp(1) = _tmp2 + eps;
  double _max_observed = tmp.maxCoeff();

  // return NAs if parameter for gamma is smaller than numerical zero
  if(_min_observed < 1e-7 || !std::isfinite(_min_observed) ||
     _max_observed > 1e7 || !std::isfinite(_max_observed)){//1.175494e-38){
    jac = Eigen::MatrixXd::Constant(p, p, NAN);
    return jac;
  }

  Eigen::VectorXd d_mu(p-1);
  Eigen::VectorXd yt = y(x,theta.head(p-1),theta(p-1),seed);
  Eigen::MatrixXd dyt = d_y(x,theta.head(p-1),theta(p-1),seed);
  Eigen::RowVectorXd lambda(p);

  for(unsigned int i(0);i<n;++i){
    t4 = yt(i) - yt(i) * yt(i);
    lambda = dyt.row(i) / t4;
    mu = logistic(eta(i));
    t3 = mu - mu * mu;
    d_mu = t3 * x.row(i);
    jac.topRows(p-1) -= phi * d_mu * lambda;
    jac.row(p-1) += yt(i) * lambda - mu * lambda;
  }
  jac.row(p-1) *= phi;

  return jac / n;
}

// ------------------
// Maximum likelihood
// ------------------

class mle: public Numer::MFuncGrad
{
private:
  const Eigen::VectorXd y;
  const Eigen::MatrixXd x;

public:
  mle(const Eigen::VectorXd y_, const Eigen::MatrixXd x_) : y(y_), x(x_) {}
  double f_grad(Numer::Constvec& theta, Numer::Refvec gr);
};

double mle::f_grad(
    Numer::Constvec& theta,
    Numer::Refvec gr
){
  const double of = ll(theta,y,x);
  gr = d_ll(theta,y,x);
  return of;
}

// [[Rcpp::export]]
Rcpp::List optim_mle(
    Eigen::VectorXd& theta,
    Eigen::VectorXd& y,
    Eigen::MatrixXd& x,
    int maxit = 300,
    double eps_f = 1e-6,
    double eps_g = 1e-6
){
  double fopt;
  mle f(y,x);
  int status = Numer::optim_lbfgs(f, theta, fopt, maxit, eps_f, eps_g);
  return Rcpp::List::create(
    Rcpp::Named("par") = theta,
    Rcpp::Named("value") = fopt,
    Rcpp::Named("conv") = status
  );
}

// [[Rcpp::export]]
Rcpp::List of_mle(
    Eigen::VectorXd& theta,
    Eigen::VectorXd& y,
    Eigen::MatrixXd& x
){
  double fopt;
  unsigned int p = theta.size();
  Eigen::VectorXd gr(p);
  mle f(y,x);
  fopt = f.f_grad(theta,gr);
  return Rcpp::List::create(
    Rcpp::Named("f") = fopt,
    Rcpp::Named("grad") = gr
  );
}

// ------------------
// Switched Z-estimator
// ------------------

class swiz: public Numer::MFuncGrad
{
private:
  const Eigen::VectorXd pi;
  const Eigen::MatrixXd x;
  const unsigned int seed;

public:
  swiz(const Eigen::VectorXd pi_, const Eigen::MatrixXd x_, const unsigned int seed_) :
  pi(pi_), x(x_), seed(seed_) {}
  double f_grad(Numer::Constvec& theta, Numer::Refvec gr);
};

double swiz::f_grad(
    Numer::Constvec& theta,
    Numer::Refvec gr
){
  Eigen::VectorXd p = psi(theta,pi,x,seed);
  Eigen::MatrixXd dp = d_psi(theta,pi,x,seed);
  const double of = p.squaredNorm() / 0.2e1;
  gr = dp.transpose() * p;
  return of;
}

Rcpp::List psi_dpsi(
    Eigen::VectorXd& theta,
    Eigen::VectorXd& pi,
    Eigen::MatrixXd& x,
    unsigned int seed
){
  Eigen::VectorXd p = psi(theta,pi,x,seed);
  Eigen::MatrixXd dp = d_psi(theta,pi,x,seed);
  return Rcpp::List::create(
    Rcpp::Named("psi") = p,
    Rcpp::Named("d_psi") = dp
  );
}

// [[Rcpp::export]]
Rcpp::List of_swiz(
    Eigen::VectorXd& theta,
    Eigen::VectorXd& pi,
    Eigen::MatrixXd& x,
    unsigned int seed
){
  double fopt;
  unsigned int p = theta.size();
  Eigen::VectorXd gr(p);
  swiz f(pi,x,seed);
  fopt = f.f_grad(theta,gr);
  return Rcpp::List::create(
    Rcpp::Named("f") = fopt,
    Rcpp::Named("grad") = gr
  );
}

// [[Rcpp::export]]
double swiz_fn(
    Eigen::VectorXd& theta,
    Eigen::VectorXd& pi,
    Eigen::MatrixXd& x,
    unsigned int seed
){
  Eigen::VectorXd p = psi(theta,pi,x,seed);
  return p.squaredNorm() / 0.2e1;
}

// [[Rcpp::export]]
Eigen::VectorXd swiz_gr(
    Eigen::VectorXd& theta,
    Eigen::VectorXd& pi,
    Eigen::MatrixXd& x,
    unsigned int seed
){
  Eigen::VectorXd p = psi(theta,pi,x,seed);
  Eigen::MatrixXd dp = d_psi(theta,pi,x,seed);
  return dp.transpose() * p;
}

// [[Rcpp::export]]
Rcpp::List optim_swiz(
    Eigen::VectorXd& theta,
    Eigen::VectorXd& pi,
    Eigen::MatrixXd& x,
    unsigned int seed
){
  double fopt;
  swiz f(pi,x,seed);
  int maxit = 300;
  double eps_f = 1e-8;
  double eps_g = 1e-8;
  int status = Numer::optim_lbfgs(f, theta, fopt, maxit, eps_f, eps_g);
  return Rcpp::List::create(
    Rcpp::Named("par") = theta,
    Rcpp::Named("value") = fopt,
    Rcpp::Named("conv") = status
  );
}

// [[Rcpp::export]]
Eigen::MatrixXd swiz_dist(
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
    unsigned int se = seed + i;
    double fopt;
    swiz f(pi,x,se);
    int maxit = 300;
    double eps_f = 1e-8;
    double eps_g = 1e-8;
    Numer::optim_lbfgs(f, theta, fopt, maxit, eps_f, eps_g);
    boot.row(i) = theta;
  }

  return boot;
}

// ------------------
// Parametric Bootstrap
// ------------------

// [[Rcpp::export]]
Eigen::MatrixXd par_bootstrap_mle(
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
    unsigned int se = seed + i;
    Eigen::VectorXd yt = y(x,start.head(p-1),start(p-1),se);
    double fopt;
    mle f(yt,x);
    int maxit = 300;
    double eps_f = 1e-8;
    double eps_g = 1e-8;
    Numer::optim_lbfgs(f, theta, fopt, maxit, eps_f, eps_g);
    boot.row(i) = theta;
  }

  return boot;
}

// [[Rcpp::export]]
Eigen::MatrixXd par_boott(
    Eigen::VectorXd& theta,
    Eigen::MatrixXd& boot,
    Eigen::MatrixXd& x,
    unsigned int B,
    unsigned int seed,
    unsigned int ncores
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
    new_boot = par_bootstrap_mle(start,x,B,se(i),1);
    Eigen::MatrixXd cov = covariance(new_boot);
    Eigen::ArrayXd sd = cov.diagonal().array().sqrt();
    boott.row(i) = (start - theta).array() / sd;
  }

  return boott;
}
