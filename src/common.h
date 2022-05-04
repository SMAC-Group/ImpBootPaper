#ifndef _COMMON
#define _COMMON
// ------------------
// Headers
// ------------------
// [[Rcpp::depends(RcppEigen,RcppNumerical,BH)]]

// Enable C++14 via the Rcpp plugin
// [[Rcpp::plugins("cpp14")]]

// Libraries
#include <RcppEigen.h>
#include <random>
#include <math.h>
#include <numeric>
#include <RcppNumerical.h>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/distributions/gamma.hpp>

// Add a flag to enable OpenMP at compile time
// [[Rcpp::plugins(openmp)]]

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

// Change error handling for Boost functions
// Define a specific policy:
typedef boost::math::policies::policy<
  boost::math::policies::digits10<5>,
  boost::math::policies::overflow_error<boost::math::policies::ignore_error>
> my_policy;

// User written headers
#include "misc.h"

#endif
