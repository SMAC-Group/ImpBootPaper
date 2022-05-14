#' @title Control for bootstraps
#' @param B number of bootstrap replicates
#' @param seed seed for random number generator
#' @param nc number of cores (for parallel computing)
#' @param B2 number of second-layer bootstrap replicates (only for studentized bootstrap)
#' @param robust boolean, if \code{TRUE} uses the \code{\link[stats]{mad}}
#' estimator instead of \code{\link[stats]{sd}} (only for studentized bootstrap)
#' @export
boot.control <- function(B = 1e4, seed = 123L, nc = 1L, B2 = 100, robust = FALSE){
  # verification
  if(as.integer(B)<=0) stop("'B' must be a positive integer")
  if(as.integer(seed)<=0) stop("'seed' must be a positive integer")
  if(as.integer(nc)<=0) stop("'nc' must be a positive integer")
  if(as.integer(B2)<=0) stop("'B2' must be a positive integer")
  if(!is.logical(robust)) stop("'robust' must be logical")
  list(B=B,seed=seed,nc=nc,B2=B2,robust=robust)
}

#' @title Parametric bootstrap
#' @param y       observations
#' @param x       matrix of design
#' @param initial MLE on observations
#' @param model   either "\code{lomax}" for Lomax distribution,
#' "\code{treg}" for t regression, "\code{betareg}" for beta regression
#' @param control  see \code{\link{boot.control}}
#' @return a \code{matrix} of bootstrap estimates
#' @export
parboot <- function(y, x=NULL, initial, model, control=boot.control()){
  # verification
  control <- do.call("boot.control",control)
  if(!is.numeric(y)) stop("'y' must be numeric")
  if(!is.null(x)) if(!is.numeric(x)) stop("'x' must be numeric")
  if(!is.null(x)) if(!is.matrix(x)) stop("'x' must be a matrix")
  p <- length(initial)
  if(!model %in% c("lomax","treg","betareg")) stop(sprintf("'%s' is not implemented",model))

  # computation
  if(model == "lomax") boot <- par_bootstrap_mle_lomax(initial,length(y),control$B,control$seed,control$nc)
  if(model == "treg") boot <- par_bootstrap_mle_treg(initial,x,control$B,control$seed,control$nc)
  if(model == "betareg") boot <- par_bootstrap_mle_betareg(initial,x,control$B,control$seed,control$nc)
  boot
}

#' @title Implicit bootstrap distribution
#' @param y       observations
#' @param x       matrix of design
#' @param initial MLE on observations
#' @param model   either "\code{lomax}" for Lomax distribution,
#' "\code{treg}" for t regression, "\code{betareg}" for beta regression
#' @param control  see \code{\link{boot.control}}
#' @return a \code{matrix} of bootstrap estimates
#' @export
impboot <- function(y, x=NULL, initial, model, control=boot.control()){
  # verification
  control <- do.call("boot.control",control)
  if(!is.numeric(y)) stop("'y' must be numeric")
  if(!is.null(x)) if(!is.numeric(x)) stop("'x' must be numeric")
  if(!is.null(x)) if(!is.matrix(x)) stop("'x' must be a matrix")
  p <- length(initial)
  if(!model %in% c("lomax","treg","betareg")) stop(sprintf("'%s' is not implemented",model))

  # computation
  if(model == "lomax") boot <- swiz_dist_lomax(initial,length(y),control$B,control$seed,control$nc)
  if(model == "treg") boot <- swiz_dist_treg(initial,x,control$B,control$seed,control$nc)
  if(model == "betareg") boot <- swiz_dist_betareg(initial,x,control$B,control$seed,control$nc)
  boot
}

#' @title Simulation from a parametric model
#' @param n       number of observations
#' @param x       matrix of design
#' @param theta   a list of model parameters (see 'Details')
#' @param model   either "\code{lomax}" for Lomax distribution,
#' "\code{treg}" for t regression, "\code{betareg}" for beta regression
#' @param seed    seed for random number generation
#' @details
#' The list of \code{theta} depends on the model. The names of the element of the list
#' must follow this convention:
#' \itemize{
#'    \item \code{lomax}: positive parameters \code{b}
#'    and \code{q};
#'    \item \code{treg}: the regression coefficients \code{beta},
#'    the variance \code{sig2}\eqn{>0} and the
#'    degrees of freedom \code{nu}\eqn{>0};
#'    \item \code{betareg}: the regression coefficients \code{beta}
#'    and the precision parameter (positive) \code{phi}.
#' }
#' @return a \code{matrix} of bootstrap estimates
#' @export
rmodel <- function(n=NULL, x=NULL, theta, model, seed){
  # verification
  if(is.null(n) && is.null(x)) stop("at least 'n' or 'x' must be specify")
  if(!is.null(n)) if(as.integer(n)<0) stop("'n' must be positive")
  if(!is.null(x)) if(!is.numeric(x)) stop("'x' must be numeric")
  if(!is.null(x)) if(!is.matrix(x)) stop("'x' must be a matrix")
  if(!is.null(x) && !is.null(n)) if(n != nrow(x)) stop("''")
  if(!is.list(theta)) stop("'theta' must be provided as a list")
  if(!model %in% c("lomax","treg","betareg")) stop(sprintf("'%s' is not implemented",model))
  if(as.integer(seed)<0) stop("'seed' must be positive")

  # model specific verification
  if(model == "lomax") {
    if(length(theta) != 2) stop("'theta' must be of length 2")
    if(!all(names(theta)%in%c("b","q"))) stop("'theta' name error")
    if(as.double(theta$b)<0 || as.double(theta$q)<0) stop("'b' and 'q' must be positive")
    if(is.null(n)) stop("'n' must be provided")
  }
  if(model == "treg"){
    if(length(theta) != 3) stop("'theta' must be of length 3")
    if(!all(names(theta)%in%c("beta","sig2","nu"))) stop("'theta' name error")
    if(is.null(x)) stop("'x' must be provided")
    p <- ncol(x)
    if(length(theta$beta) != p) stop("'x' and 'beta' dimension mismactch")
    if(as.double(theta$sig2)<0) stop("'sig2' must be positive")
    if(as.double(theta$nu)<0) stop("'nu' must be positive")
  }
  if(model == "betareg"){
    if(length(theta) != 2) stop("'theta' must be of length 2")
    if(!all(names(theta)%in%c("beta","phi"))) stop("'theta' name error")
    if(is.null(x)) stop("'x' must be provided")
    p <- ncol(x)
    if(length(theta$beta) != p) stop("'x' and 'beta' dimension mismactch")
    if(as.double(theta$phi)<0) stop("'phi' must be positive")
  }

  # computation
  if(model == "lomax") simu <- rlomax(n, theta$b, theta$q, seed)
  if(model == "treg") simu <- rtreg(x, theta$beta, theta$sig2, theta$nu, seed)
  if(model == "betareg") simu <- rbetareg(x, theta$beta, theta$phi, seed)
  simu
}
