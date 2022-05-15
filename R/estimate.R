#' @title Control for BFGS
#' @param maxit max number of iteration
#' @param eps_f tolerance for objective function
#' @param eps_g tolerance for gradient
#' @seealso \code{\link[RcppNumerical]{fastLR}}
#' @export
bfgs.control <- function(maxit = 500, eps_f = 1e-10, eps_g = 1e-10){
  # verification
  if(as.integer(maxit)<=0) stop("'maxit' must be a positive integer")
  if(as.numeric(eps_f)<0) stop("'eps_f' must be nonnegative")
  if(as.numeric(eps_g)<0) stop("'eps_g' must be nonnegative")
  list(maxit=maxit,eps_f=eps_f,eps_g=eps_g)
}

#' @title Maximum likelihood estimator
#' @param y       observations
#' @param x       matrix of design
#' @param starting_value starting value for L-BFGS optimization routine
#' @param model   either "\code{lomax}" for Lomax distribution,
#' "\code{treg}" for t regression, "\code{betareg}" for beta regression
#' @param control control for estimation (see \code{\link{bfgs.control}})
#' @return a \code{list} with \code{par}
#' @export
mle_estimate <- function(y, x=NULL, starting_value, model, control=bfgs.control()){
  # verification
  control <- do.call("bfgs.control",control)
  if(!is.numeric(y)) stop("'y' must be numeric")
  if(!is.null(x)) if(!is.numeric(x)) stop("'x' must be numeric")
  if(!is.null(x)) if(!is.matrix(x)) stop("'x' must be a matrix")
  if(!is.null(x)) p <- ncol(x)
  if(!model %in% c("lomax","treg","betareg")) stop(sprintf("'%s' is not implemented",model))
  if(model == "lomax") if(length(starting_value)!=2) stop("'starting_value' has not the right dimension")
  if(model == "treg") if(length(starting_value)!=p+2) stop("'starting_value' has not the right dimension")
  if(model == "betareg") if(length(starting_value)!=p+1) stop("'starting_value' has not the right dimension")

  # computation
  if(model == "lomax") est <- optim_mle_lomax(starting_value, y, control$maxit, control$eps_f, control$eps_g)
  if(model == "treg") est <- optim_mle_treg(starting_value, y, x, control$maxit, control$eps_f, control$eps_g)
  if(model == "betareg") est <- optim_mle_betareg(starting_value, y, x, control$maxit, control$eps_f, control$eps_g)
  est
}
