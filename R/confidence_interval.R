#' @title Verify the position of a parameter wrt a closed interval
#' @param ci a matrix with two columns and same number of rows as \code{theta}
#' @param theta a vector of parameters
#' @export
check_ci <- function(ci,theta){
    if(any(is.na(ci))) return(NA_character_)
    if(!is.matrix(ci)) ci <- as.matrix(ci)
    if(ncol(ci)!=2) ci <- t(ci)
    if(nrow(ci)!=length(theta)) stop("verify 'ci' and 'theta' dimensions match")
    apply(cbind(ci,theta), 1, function(x){
        if(x[3] <= min(x[1:2])){"left"}
        else if(x[3] >= max(x[1:2])){"right"}else{"center"}
    })
}

#' @title Parametric BCa interval
#' @param y       observations
#' @param x       matrix of design
#' @param initial MLE
#' @param boot    a \code{matrix} of parametric bootstrap estimates,
#' if \code{NULL} it is estimated internally
#' @param model   either "\code{lomax}" for Lomax distribution,
#' "\code{treg}" for t regression, "\code{betareg}" for beta regression
#' @param alpha   levels
#' @param control  only used if \code{boot} is \code{NULL}, see \code{\link{boot.control}}
#' @return a \code{matrix} levels in rows and MLE in columns
#' @export
#' @seealso \code{\link[bootstrap]{bcanon}}
#' @importFrom stats qnorm pnorm quantile
BCa_interval <- function(y, x=NULL, initial, boot=NULL, model, alpha=c(.025, .975), control=boot.control()){
  # verification
  control <- do.call("boot.control",control)
  if(!is.numeric(y)) stop("'y' must be numeric")
  if(!is.null(x)) if(!is.numeric(x)) stop("'x' must be numeric")
  if(!is.null(x)) if(!is.matrix(x)) stop("'x' must be a matrix")
  p <- length(initial)
  if(!is.null(boot)) if(!is.matrix(boot)) stop("'boot' must be a matrix")
  if(!is.null(boot)) if(ncol(boot)!=p) stop("'boot' dimension mistaches 'initial'")
  if(!model %in% c("lomax","treg","betareg")) stop(sprintf("'%s' is not implemented",model))
  if(any(alpha<=0) || any(alpha>=1)) stop("'alpha' must lie in (0,1)")

  # computation
  if(is.null(boot)){
    if(model == "lomax") boot <- par_bootstrap_mle_lomax(initial,length(y),control$B,control$seed,control$nc)
    if(model == "treg") boot <- par_bootstrap_mle_treg(initial,x,control$B,control$seed,control$nc)
    if(model == "betareg") boot <- par_bootstrap_mle_betareg(initial,x,control$B,control$seed,control$nc)
  }
  if(model == "lomax") acc <- acceleration_lomax(initial, y)
  if(model == "treg") acc <- acceleration_treg(initial, y, x)
  if(model == "betareg") acc <- acceleration_betareg(initial, y, x)
  zalpha <- qnorm(alpha)
  z0 <- qnorm(rowMeans(t(boot)<initial))
  ci <- matrix(nrow = length(alpha), ncol = p)
  for(i in seq_len(p)){
    tt <- pnorm(z0[i] + (z0[i] + zalpha) / (1.0 - acc[i] * (z0[i] + zalpha)))
    ci[,i] <- quantile(boot[,i], probs = tt, na.rm=TRUE)
  }
  ci
}

#' @title Percentile parametric bootstrap interval
#' @param y       observations
#' @param x       matrix of design
#' @param initial MLE on observations
#' @param boot    a \code{matrix} of parametric bootstrap estimates,
#' if \code{NULL} it is estimated internally
#' @param model   either "\code{lomax}" for Lomax distribution,
#' "\code{treg}" for t regression, "\code{betareg}" for beta regression
#' @param alpha   levels
#' @param control  see \code{\link{boot.control}}
#' @return a \code{matrix} levels in rows and MLE in columns
#' @export
parboot_interval <- function(y, x=NULL, initial, boot=NULL, model, alpha=c(.025, .975), control=boot.control()){
  # verification
  control <- do.call("boot.control",control)
  if(!is.numeric(y)) stop("'y' must be numeric")
  if(!is.null(x)) if(!is.numeric(x)) stop("'x' must be numeric")
  if(!is.null(x)) if(!is.matrix(x)) stop("'x' must be a matrix")
  p <- length(initial)
  if(!is.null(boot)) if(!is.matrix(boot)) stop("'boot' must be a matrix")
  if(!is.null(boot)) if(ncol(boot)!=p) stop("'boot' dimension mistaches 'initial'")
  if(!model %in% c("lomax","treg","betareg")) stop(sprintf("'%s' is not implemented",model))
  if(any(alpha<=0) || any(alpha>=1)) stop("'alpha' must lie in (0,1)")

  # computation
  if(is.null(boot)){
    if(model == "lomax") boot <- par_bootstrap_mle_lomax(initial,length(y),control$B,control$seed,control$nc)
    if(model == "treg") boot <- par_bootstrap_mle_treg(initial,x,control$B,control$seed,control$nc)
    if(model == "betareg") boot <- par_bootstrap_mle_betareg(initial,x,control$B,control$seed,control$nc)
  }
  ci <- apply(boot,2,function(x)quantile(x,probs=alpha,na.rm=T))
  ci
}

#' @title Studentized percentile parametric bootstrap interval
#' @param y       observations
#' @param x       matrix of design
#' @param initial MLE on observations
#' @param model   either "\code{lomax}" for Lomax distribution,
#' "\code{treg}" for t regression, "\code{betareg}" for beta regression
#' @param alpha   levels
#' @param boot    a \code{matrix} of parametric bootstrap estimates,
#' if \code{NULL} it is estimated internally
#' @param control  see \code{\link{boot.control}}
#' @return a \code{matrix} levels in rows and MLE in columns
#' @importFrom stats mad sd
#' @export
parboott_interval <- function(y, x=NULL, initial, model, alpha=c(.025, .975), boot=NULL, control=boot.control()){
  # verification
  control <- do.call("boot.control",control)
  if(!is.numeric(y)) stop("'y' must be numeric")
  if(!is.null(x)) if(!is.numeric(x)) stop("'x' must be numeric")
  if(!is.null(x)) if(!is.matrix(x)) stop("'x' must be a matrix")
  p <- length(initial)
  if(!model %in% c("lomax","treg","betareg")) stop(sprintf("'%s' is not implemented",model))
  if(!is.null(boot)) if(!is.matrix(boot)) stop("'boot' must be a matrix")
  if(!is.null(boot)) if(ncol(boot)!=p) stop("'boot' dimension mistaches 'initial'")
  if(any(alpha<=0) || any(alpha>=1)) stop("'alpha' must lie in (0,1)")

  # computation
  if(is.null(boot)){
    if(model == "lomax") boot <- par_bootstrap_mle_lomax(initial,length(y),control$B,control$seed,control$nc)
    if(model == "treg") boot <- par_bootstrap_mle_treg(initial,x,control$B,control$seed,control$nc)
    if(model == "betareg") boot <- par_bootstrap_mle_betareg(initial,x,control$B,control$seed,control$nc)
  }
  if(model == "lomax") boott <- par_boott_lomax(initial,boot,length(y),control$B2,control$seed,control$nc,control$robust)
  if(model == "treg") boott <- par_boott_treg(initial,boot,x,control$B2,control$seed,control$nc,control$robust)
  if(model == "betareg") boott <- par_boott_betareg(initial,boot,x,control$B2,control$seed,control$nc,control$robust)
  if(!control$robust) ci <- t(initial - apply(boot,2,sd,na.rm=T) * t(apply(boott,2,function(x)quantile(x,probs=alpha[length(alpha):1],na.rm=T))))
  if(control$robust) ci <- t(initial - apply(boot,2,mad,na.rm=T) * t(apply(boott,2,function(x)quantile(x,probs=alpha[length(alpha):1],na.rm=T))))
  ci
}

#' @title Implicit bootstrap interval
#' @param y       observations
#' @param x       matrix of design
#' @param initial MLE on observations
#' @param model   either "\code{lomax}" for Lomax distribution,
#' "\code{treg}" for t regression, "\code{betareg}" for beta regression
#' @param alpha   levels
#' @param control  see \code{\link{boot.control}}
#' @return a \code{matrix} levels in rows and MLE in columns
#' @export
impboot_interval <- function(y, x=NULL, initial, model, alpha=c(.025, .975), control=boot.control()){
  # verification
  control <- do.call("boot.control",control)
  if(!is.numeric(y)) stop("'y' must be numeric")
  if(!is.null(x)) if(!is.numeric(x)) stop("'x' must be numeric")
  if(!is.null(x)) if(!is.matrix(x)) stop("'x' must be a matrix")
  p <- length(initial)
  if(!model %in% c("lomax","treg","betareg")) stop(sprintf("'%s' is not implemented",model))
  if(any(alpha<=0) || any(alpha>=1)) stop("'alpha' must lie in (0,1)")

  # computation
  if(model == "lomax") boot <- swiz_dist_lomax(initial,length(y),control$B,control$seed,control$nc)
  if(model == "treg") boot <- swiz_dist_treg(initial,x,control$B,control$seed,control$nc)
  if(model == "betareg") boot <- swiz_dist_betareg(initial,x,control$B,control$seed,control$nc)
  ci <- apply(boot,2,function(x)quantile(x,probs=alpha,na.rm=T))
  ci
}

#' @title Asymptotic interval with parametric bootstrap covariance estimate
#' @param y       observations
#' @param x       matrix of design
#' @param initial MLE on observations
#' @param model   either "\code{lomax}" for Lomax distribution,
#' "\code{treg}" for t regression, "\code{betareg}" for beta regression
#' @param alpha   levels
#' @param boot    a \code{matrix} of parametric bootstrap estimates,
#' if \code{NULL} it is estimated internally
#' @param control  see \code{\link{boot.control}}
#' @return a \code{matrix} levels in rows and MLE in columns
#' @export
asymptotic_interval <- function(y, x=NULL, initial, model, alpha=c(.025, .975), boot=NULL, control=boot.control()){
  # verification
  control <- do.call("boot.control",control)
  if(!is.numeric(y)) stop("'y' must be numeric")
  if(!is.null(x)) if(!is.numeric(x)) stop("'x' must be numeric")
  if(!is.null(x)) if(!is.matrix(x)) stop("'x' must be a matrix")
  p <- length(initial)
  if(!model %in% c("lomax","treg","betareg")) stop(sprintf("'%s' is not implemented",model))
  if(!is.null(boot)) if(!is.matrix(boot)) stop("'boot' must be a matrix")
  if(!is.null(boot)) if(ncol(boot)!=p) stop("'boot' dimension mistaches 'initial'")
  if(any(alpha<=0) || any(alpha>=1)) stop("'alpha' must lie in (0,1)")

  # computation
  if(is.null(boot)){
    if(model == "lomax") boot <- par_bootstrap_mle_lomax(initial,length(y),control$B,control$seed,control$nc)
    if(model == "treg") boot <- par_bootstrap_mle_treg(initial,x,control$B,control$seed,control$nc)
    if(model == "betareg") boot <- par_bootstrap_mle_betareg(initial,x,control$B,control$seed,control$nc)
  }
  if(!control$robust) ci <- t(initial + apply(boot,2,sd,na.rm=T) * matrix(qnorm(p = alpha),nrow=2,ncol=length(initial),byrow=T))
  if(control$robust) ci <- t(initial + apply(boot,2,mad,na.rm=T) * matrix(qnorm(p = alpha),nrow=2,ncol=length(initial),byrow=T))
  ci
}
