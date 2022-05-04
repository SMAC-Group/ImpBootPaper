# ----------------
# General function
# ----------------

#' @title Verify the position of a parameter wrt a closed interval
#' @param ci a matrix with two columns and same number of rows as \code{theta}
#' @param theta a vector of parameters
#' @export
check_ci <- function(ci,theta){
    if(!is.matrix(ci)) ci <- as.matrix(ci)
    if(ncol(ci)!=2) ci <- t(ci)
    if(nrow(ci)!=length(theta)) stop("verify 'ci' and 'theta' dimensions match")
    apply(cbind(ci,theta), 1, function(x){
        if(x[3] <= min(x[1:2])){"left"}
        else if(x[3] >= max(x[1:2])){"right"}else{"center"}
    })
}

# ----------------
# Lomax
# ----------------
#' @title BCa interval for Lomax distribution
#' @param boot parametric bootstrap distribution
#' @param y observations
#' @param initial MLE
#' @param alpha level
#' @param which binary, 0:alpha 1:lambda
#' @export
#' @importFrom stats qnorm pnorm quantile
BCa_interval_lomax <- function(boot, y, initial, alpha=c(.025, .975), which){
  z0 <- qnorm(mean(boot[,which+1]<initial[which+1]))
  a <- acceleration_lomax(initial, y, which)
  alpha_cor <- pnorm(z0 + (z0 + qnorm(alpha)) / (1.0 - a * (z0 + qnorm(alpha))))
  quantile(boot[,which+1],probs=alpha_cor)
}


