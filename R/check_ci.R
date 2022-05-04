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
