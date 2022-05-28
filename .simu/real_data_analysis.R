library(ImpBootPaper)
library(idar)
library(MASS)
library(VGAM)
data("codex")

robscale <- function(x, center = TRUE, scale = TRUE) UseMethod("robscale")
robscale.default <- function(x, center = TRUE, scale = TRUE){
  ifelse(center,m<-median(x, na.rm = TRUE),m<-0)
  ifelse(scale,s<-mad(x, na.rm=TRUE),s<-1)
  (x - m) / s
}
robscale.matrix <- function(x, center = TRUE, scale = TRUE){
  apply(x,MARGIN=2,FUN=robscale.default)
}

# variable selection
selected_variables <- c("gender", "bmi", "t1_2", "age")
interaction <- list(c("gender", "bmi"))
y <- robscale(codex$auc)
# y <- scale(codex$auc)
intercept <- rep(1, nrow(codex))
x <- as.matrix(subset(codex, select=selected_variables))
x[,selected_variables[-1]] <- robscale(x[,selected_variables[-1]])
if(length(interaction)!=0){
  inter <- sapply(interaction, FUN=function(y) apply(x[,y],1,prod))
  colnames(inter) <- sapply(interaction, paste, collapse = "_")
  x <- cbind(x, inter)
}

# modeling
fit_lm <- lm(y~x)
fit_lmrob <- rlm(y~x)
start <- c(coef(fit_lmrob), log(sigma(fit_lmrob)), log(30))
fit_vgam <- vglm(y ~ x, family = studentt3(ldf = "loglink"), coefstart = start)
x <- cbind(intercept,x)
param <- list(beta = coef(fit_lmrob),
              sig2 = sigma(fit_lmrob)^2,
              nu = 30)
control <- list(seed = 4987, nc = 1L)

fit_sample <- mle_estimate(y, x, unlist(param), "treg")
theta_hat <- fit_sample$par

res <- cbind(coef(fit_lm),coef(fit_lmrob),coef(fit_vgam,matrix=T)[,1],head(theta_hat,-2))
colnames(res) <- c("LM", "RLM", "VGAM", "treg")
res

extra <- cbind(c(exp(coef(fit_vgam,matrix=T)[1,2])^2,exp(coef(fit_vgam,matrix=T)[1,3])),tail(theta_hat,2))
colnames(extra) <- c("VGAM", "treg")
rownames(extra) <- c("sig2", "nu")
extra

resid <- c(y - tcrossprod(head(theta_hat,-2), x))
qqnorm(resid)
qqline(resid, col = 2)
qqplot(resid, sqrt(theta_hat[length(theta_hat)-1]) * rt(length(y), df = theta_hat[length(theta_hat)]))

# confidence intervals
class_lm <- confint(fit_lm)
class_vgam <- confint(fit_vgam)

implicit <- t(impboot_interval(y, x, theta_hat, "treg", control = control))
rownames(implicit) <- c(colnames(x), "sig2", "nu")
implicit

boot <- parboot(y, x, theta_hat, "treg", control = control)
percentile <- t(parboot_interval(y, x, theta_hat, boot, "treg"))
rownames(percentile) <- c(colnames(x), "sig2", "nu")
percentile

boott <- t(parboott_interval(y, x, theta_hat, "treg", boot = boot, control = control))
rownames(boott) <- c(colnames(x), "sig2", "nu")
boott
bca <- t(BCa_interval(y=y, x=x, initial=theta_hat, model="treg", boot = boot, control = control))
rownames(bca) <- c(colnames(x), "sig2", "nu")
bca
classic <- t(asymptotic_interval(y, x, theta_hat, "treg", boot = boot, control = control))
rownames(classic) <- c(colnames(x), "sig2", "nu")
classic

