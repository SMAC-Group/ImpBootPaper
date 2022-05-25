library(ImpBootPaper)
library(idar)
data("codex")

y <- scale(codex$auc)
intercept <- rep(1, nrow(codex))
x <- as.matrix(cbind(intercept, subset(codex, select=c("gender", "age", "bmi"))))
x[,3:4] <- scale(x[,3:4])
param <- list(beta = rep(0,4),
              sig2 = 4,
              nu = 30)

fit_sample <- mle_estimate(y, x, unlist(param), "treg")
theta_hat <- fit_sample$par
implicit <- impboot_interval(y, x, theta_hat, "treg", control = list(seed = 498732, nc = 4L))
boot <- parboot(y, x, theta_hat, "treg", control = list(seed = 498732, nc = 4L))
percentile <- parboot_interval(y, x, theta_hat, boot, "treg")
boott <- parboott_interval(y, x, theta_hat, "treg", boot = boot, control = list(seed = 498732, nc = 4L))
bca <- BCa_interval(y=y, x=x, initial=theta_hat, model="treg", boot = boot, control = list(seed = 498732))
classic <- asymptotic_interval(y, x, theta_hat, "treg", boot = boot, control = list(seed = 498732, nc = 4L))
