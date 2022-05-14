#------------------------
# Results
#------------------------
library(ImpBootPaper)
library(ggplot2)
method_name <- c("Implicit", "Percentile", "Studentized", "Studentized (rob)", "BCa", "Asymptotic", "Asymptotic (rob)")

#------------------------
# Different parameters
#------------------------
n <- 50
param <- seq(.5,5,by=.5)
model <- "lomax"
coverages1 <- right_coverages1 <- ci_length1 <- coverages2 <- right_coverages2 <- ci_length2 <- matrix(nrow = length(method_name), ncol = length(param))

for(i in seq_along(param)){
  load(file = paste0(".simu/data/",model,"_n_",n,"_param_",param[i],".rds"))
  coverages1[,i] <- sapply(res[1:7], function(x) mean(apply(x,1,function(y){check_ci(ci=y[1:2],theta=param[i])})=="center",na.rm=T))
  right_coverages1[,i] <- sapply(res[1:7], function(x) mean(apply(x,1,function(y){check_ci(ci=y[1:2],theta=param[i])})=="right",na.rm=T))
  ci_length1[,i] <- sapply(res[1:7], function(x) median(apply(x,1,function(y){diff(y[1:2])}),na.rm=T))
  coverages2[,i] <- sapply(res[1:7], function(x) mean(apply(x,1,function(y){check_ci(ci=y[3:4],theta=1)})=="center",na.rm=T))
  right_coverages2[,i] <- sapply(res[1:7], function(x) mean(apply(x,1,function(y){check_ci(ci=y[3:4],theta=1)})=="right",na.rm=T))
  ci_length2[,i] <- sapply(res[1:7], function(x) median(apply(x,1,function(y){diff(y[3:4])}),na.rm=T))
}

library(ggplot2)
# coverage
ds1 <- data.frame(
  Coverage = c(coverages1),
  Method = rep(method_name, length(param)),
  Parameter = rep(param, each = length(method_name))
)
ds2 <- data.frame(
  Coverage = c(coverages2),
  Method = rep(method_name, length(param)),
  Parameter = rep(param, each = length(method_name))
)

ci <- data.frame(
  ymin = binom.test(x=9500,n=1e4,p=.95,conf.level = .99)$conf.int[1],
  ymax = binom.test(x=9500,n=1e4,p=.95,conf.level = .99)$conf.int[2],
  xmin = -Inf,
  xmax = Inf
)

p <- ggplot(ds1, aes(x = Parameter, y = Coverage, col = Method)) +
  geom_rect(data = ci, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill = "grey80", alpha = 0.8, inherit.aes = FALSE) +
  geom_hline(yintercept = 0.95, size = 1, color = "grey50") +
  geom_point(size = 4) +
  geom_path() +
  theme_bw() +
  theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Parameter: b\nSample size: 50") +
  xlab("b\n(q=1)")
p

p <- ggplot(ds2, aes(x = Parameter, y = Coverage, col = Method)) +
  geom_rect(data = ci, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill = "grey80", alpha = 0.8, inherit.aes = FALSE) +
  geom_hline(yintercept = 0.95, size = 1, color = "grey50") +
  geom_point(size = 4) +
  geom_path() +
  theme_bw() +
  theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Parameter: q\nSample size: 50") +
  xlab("b\n(q=1)")
p

# right

ds1 <- data.frame(
  Coverage = c(right_coverages1),
  Method = rep(method_name, length(param)),
  Parameter = rep(param, each = length(method_name))
)
ds2 <- data.frame(
  Coverage = c(right_coverages2),
  Method = rep(method_name, length(param)),
  Parameter = rep(param, each = length(method_name))
)

ci <- data.frame(
  ymin = binom.test(x=9500,n=1e4,p=.95,conf.level = .99)$conf.int[1],
  ymax = binom.test(x=9500,n=1e4,p=.95,conf.level = .99)$conf.int[2],
  xmin = -Inf,
  xmax = Inf
)

p <- ggplot(ds1, aes(x = Parameter, y = Coverage, col = Method)) +
  # geom_rect(data = ci, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill = "grey80", alpha = 0.8, inherit.aes = FALSE) +
  geom_hline(yintercept = 0.025, size = 1, color = "grey50") +
  geom_point(size = 4) +
  geom_path() +
  theme_bw() +
  theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Parameter: b\nSample size: 50") +
  xlab("b\n(q=1)")
p

p <- ggplot(ds2, aes(x = Parameter, y = Coverage, col = Method)) +
  # geom_rect(data = ci, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill = "grey80", alpha = 0.8, inherit.aes = FALSE) +
  geom_hline(yintercept = 0.025, size = 1, color = "grey50") +
  geom_point(size = 4) +
  geom_path() +
  theme_bw() +
  theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Parameter: q\nSample size: 50") +
  xlab("b\n(q=1)")
p

# confidence interval length
ds1 <- data.frame(
  Length = c(ci_length1[-c(3,6),]),
  Method = rep(method_name[-c(3,6)], length(param)),
  Parameter = rep(param, each = length(method_name[-c(3,6)]))
)
ds2 <- data.frame(
  Length = c(ci_length2[-c(3,6),]),
  Method = rep(method_name[-c(3,6)], length(param)),
  Parameter = rep(param, each = length(method_name[-c(3,6)]))
)

p <- ggplot(ds1, aes(x = Parameter, y = Length, col = Method)) +
  geom_point(size = 4) +
  geom_path() +
  theme_bw() +
  theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Parameter: b\nSample size: 50") +
  xlab("b\n(q=1)") +
  ylab("Median confidence interval length")
p

p <- ggplot(ds2, aes(x = Parameter, y = Length, col = Method)) +
  geom_point(size = 4) +
  geom_path() +
  theme_bw() +
  theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Parameter: q\nSample size: 50") +
  xlab("b\n(q=1)") +
  ylab("Median confidence interval length")
p

#------------------------
# Different sample size
#------------------------
n <- c(seq(25,200,by=25), seq(250,500,by=50))
param <- 1
model <- "lomax"
coverages1 <- ci_length1 <- coverages2 <- ci_length2 <- matrix(nrow = length(method_name), ncol = length(n))

for(i in seq_along(n)){
  load(file = paste0(".simu/data/",model,"_n_",n[i],"_param_",param,".rds"))
  coverages1[,i] <- sapply(res[1:7], function(x) mean(apply(x,1,function(y){check_ci(ci=y[1:2],theta=param)})=="center",na.rm=T))
  ci_length1[,i] <- sapply(res[1:7], function(x) median(apply(x,1,function(y){diff(y[1:2])}),na.rm=T))
  coverages2[,i] <- sapply(res[1:7], function(x) mean(apply(x,1,function(y){check_ci(ci=y[3:4],theta=1)})=="center",na.rm=T))
  ci_length2[,i] <- sapply(res[1:7], function(x) median(apply(x,1,function(y){diff(y[3:4])}),na.rm=T))
}

# coverage
ds1 <- data.frame(
  Coverage = c(coverages1),
  Method = rep(method_name, length(n)),
  Sample = rep(n, each = length(method_name))
)
ds2 <- data.frame(
  Coverage = c(coverages2),
  Method = rep(method_name, length(n)),
  Sample = rep(n, each = length(method_name))
)

ci <- data.frame(
  ymin = binom.test(x=9500,n=1e4,p=.95,conf.level = .99)$conf.int[1],
  ymax = binom.test(x=9500,n=1e4,p=.95,conf.level = .99)$conf.int[2],
  xmin = -Inf,
  xmax = Inf
)

p <- ggplot(ds1, aes(x = Sample, y = Coverage, col = Method)) +
  geom_rect(data = ci, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill = "grey80", alpha = 0.8, inherit.aes = FALSE) +
  geom_hline(yintercept = 0.95, size = 1, color = "grey50") +
  geom_point(size = 4) +
  geom_path() +
  theme_bw() +
  theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Parameter: b (b=1, q=1)") +
  xlab("Sample size")
p

p <- ggplot(ds2, aes(x = Sample, y = Coverage, col = Method)) +
  geom_rect(data = ci, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill = "grey80", alpha = 0.8, inherit.aes = FALSE) +
  geom_hline(yintercept = 0.95, size = 1, color = "grey50") +
  geom_point(size = 4) +
  geom_path() +
  theme_bw() +
  theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Parameter: q (b=1, q=1)") +
  xlab("Sample size")
p

# confidence interval length
ci_length1[c(3,6),1:2] <- NA_real_
ci_length1[2,1] <- NA_real_
ds1 <- data.frame(
  Length = c(ci_length1),
  Method = rep(method_name, length(param)),
  Sample = rep(n, each = length(method_name))
)
ci_length2[c(3,6),1:2] <- NA_real_
ci_length2[2,1] <- NA_real_
ds2 <- data.frame(
  Length = c(ci_length2),
  Method = rep(method_name, length(param)),
  Sample = rep(n, each = length(method_name))
)

p <- ggplot(ds1, aes(x = Sample, y = Length, col = Method)) +
  geom_point(size = 4) +
  geom_path() +
  theme_bw() +
  theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Parameter: b (b=1, q=1)") +
  xlab("Sample size") +
  ylab("Median confidence interval length")
p

p <- ggplot(ds2, aes(x = Sample, y = Length, col = Method)) +
  geom_point(size = 4) +
  geom_path() +
  theme_bw() +
  theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Parameter: q (b=1, q=1)") +
  xlab("Sample size") +
  ylab("Median confidence interval length")
p
