#------------------------
# Results
#------------------------
library(ImpBootPaper)
library(ggplot2)
method_name <- c("Implicit", "Percentile", "Studentized", "Studentized (rob)", "BCa")#, "Asymptotic", "Asymptotic (rob)")
model <- "treg"

#------------------------
# Different parameters
#------------------------
nu <- 3
k <- seq(5,40,by=5)
coverages <- ci_length <- matrix(nrow = length(method_name), ncol = length(k))
id <- 3:4

for(i in seq_along(k)){
  load(file = paste0(".simu/data/",model,"_n_250_p_",k[i],"_nu_",nu,".rds"))
  coverages[,i] <- sapply(res[1:5], function(x) mean(apply(x,1,function(y){check_ci(ci=y[id],theta=-1)})=="center",na.rm=T))
  ci_length[,i] <- sapply(res[1:5], function(x) median(apply(x,1,function(y){diff(y[id])}),na.rm=T))
}

ds <- data.frame(
  Coverage = c(coverages),
  Method = rep(method_name, length(k)),
  Dimension = rep(k, each = length(method_name))
)

ci <- data.frame(
  ymin = binom.test(x=9500,n=1e4,p=.95,conf.level = .99)$conf.int[1],
  ymax = binom.test(x=9500,n=1e4,p=.95,conf.level = .99)$conf.int[2],
  xmin = -Inf,
  xmax = Inf
)

p <- ggplot(ds, aes(x = Dimension, y = Coverage, col = Method)) +
  geom_rect(data = ci, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill = "grey80", alpha = 0.8, inherit.aes = FALSE) +
  geom_hline(yintercept = 0.95, size = 1, color = "grey50") +
  geom_point(size = 4) +
  geom_path() +
  theme_bw() +
  theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5)) +
  ggtitle(paste0("Parameter: beta1\nSample size: 250\nDegrees of freedom: ",nu)) +
  xlab("p")
p

ds <- data.frame(
  Length = c(ci_length),
  Method = rep(method_name, length(k)),
  Dimension = rep(k, each = length(method_name))
)

p <- ggplot(ds, aes(x = Dimension, y = Length, col = Method)) +
  geom_point(size = 4) +
  geom_path() +
  theme_bw() +
  theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5)) +
  ggtitle(paste0("Parameter: beta1\nSample size: 250\nDegrees of freedom: ",nu)) +
  xlab("p") +
  ylab("Median confidence interval length")
p
