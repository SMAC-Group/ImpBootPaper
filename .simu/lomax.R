#---------------------------------
# Estimation of Lomax distribution
#---------------------------------
library(ImpBootPaper)

#------------------------
# Simulation setting
#------------------------
assumed_model <- as.character(Sys.getenv("MODEL"))
n <- as.integer(Sys.getenv("N")) # sample size
param <- as.double(Sys.getenv("PARAM"))
MC <- 1e4 # Monte Carlo simulations
ncores <- 1L
theta <- c(param, 1.0) # true parameter
param <- list(b = theta[1], q = theta[2])
x <- NULL
p <- 2
alpha <- c(.025, .975)
S <- 1e4 # number of values for approximating the distribution
S2 <- 1e2 # number of second layer of bootstrap for studentized bootstrap
# seeds for random generator
set.seed(round(1e3 * theta[1]))
random_number1 <- sample.int(1e6,1)
set.seed(7483 * n)
random_number2 <- sample.int(1e6,1)
set.seed(98323 + random_number1 + random_number2) # random hand typing algorithm
se1 <- sample.int(1e7, MC)
se2 <- sample.int(1e7, MC)
res <- list(
  bca = matrix(nrow = MC, ncol = 2 * p)
)

#------------------------
# slurm specs
#------------------------
n_array <- 1000
ind <- matrix(seq_len(MC), nr=n_array, byr=TRUE)
id_slurm <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

#------------------------
# Simulation
#------------------------
for(m in na.omit(ind[id_slurm,])){
  # Sample computation
  y <- rmodel(n, x, param, assumed_model, se1[m])

  # initial estimator
  fit_sample <- NULL
  try(fit_sample <- mle_estimate(y, x, unlist(param), assumed_model), silent=TRUE)
  if(is.null(fit_sample)) next
  if(any(is.na(fit_sample$par))) next
  theta_hat <- fit_sample$par

  # BCa bootstrap
  res$bca[m,] <- BCa_interval(y=y, initial=theta_hat, model=assumed_model, control = list(seed = se2[m]))

 # save results
  save(res, file=paste0("tmp/",assumed_model,"_id_",id_slurm,".rds"))
  cat(m,"\n")
}

