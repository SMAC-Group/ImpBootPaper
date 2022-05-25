#---------------------------------
# Simulation study
#---------------------------------
library(ImpBootPaper)

#------------------------
# Simulation setting
#------------------------
model <- as.character(Sys.getenv("MODEL"))
if(!model %in% c("lomax","treg","betareg")) stop(sprintf("'%s' is not implemented",model))
n <- as.integer(Sys.getenv("N")) # sample size
p <- as.integer(Sys.getenv("P"))
if(model == "lomax"){
  param <- list(b = as.double(Sys.getenv("b")),
                q = as.double(Sys.getenv("q")))
}
if(model == "treg"){
  beta <- na.omit(as.double(Sys.getenv(paste0("beta",0:10))))
  beta <- c(beta, rep(0,p-length(beta)+1))
  param <- list(beta = beta,
                sig2 = as.double(Sys.getenv("sig2")),
                nu = as.double(Sys.getenv("nu")))
}
if(model == "betareg"){
  beta <- na.omit(as.double(Sys.getenv(paste0("beta",0:10))))
  beta <- c(beta, rep(0,p-length(beta)+1))
  param <- list(beta = beta,
                phi = as.double(Sys.getenv("phi")))
}
k <- length(unlist(param))

MC <- 1e4 # Monte Carlo simulations
# seeds for random generator
set.seed(round(1e3 * param[[1]]))
random_number1 <- sample.int(1e6,1)
set.seed(7483 * n)
random_number2 <- sample.int(1e6,1)
set.seed(98323 + random_number1 + random_number2) # random hand typing algorithm
se1 <- sample.int(1e7, MC)
se2 <- sample.int(1e7, MC)
se3 <- sample.int(1e7, MC)
res <- list(
  implicit = matrix(nrow = MC, ncol = 2 * k),
  percentile = matrix(nrow = MC, ncol = 2 * k),
  boott = matrix(nrow = MC, ncol = 2 * k),
  boottrob = matrix(nrow = MC, ncol = 2 * k),
  bca = matrix(nrow = MC, ncol = 2 * k),
  classic = matrix(nrow = MC, ncol = 2 * k),
  classicrob = matrix(nrow = MC, ncol = 2 * k),
  time = matrix(nrow = MC, ncol = 5),
  initial = matrix(nrow = MC, ncol = k)
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
  x <- NULL
  if(model == "treg" || model == "betareg"){
    set.seed(se3[m])
    x <- cbind(1, matrix(rnorm(n*p),nc=p))
  }
  y <- rmodel(n, x, param, model, se1[m])

  # initial estimator
  fit_sample <- NULL
  try(fit_sample <- mle_estimate(y, x, unlist(param), model), silent=TRUE)
  if(is.null(fit_sample)) next
  if(any(is.na(fit_sample$par))) next
  theta_hat <- fit_sample$par
  res$initial[m,] <- theta_hat

  # JIME
  t1 <- Sys.time()
  res$implicit[m,] <- impboot_interval(y, x, theta_hat, model, control = list(seed = se2[m]))
  t2 <- Sys.time()
  res$time[m,1] <- difftime(t2, t1, units = "secs")

  # percentile bootstrap
  t1 <- Sys.time()
  boot <- parboot(y, x, theta_hat, model, control = list(seed = se2[m]))
  res$percentile[m,] <- parboot_interval(y, x, theta_hat, boot, model)
  t2 <- Sys.time()
  res$time[m,2] <- difftime(t2, t1, units = "secs")

  # studentized bootstrap (with second bootstrap)
  t1 <- Sys.time()
  res$boott[m,] <- parboott_interval(y, x, theta_hat, model, boot = boot, control = list(seed = se2[m]))
  t2 <- Sys.time()
  res$time[m,3] <- difftime(t2, t1, units = "secs")

  # studentized bootstrap (with robust variance)
  t1 <- Sys.time()
  res$boottrob[m,] <- parboott_interval(y, x, theta_hat, model, boot = boot, control = list(seed = se2[m], robust = TRUE))
  t2 <- Sys.time()
  res$time[m,4] <- difftime(t2, t1, units = "secs")

  # BCa bootstrap
  t1 <- Sys.time()
  res$bca[m,] <- BCa_interval(y, x, theta_hat, model, boot = boot, control = list(seed = se2[m]))
  t2 <- Sys.time()
  res$time[m,5] <- difftime(t2, t1, units = "secs")

  # classic
  res$classic[m,] <- asymptotic_interval(y, x, theta_hat, model, boot = boot, control = list(seed = se2[m]))
  res$classicrob[m,] <- asymptotic_interval(y, x, theta_hat, model, boot = boot, control = list(seed = se2[m], robust = TRUE))

  # save results
  save(res, file=paste0("tmp/",model,"_id_",id_slurm,".rds"))
  cat(m,"\n")
}

