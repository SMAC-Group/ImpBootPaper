# -----------
# Simulations
# -----------
# general setting

# simulation specifics
MC <- 10000 # number of simulations

model <- as.character(Sys.getenv("MODEL"))
n <- as.integer(Sys.getenv("N")) # sample size
param <- as.double(Sys.getenv("PARAM"))
p <- 2

implicit = matrix(nrow = MC, ncol = 2 * p)
percentile = matrix(nrow = MC, ncol = 2 * p)
boott = matrix(nrow = MC, ncol = 2 * p)
boottrob = matrix(nrow = MC, ncol = 2 * p)
bca = matrix(nrow = MC, ncol = 2 * p)
classic = matrix(nrow = MC, ncol = 2 * p)
classicrob = matrix(nrow = MC, ncol = 2 * p)
time = matrix(nrow = MC, ncol = 5)
initial = matrix(nrow = MC, ncol = p)
implicit_exact = matrix(nrow = MC, ncol = p)
percentile_exact = matrix(nrow = MC, ncol = p)

##------------------ Slurm specs --------------
n_array <- 1000
ind <- matrix(seq_len(MC),nr=n_array,byr=T)

for(i in seq_len(n_array)){
  if(!file.exists(file=paste0("tmp/",model,"_id_",i,".rds"))) next
  load(file=paste0("tmp/",model,"_id_",i,".rds"))
  implicit[ind[i,],] <- res$implicit[ind[i,],]
  percentile[ind[i,],] <- res$percentile[ind[i,],]
  boott[ind[i,],] <- res$boott[ind[i,],]
  boottrob[ind[i,],] <- res$boottrob[ind[i,],]
  bca[ind[i,],] <- res$bca[ind[i,],]
  classic[ind[i,],] <- res$classic[ind[i,],]
  classicrob[ind[i,],] <- res$classicrob[ind[i,],]
  time[ind[i,],] <- res$time[ind[i,],]
  initial[ind[i,],] <- res$initial[ind[i,],]
  implicit_exact[ind[i,],] <- res$implicit_exact[ind[i,],]
  percentile_exact[ind[i,],] <- res$percentile_exact[ind[i,],]
}

res <- list(
  implicit = implicit,
  percentile = percentile,
  boott = boott,
  boottrob = boottrob,
  bca = bca,
  classic = classic,
  classicrob = classicrob,
  time = time,
  initial = initial,
  implicit_exact = implicit_exact,
  percentile_exact = percentile_exact
)

save(res, file=paste0("data/",model,"_n_",n,"_param_",param,".rds"))
