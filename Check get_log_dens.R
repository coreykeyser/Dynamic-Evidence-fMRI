# Exploring Density
library(tidyverse) # Library for quick data manipulation and exploration
library(doParallel)# For parallel multicore script execution
library(dqrng)     # For faster RNG
library(compiler)  # For compiling R code
source(('fast rlcca.R')) # Contains random Leaky Contious Competing
# Accumulator simulation process
source(('pseudolikelihood(trialByTrial).r')) # Contains likelihood and prior functions
source(("DE Try 2.R"))   # Contains DEMCMC Samplers


# Compile Function to machine code for faster execution
rlcca <- cmpfun(rlcca)
pseudolikelihood <- cmpfun(pseudolikelihood)
get_log_dens <- cmpfun(get_log_dens)

# Make a cluster for parallel execution
print("Cluster")
# 3 is the number of cores
# stopImplicitCluster()
# cl <- makePSOCKcluster(7)
# registerDoParallel(cl)

# Set seed for reproducibility
set.seed(123)

# Get data from experiment
data = read_csv(("revisedEvidenceAccumulationMaster.csv"))

# Get one subjects data
sub.name <- "Corey Responses"
sub.data <- select(data, Trial, Coherence, sub.name)
colnames(sub.data)[3] <- "response"

# Rename varaibles for later use
sub.data  <- sub.data %>%
  mutate(
    drift1 = ifelse(Coherence > 0, Coherence, 0),
    drift2 = ifelse(Coherence < 0,-Coherence, 0),
    Iter = seq_along(Coherence)
  )

# Extract drifts to use in paramter recovery
drift = cbind(sub.data$drift1, sub.data$drift2)
numtrials = dim(drift)[1]/60
# drift = drift[1:125, ] # I am only running over part of the drift data for
# computational ease (Parameter Recovery Only)

# true Params with the estimated parameters first (clear why below)
true = tibble(
  L = log(.4),
  # lateral inhibition
  K = log(.1),
  # leakage
  t0 = .0,
  # nondecision time
  xi = 1,
  # this is xi/sqrt(dt/tau)
  startx = c(1, 1),
  # starting point
  I0 = .1,
  # ?
  delta_t = .1,
  # change in time constant
  dt = .1,
  # step size
  maxtime = dim(drift)[1]/numtrials,
  # total amount of time before terminating
  eta = .1,
  # noise added to accumulation process (add more later)
  n.items = 2,
  # number of choices
  tau = 1000        # some scaling constant (?)
)
drift=drift+true$I0

# Get free parameter names and amount
params = param.names = names(true)[2]
n.params = length(param.names)
# Size of proposed dataset in pseudolikelihood
n.obs = 500

print("Make Data")
obs.data=array(NA, dim=c(n.obs,numtrials,60))
for(s in 1:n.obs){
  sub.data = foreach(
    i = 1:numtrials,
    .combine = rbind,
    .export = c("rlcca", "true", "drift"),
    .packages = c("dqrng")
  ) %dopar% {
    # Simulate one trial
    temp = rlcca(
      n.items = true$n.items[1],
      max.time = true$maxtime[1],
      startx = true$startx,
      drift = drift[(((i-1)*60)+1):(i*60),],
      K = true$K[1],
      L = true$L[1],
      eta = true$eta[1],
      dt = true$dt[1],
      tau = true$tau[1],
      t0 = true$t0[1]
    )
    # return difference
    temp[, 1] - temp[, 2]
  }
  obs.data[s,,]=sub.data
}
obs.data

get_log_dens=function(params, data){
  # Check that proposal is finite
  if(all(is.finite(params))){
    # Container for simulated data
    sim.data=matrix(NA, nrow = n.obs, ncol = length(obs.data))
    # Simulate datasets for Pseudolike in Parallel
    sim.data=NULL
    logprobtrial=numeric(numtrials)
    for(n in 1:numtrials){
      print(n)
      for(i in 1:n.obs){
        temp=rlcca(n.items = true$n.items[1], max.time = true$maxtime[1], startx = true$startx, 
                   drift = drift[(((n-1)*60)+1):(n*60),], K = params[,'K'], L = params[,'L'], eta = true$eta[1], dt = true$dt[1],
                   tau = true$tau[1], t0 = params['t0'])
        
        # return difference
        if(all(is.finite(temp))){
          sim.data=rbind(sim.data, temp[,1]-temp[,2])
        }else{return(-Inf)}
      }
      # print("Sims Done!")
      # Return pseudolike of observed data under simulated prosposed data
      logprobtrial[n]=(pseudolikelihood(data.mx=sim.data,
                                        obs.data=obs.data[,n,]))
      print(logprobtrial[n])
      
    }
    if(all(is.finite(logprobtrial))){return(sum(logprobtrial))}
    else{return(-Inf)}
  }else{return(-Inf) # If not finite, return worse density possible
  }
}

params=true[1,]
params=as.matrix(params)
params[,'K']
get_log_dens(params = as.matrix(params), data = obs.data)
