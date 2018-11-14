library(tidyverse) # Library for quick data manipulation and exploration
library(doParallel)# For parallel multicore script execution
library(dqrng)     # For faster RNG
library(compiler)  # For compiling R code
source(('fast rlcca.R')) # Contains random Leaky Contious Competing
# Accumulator simulation process
source(('pseudolikelihood.r')) # Contains likelihood and prior functions
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
data = read_csv(("evidenceAccumulationMaster.csv"))

# Get one subjects data
sub.name <- "Chris Responses"
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
drift = drift[1:60, ] # I am only running over part of the drift data for computational ease (Parameter Recovery Only)

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
  maxtime = dim(drift)[1],
  # total amount of time before terminating
  eta = .1,
  # noise added to accumulation process (add more later)
  n.items = 2,
  # number of choices
  tau = 1000        # some scaling constant (?)
)
drift=drift+true$I0

# Get free parameter names and amount
params = param.names = names(true)[1:2]
n.params = length(param.names)
# Size of proposed dataset in pseudolikelihood
n.obs = 5000

print("Make Data")
obs.data = foreach(
  i = 1:n.obs,
  .combine = rbind,
  .export = c("rlcca", "true", "drift"),
  .packages = c("dqrng")
) %dopar% {
  # Simulate one trial
  temp = rlcca(
    n.items = true$n.items[1],
    max.time = true$maxtime[1],
    startx = true$startx,
    drift = drift,
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

# Hyper
# Number of chains
K = n.params * 2 + 1
# Length of burnin process
burnin = 100
# Length of all samples
N = 500

# containers
density = matrix(NA, K, N)
Theta = array(NA, dim = c(K, n.params, N))
dimnames(Theta)[2] = list(param.names)

# Initialize Thetas
############
print("Init")
t0 = Sys.time()
log.dens.like = get_log_dens # get likelihood function
out = initializeChain(Theta, K, density, log.dens.like, data, param.names)
Theta = out$Theta
density = out$density
t1 = Sys.time()
InitTime = t1 - t0

# burnin
############
print("Burnin")
t0 = Sys.time()
mig_prob = .05 # Probability of theta migration at each burnin sample
out = BURNIN(
  n = 2,
  N = burnin,
  n.params,
  K,
  log.dens.like,
  data,
  prior,
  mig_prob,
  param.names
)
log.dens.like = get_log_dens # get likelihood function
Theta = out$Theta
density = out$density
burnTheta = Theta
burndensity = density
t1 = Sys.time()
burnTime = t1 - t0

# Vanilla
##########
# Draws posterior samples Turner 2014 (?)
print("Vanilla")
t0 = Sys.time()
out = DEMCMC(n = burnin, N, n.params, K, log.dens.like, data, prior)
# Theta = out$Theta
# density = out$density
vanillaTheta = Theta
vanillaburndensity = density
t1 = Sys.time()
vanillaTime = t1 - t0

# Save file for laters
save.image("Posterior Sample Data.RData")
