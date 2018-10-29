library(abind)
library(doParallel)
library(tidyverse)
library(here)
library(compiler)
library(dqrng)
source(here("LCCA.R"))
rlcca <- cmpfun(rlcca)

cl <- makePSOCKcluster(3)
registerDoParallel(cl)

acomb2 <- function(...) abind(..., along=2)
acomb3 <- function(...) abind(..., along=3)
acomb4 <- function(...) abind(..., along=4)

# dynamic models
thresh=100              # threshold
L=.4   						# lateral inhibition
K=.1							# leakage
nu=0                  # feed forward inhibition
xi=1                  # this is xi/sqrt(dt/tau)
tau=.15               # nondecision time
# drift = drift          # drift rate
startx=c(1,1)         # starting point
I0=.1 

delta_t=.1            # change in time constant
dt=.1             # step size
maxtime=5           # total amount of time before terminating

n.data <- dim(data)[1]
n.items <- 2 # left and right
n.sim <- 10


data=read_csv(here("evidenceAccumulationMaster.csv"))

sub.name <- "Chris Responses"
sub.data <- select(data, Trial, Coherence, sub.name)
colnames(sub.data)[3] <- "response"

sub.data  <- sub.data %>%
  mutate(drift1 = ifelse(Coherence>0, Coherence, 0),
         drift2 = ifelse(Coherence<0, -Coherence, 0),
         Iter = seq_along(Coherence))

sse.lcca <- function(x, data){
  startx <- c(1, 1)
  n.data <- dim(data)[1]
  n.items <- 2 # left and right
  n.sim <- 10
  dt <- .1
  eta <- .000001
  nu <- 0
  xi=1
  I0=.1 
  delta_t=.1            # change in time constant
  
  sim.array <- array(NA, dim = c(n.sim, n.data, n.items))
  for(i in 1:n.sim){
    sim.array[i, , ] <- rlcca(n.items = 2, max.time = n.data, startx = startx,
                              drift = cbind(data$drift1, data$drift2), 
                              K = x["K"], L = x["L"], nu = nu, eta = eta, thresh = 1000, 
                              dt = dt, tau = x["tau"], t0 = .2)$state
  }
  mean.accum <- apply(sim.array, c(2, 3), mean)
  mean.accum[, 1] <- normalit(mean.accum[, 1])
  mean.accum[, 2] <- normalit(mean.accum[, 2])
  mean.accum.diff <- mean.accum[, 1] - mean.accum[, 2]
  response <- normalit(data$response)
  sum((mean.accum.diff - response)^2)
}


x <- NULL
x$K <- K
x$L <- L
x$tau <- tau
tmp1 <- optim(unlist(x), sse.tmp, data=sub.data, control=list("maxit"=1000))

# test --------------------------------------------------------------------
data=read_csv(here("evidenceAccumulationMaster.csv"))

sub.name <- "Chris Responses"
sub.data <- select(data, Trial, Coherence, sub.name)
colnames(sub.data)[3] <- "response"

sub.data  <- sub.data %>%
  mutate(drift1 = ifelse(Coherence>0, Coherence, 0),
         drift2 = ifelse(Coherence<0, -Coherence, 0),
         Iter = seq_along(Coherence))
# K <- unname(tmp1$par[1])
# L <- unname(tmp1$par[2])
# tau <- unname(tmp1$par[3])
K <- x$K
L <- x$L
tau <- x$tau
sim.array <- array(NA, dim = c(n.data, n.items, n.sim))
sim       <- array(NA, dim = c(n.data, n.items))
t0=Sys.time()
sim.array=foreach(iter=1:n.sim, .combine = 'acomb3', .packages = c("dqrng", "tidyverse"), .multicombine=F)%dopar%{
  sim <- rlcca(n.items = 2, max.time = n.data, startx = startx,
                            drift = bind_cols(list("drift1"=sub.data$drift1, "drift2"=sub.data$drift2)),
                            K = K, L = L, nu = nu, eta = .000001, thresh = 1000, dt = dt, tau = tau, t0 = .2)$state
}
t1=Sys.time()
t1-t0

mean.accum <- apply(sim.array, c(2, 3), mean)
colnames(mean.accum) <- c("mean_accum1", "mean_accum2")
sub.data <- cbind(sub.data, mean.accum)
sub.data <- sub.data %>% mutate(Response_norm = normalit(response),
                                mean_accum1_norm = normalit(mean_accum1),
                                mean_accum2_norm = normalit(mean_accum2))

ggplot(sub.data) +
  geom_line(aes(Iter, mean_accum1_norm - mean_accum2_norm), col = "red") +
  geom_line(aes(Iter, 2*(Response_norm-.5)), col = "black") +
  ylab("Accumulators")


