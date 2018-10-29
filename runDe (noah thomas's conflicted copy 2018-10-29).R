library(tidyverse)
library(here())
library(doParallel)
library(dqrng)
library(compiler)
source(here('fast rdmc.R'))
source(here('pseudolikelihood.r'))
source(here('get.data.r'))
source(here("DE Try 2.R"))

rmdc <- cmpfun(rdmc)


print("Cluster")
cl <- makePSOCKcluster(3)
registerDoParallel(cl)
set.seed(123)

true=tibble(A = log(3), tau = log(4), mu.c = log(.5), b.1 = log(4))
param.names=names(true)
n.params=length(param.names)

print("Make Data")
obs.data=rdmc(N = 1000, A = log(3), tau = log(4), mu.c = log(.5), b.1 = log(4), sigma=1)
obs.rt.1 <- obs.data[obs.data[,2]==1,1]
obs.rt.2 <- obs.data[obs.data[,2]==-1,1]
data=list(rt1=obs.rt.1,
            rt2=obs.rt.2)

# Hyper
K=n.params*2+1
burnin=1000*.1
N=5000*.1
n.obs = 500

# containers
density=matrix(NA, K, N)
Theta=array(NA, dim=c(K,n.params, N))
dimnames(Theta)[2]=list(names(true))

# Init 
############
print("Init")
t0=Sys.time()
log.dens.like=get_log_dens
out=initializeChain(Theta, K, density, log.dens.like, data, param.names)
Theta=out$Theta
density=out$density
t1=Sys.time()
InitTime=t1-t0

# burnin 
############
print("Burnin")
t0=Sys.time()
mig_prob=.5
params=(names(true))
out=BURNIN(burnin, n.params, K, Theta, log.dens.like, data, prior, mig_prob, params)
Theta=out$Theta
density=out$density
burnTheta=Theta
burndensity=density
t1=Sys.time()
burnTime=t1-t0


# Vanilla
##########
print("Vanilla")
t0=Sys.time()
out=DEMCMC(n=burnin, N, n.params, K, Theta, log.dens.like, data, prior)
Theta=out$Theta
density=out$density
vanillaTheta=Theta
vanillaburndensity=density
t1=Sys.time()
vanillaTime=t1-t0

# Crossover Probabilitically
###########
print("Crossover")
t0=Sys.time()
Theta=burnTheta
density=burndensity
CR=.9
out=DEMCMCcrossover(n=burnin, N, n.params, K, Theta, log.dens.like, data, prior, CR)
crossoverTheta=out$Theta
crossoverdensity=out$density
t1=Sys.time()
crossTime=t1-t0

# thisTheta=vanillaTheta
# ts.plot(t(exp(thisTheta[,"A",])))
# ts.plot(t(exp(thisTheta[,"b.1",])))
# ts.plot(t(exp(thisTheta[,"tau",])))
# ts.plot(t(exp(thisTheta[,"mu.c",])))

timeTB=tibble(InitTime, burnTime, vanillaTime, crossTime)
print(timeTB)

save.image("Noah Data Giv2.RData")
