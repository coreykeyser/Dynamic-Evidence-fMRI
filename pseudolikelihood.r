library(doParallel)
# Compute the pseudolikelihood for two sets of observations (RTs and responses)
# for a two-choice task (obs.rt.1 and obs.rt.2, where the length of each set is
# the number of each responses.
# The likelihood is computed from an estimated density constructed from a
# simulated data set data.mx.
# The data set data.mx is an Nx2 matrix, where data.mx[,1] are the RTs and
# data.mx[,2] are the responses (1, -1).
pseudolikelihood <- function(data.mx, obs.rt.1, obs.rt.2, log=TRUE) {
    # Compute the number of each response
    n.1 <- min(0,sum(data.mx[,2]==1))
    n.2 <- min(0,sum(data.mx[,2]==-1))

    # Initialize the likelihood vectors to -750 (so the default likelihood
    # value is exp(-750) which is approximately zero.
    like.1 <- rep(-750,times=length(obs.rt.1))
    like.2 <- rep(-750,times=length(obs.rt.2))
    lp.1 <- lp.2 <- -750

    # If there are sufficient observations to construct a density estimate
    # (2 or more), compute the likelihoods for each response.
    if(n.1>1) {
        lp.1 <- log(n.1/(n.1+n.2)) # log the proportion of responses = 1
	# Compute the kernel density estimate
        dens.1 <- density(data.mx[data.mx[,2]==1,1],kernel="epanechnikov",
	                  from=0)
        # Convert it to a function
        func.1 <- approxfun(dens.1$x,dens.1$y,rule=2)
	# Use the function to compute the likelihood of the observed RTs
        like.1 <- log(func.1(obs.rt.1))
	# If there are any zeros, convert the log to -750
        like.1[like.1==-Inf] <- -750
        }

    # See notes for response = 1 above
    if(n.2>1) {
        lp.2 <- log(n.2/(n.1+n.2))
        dens.2 <- density(data.mx[data.mx[,2]==-1,1],kernel="epanechnikov",from=0)
        func.2 <- approxfun(dens.2$x,dens.2$y,rule=2)
        like.2 <- log(func.2(obs.rt.2))
        like.2[like.2==-Inf] <- -750
        }
    # Add likelihood components together.  Note the proportion terms that
    # rescale the conditional likelihoods to joints.
    like <- sum(like.1)+sum(like.2)+n.1*lp.1+n.2*lp.2

    # Convert from log back to likelihood if necessary and return
    if (log==FALSE) like <- exp(like)
    return(like)
}    

get_log_dens=function(params, data){
  simData=foreach(i = 1:(n.obs/10), .combine = rbind, .export = c("rdmc"), .packages = c("dqrng")) %dopar% {
   out=rdmc(10, A=as.numeric(params["A"]), tau=as.numeric(params["tau"]), mu.c=as.numeric(params["mu.c"]),
       sigma=1, b.1=as.numeric(params["b.1"]), t.max=1000)
   out
  }
  sum(pseudolikelihood(data.mx=simData, 
                       obs.rt.1=data$rt1, 
                       obs.rt.2=data$rt2))
}

samplePrior=function(){
  A=log(rnorm(1,3,1))
  tau=log(rnorm(1,4,1))
  mu.c=log(rnorm(1,.5,.1))
  b.1=log(rnorm(1,4,1))
  out=tibble(A, tau, mu.c, b.1)
}

prior=function(params){
  prior=max(dnorm(exp(params["A"]),3,1,log=TRUE),-750)+max(dnorm(exp(params["tau"]),4,1,log=TRUE),-750)+
    max(dnorm(exp(params["mu.c"]),.5,.1,log=TRUE),-750)+max(dnorm(exp(params["b.1"]),4,1,log=TRUE),-750)
}
