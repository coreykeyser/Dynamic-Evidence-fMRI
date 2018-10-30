library(doParallel)
# Compute the pseudolikelihood for difference of accumulators though time
# The likelihood is computed from an estimated density constructed from a
# simulated data set data.mx.
# The data set data.mx is an Nxtime matrix, where dimension 1 are the simulated
# data from this proposal of theta* and dimenstion 2 is the state of the difference
# of the accumulators at time i
pseudolikelihood <- function(data.mx, obs.data, log=TRUE) {
  # Initialize the likelihood vectors to -750 for length of time
  # (so the default likelihood value is exp(-750) which is 
  # approximately zero).
  like <- rep(-750,times=length(obs.data))
  
  # Compute the likelihood at each time point i
  for(i in 1:length(obs.data)){
    # Compute the kernel density estimate
    dens.1 <- density(data.mx[,i], kernel="epanechnikov",
                      from=0)
    # Convert it to a function
    func.1 <- approxfun(dens.1$x,dens.1$y,rule=2)
    # Use the function to compute the likelihood of the observed RTs
    like[i] <- log(func.1(obs.data[i]))
  }
  # If there are any zeros, convert the log to -750
  like[like==-Inf] <- -750
  
  # Add log likelihood components together.  Note the proportion terms that
  # rescale the conditional likelihoods to joints.
  like <- sum(like)
  
  # Convert from log back to likelihood if necessary and return
  if (log==FALSE) like <- exp(like)
  return(like)
}    

# Wrapper function to us in sampler
# I like to make mine so that it always takes 
# Theta* and data 
get_log_dens=function(params, data){
  # Check that proposal is finite
  if(all(is.finite(params))){
    # Container for simulated data
    sim.data=array(NA, dim=c(n.obs,length(obs.data)))
    # Simulate datasets for Pseudolike in Parallel
    sim.data=foreach(i=1:n.obs, .combine = rbind, 
                     .export = c("rlcca", "true","drift"), 
                     .packages = c("dqrng"))%dopar%{
                       # Simulate one trial
                       temp=rlcca(n.items = true$n.items[1], max.time = true$maxtime[1], startx = true$startx, 
                                  drift = drift, K = params['K'], L = params['L'], eta = true$eta[1], dt = true$dt[1],
                                  tau = true$tau[1], t0 = params['t0'])
                       # return difference
                       temp[,1]-temp[,2]
                     }
    # Return pseudolike of observed data under simulated prosposed data
    return(pseudolikelihood(data.mx=sim.data,
                         obs.data=obs.data))
  }else{return(-Inf) # If not finite, return worse density possible
    }
}

# Sample a "prior" distribution of the varaible in model
  # Essentially, this is how likely you think that parameters are
samplePrior=function(){
  # Use log tranform because these will later be eponentiated
  # This is fairly "tight" which makes the Initialization easier
  L=log(rnorm(1,.4,1))  # Leak
  K=log(rnorm(1,.1,1))  # Inhibition
  t0=log(rnorm(1,.2,.1))# NonDecision time
  out=tibble(L,K,t0)
}

# Compute density under a "prior" distribution of the varaible in model
# Essentially, this is how likely you think that parameters are
prior=function(params){
  # Note that these are the same distribution as above
  # the exponentiation is used to cancel the log transform to calculate
  # the density under the prior distribution
  prior=max(dnorm(exp(params["L"]),.4,1,log=TRUE),-750)+
    max(dnorm(exp(params["K"]),.1,1,log=TRUE),-750)+
    max(dnorm(exp(params["t0"]),.2,.1,log=TRUE),-750)
}
