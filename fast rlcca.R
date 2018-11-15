require(dqrng)
# Making the LCCA code faster
# Making the LCCA in only its support space ("acceptable paramters only")
# Removed threshold and nu
# Also I think that t0 is non decision time here and I removed it

rlcca=function(n.items, max.time, startx, drift, K, L, eta, dt, tau, t0){
  #   // Params:
  #   // -------
  #   // I: input (i.e., drift rate), vector of length n.items
  #   // x: accumulator, usually start as vector of length n.items
  #   // n.items: the number of competing accumulators
  #   // max_time: max number of iterations (max_time*dt is actual time)
  #   // K: scale of leak
  #   // L: scale of lateral inhibition
  #   // eta: scale of noise
  #   // thresh: decision threshold (I often use 1.0, but this can be fit)
  #   // dt: timestep, I often use 10
  #   // tau: time constant, I've used 1000, but this can be fit
  #   // t0: non-decision time
  #   // 
  #   // Returns:
  #   // --------
  #   // t: the iteration where it crossed (will be max_time if not crossed)
  #   //
  #   // You will have to look at the values in x to determine which
  #   // accumulator crossed.  These should be set because it's a vector
  #   // passed by value.
  
  # Exponentiate variables to bound them to positive values only
  K.0 = exp(K)
  L.0 = exp(L)
  
  # Upfront calculations adn storage for the model later
  dt_tau = dt/(tau)
  sqrt_dt_tau = sqrt(dt/(tau))
  
  lx=numeric(n.items)
  x=matrix(NA,max.time,n.items)
  
  x[1,]=startx
  
  for(t in 2:(max.time)){
    # get the sum for the lateral inhibition
    sumx = 0
    for(i in 1:(n.items)){
      sumx = x[t-1,i] + sumx
    }
    
    # loop over items
    for(i in 1:(n.items)){
      # calculate the lateral inhibition
      lx[i] = sumx - x[t-1,i]
      
      # determine the change in x
      temp=dqrnorm(n = 1, mean = 0,sd = eta)
      # This is for multinode ("Mutlichoice")
      # x[t,i] = x[t-1,i] +(as.numeric(drift[t,i]) - (K.0)*(x[t-1,i]) - (L.0)*lx[i] - nu/(n.items-1)*sum(drift[t,-i]))*(dt_tau) +(temp*sqrt_dt_tau)
      
      # This is for 2 node ("2Choice")
      x[t,i] = x[t-1,i] +
        (drift[t,i] - (K.0)*(x[t-1,i]) - (L.0)*lx[i])*(dt_tau) +
        (temp*sqrt_dt_tau)
      
      x <<- x
      # make sure not below zero
      if(!is.finite(x[t,i])){
        x=NaN
        return(x)
        }
      if(x[t,i] < 0){x[t,i] = 0}
      # if(!is.finite(x[t,i]) & (x[t,i]) < 0){x[t,i] = 0}
    }
  }
  # Return states of accumulators only
  return(x)  
}


