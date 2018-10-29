

rlca_v1=function(n.items, max.time, startx, drift, K, L, nu, eta, 
                 thresh, dt, tau, t0){
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
  #   // 
  #   // Returns:
  #   // --------
  #   // t: the iteration where it crossed (will be max_time if not crossed)
  #   //
  #   // You will have to look at the values in x to determine which
  #   // accumulator crossed.  These should be set because it's a vector
  #   // passed by value.
  dt_tau = dt/(tau)
  sqrt_dt_tau = sqrt(dt/(tau))
  lx=numeric(n.items)
  crossed=0
  
  x=matrix(NA,max.time,n.items)
  x[1,]=startx
  
  for(t in 2:(max.time)){
    
    # loop over items
    for(i in 1:(n.items)){
      
      # determine the change in x
      temp=rnorm(1,0,eta)
      x[t,i] = x[t-1,i] + (drift[i] - (K)*(x[t-1,i]) - (L)*sum(x[t-1,-i]) - nu/(n.items-1)*sum(drift[-i]) )*(dt_tau) + (temp*sqrt_dt_tau)
      
      # make sure not below zero
      if((x[t,i]) < 0)x[t,i] = 0
      
      # see if crossed thresh
      if((x[t,i]) >= (thresh)){
        # crossed, but we want to make sure to calculate all vals
        crossed = 1
      }
    }
    
    if(crossed == 1)break;
  }
  
  list("time"= t + t0, "state"=x)  
}



rlcca=function(n.items, max.time, startx, drift, K, L, nu, eta, thresh, dt, tau, t0){
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
  #   // 
  #   // Returns:
  #   // --------
  #   // t: the iteration where it crossed (will be max_time if not crossed)
  #   //
  #   // You will have to look at the values in x to determine which
  #   // accumulator crossed.  These should be set because it's a vector
  #   // passed by value.
  
  dt_tau = dt/(tau)
  sqrt_dt_tau = sqrt(dt/(tau))
  lx=numeric(n.items)
  crossed=0
  
  x=matrix(NA,max.time,n.items)
  x[1,]=startx
  
  for(t in 2:(max.time)){
    # get the sum for the lateral inhibition
    sumx = 0;
    for(i in 1:(n.items)){
      sumx = x[t-1,i] + sumx
    }
    
    # loop over items
    for(i in 1:(n.items)){
      # calculate the lateral inhibition
      lx[i] = sumx - x[t-1,i];
      
      # determine the change in x
      temp=rnorm(1,0,eta)
      x[t,i] = x[t-1,i] + (as.numeric(drift[t,i]) - (K)*(x[t-1,i]) - (L)*lx[i] - nu/(n.items-1)*sum(drift[-i]) )*(dt_tau) + (temp*sqrt_dt_tau)
      
      # make sure not below zero
      if((x[t,i]) < 0)x[t,i] = 0
      
      # see if crossed thresh
      if((x[t,i]) >= (thresh)){
        # crossed, but we want to make sure to calculate all vals
        crossed = 1
      }
    }
    
    if(crossed == 1)break;
  }
  
  list("time"= t + t0, "state"=x)  
}

normalit<-function(m){
  (m - min(m))/(max(m)-min(m))
}

