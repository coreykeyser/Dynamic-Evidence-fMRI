migration=function(current_chains,density,log.dens,data, params){
  pars=dim(current_chains)[2]
  lnum1=sample(c(1:K),1)										# determine how many groups to work with
  lnum2=sort(sample(c(1:K),lnum1,replace=F))							# which groups we will work with
  thetaset=matrix(NA,lnum1,pars)									# initialize
  dimnames(thetaset)[2]=list(params)
  currentset=propset=propw=currw=numeric(lnum1)
  index=numeric(lnum1)
  for(i in 1:lnum1){
    index[i]=sample(1:K,1,replace=F)
    thetaset[i,]=current_chains[lnum2[i],] + runif(1,-.01,.01)				# create a set of these particles to swap
    # print(thetaset)
    propset[i]=log.dens(thetaset[i,],data)
    # print(density)
    # print(lnum2[i])
    currentset[i]=density[lnum2[i]]
    propw[i]=propset[i]
    currw[i]=currentset[i]
  }
  if(runif(1) < exp(propw[lnum1] - currw[1])){
    current_chains[lnum2[1],]=thetaset[lnum1,]							# swap the first with the last (creating a circle)
    density[lnum2[1]]=propset[lnum1]
  }
  if(lnum1!=1){											# make sure we are not done yet
    for(i in 1:(lnum1-1)){		
      if(runif(1) < exp(propw[i] - currw[i+1])){
        current_chains[lnum2[i+1],]=thetaset[i,]							# swap the first with the last (creating a circle)
        density[lnum2[i+1]]=propset[i]
      }}}
  list(weight=density,theta=current_chains)
}

initializeChain = function(Theta, K, density, log.dens.like, data, params){
  for(k in 1:K){
    density[k,1]=-Inf
    while(!is.finite(density[k,1])){
      Theta[k,params,1]=as.numeric(samplePrior()[params])+runif(1)
      density[k,1]=log.dens.like(params = Theta[k,params,1], data)
    }
  }
  list("Theta"=Theta, "density"=density)
}  

BURNIN=function(N, n.params, K, Theta, log.dens.like, data, prior, mig_prob, params){
  for(i in 2:N){
    # Specify gamma
    # gamma=runif(1,.5,1)
    if(mig_prob > runif(1)){
      mig=migration(Theta[,,i-1], density[,i-1], log.dens = log.dens.like, data, params)
      density[,i-1]=mig$weight
      Theta[,,i-1]=mig$theta
    }
    gamma=2.38/sqrt(2*n.params)
    if(i%%10==0){gamma=1}
    best=which.max(density[,i-1])[1]
    for(k in 1:K){
      # sample theta m and n
      mn=sample(setdiff(1:K,k),2,F)
      # sample Eps
      eps=runif(1, -.01, .01)
      # propose theta.star
      theta.star = Theta[k,,i-1] + 
        gamma*(Theta[mn[1],,i-1]-
                 Theta[mn[2],,i-1])+gamma*(Theta[k,,i-1]-Theta[best,,i-1])
        eps
      # sample alpha
      alpha = runif(1)
      dens1 = density[k,i-1]
      # dens1 = log.dens.like(Theta[k,,i-1],  
      #                       data)+
      #   prior(Theta[k,,i-1])
      # print(theta.star)
      dens2 = log.dens.like(theta.star,  
                            data)+
        prior(theta.star)
      # MH step 
      if(alpha < exp(dens2 - dens1)){
        # if accept
        # store theta_k,i <- theta*
        # print("accept")
        Theta[k,,i] <- theta.star
        density[k,i] <- dens2
      }else{
        # print("reject")
        Theta[k,,i] <- Theta[k,,i-1]      
        density[k,i] <- dens1
      }
    }
    print(i)
    # print(Theta[,,])
  }
  list("Theta"=Theta, "density"=density)
}


DEMCMC=function(n, N, n.params, K, Theta, log.dens.like, data, prior){
  accept=0
  for(i in n:N){
    # Specify gamma
    # gamma=runif(1,.5,1)
    gamma=2.38/sqrt(2*n.params)
    if(i%%10==0){gamma=1}
    for(k in 1:K){
      # sample theta m and n
      mn=sample(setdiff(1:K,k),2,F)
      # sample Eps
      eps=runif(1, -.01, .01)
      # propose theta.star
      theta.star = Theta[k,,i-1] + 
        gamma*(Theta[mn[1],,i-1]-
                 Theta[mn[2],,i-1])+
        eps
      # sample alpha
      alpha = runif(1)
      dens1 = log.dens.like(Theta[k,,i-1],  
                            data)+
        prior(Theta[k,,i-1])
      
      dens2 = log.dens.like(theta.star,  
                            data)+
        prior(theta.star)
      # MH step 
      if(alpha < exp(dens2 - dens1)){
        # if accept
        # store theta_k,i <- theta*
        Theta[k,,i] <- theta.star
        density[k,i] <- dens2
        accept=accept+1
      }else{
        Theta[k,,i] <- Theta[k,,i-1]      
        density[k,i] <- dens1
      }
    }
    print(i)
    # print(Theta[k,,i])
  }
  list("Theta"=Theta, "density"=density, "accept"=accept)
}

DEMCMCblocked=function(n, N, n.params, K, Theta, log.dens.like, data, prior, block1, block2){
  for(i in n:N){
    # Specify gamma
    gamma=runif(1,.5,1)
    # gamma=2.38/sqrt(2*3)
    for(k in 1:K){
      # sample theta m and n
      mn=sample(setdiff(1:K,k),2,F)
      # sample Eps
      eps=runif(1, -.01, .01)
      # propose theta.star
      theta.star = Theta[k,,i-1]
      theta.star[block1] = Theta[k,block1,i-1]+ 
        gamma*(Theta[mn[1],block1,i-1]-
                 Theta[mn[2],block1,i-1])+
        eps
      # sample alpha
      alpha = runif(1)
      dens1 = log.dens.like(Theta[k,,i-1],  
                            data)+
        prior(Theta[k,,i-1])
      
      dens2 = log.dens.like(theta.star,  
                            data)+
        prior(theta.star)
      # MH step 
      if(alpha < exp(dens2 - dens1)){
        # if accept
        # store theta_k,i <- theta*
        Theta[k,block1,i] <- theta.star[block1]
        # density[k,i] <- dens2
      }else{
        Theta[k,block1,i] <- Theta[k,block1,i-1]      
        # density[k,i] <- dens1
      }
      # sample theta m and n
      mn=sample(setdiff(1:K,k),2,F)
      # sample Eps
      eps=runif(1, -.01, .01)
      # propose theta.star
      theta.star = Theta[k,,i-1]
      theta.star[block2] = Theta[k,block2,i-1]+ 
        gamma*(Theta[mn[1],block2,i-1]-
                 Theta[mn[2],block2,i-1])+
        eps
      # sample alpha
      alpha = runif(1)
      dens1 = log.dens.like(Theta[k,,i-1],  
                            data)+
        prior(Theta[k,,i-1])
      
      dens2 = log.dens.like(theta.star,  
                            data)+
        prior(theta.star)
      # MH step 
      if(alpha < exp(dens2 - dens1)){
        # if accept
        # store theta_k,i <- theta*
        Theta[k,block2,i] <- theta.star[block2]
        # density[k,i] <- dens2
      }else{
        Theta[k,block2,i] <- Theta[k,block2,i-1]      
        # density[k,i] <- dens1
      }
      density[k,i] <- log.dens.like(Theta[k,,i-1],  
                                            data)+
        prior(Theta[k,,i-1])
    }
    # print(i)
    # print(Theta[k,,i])
  }
  list("Theta"=Theta, "density"=density)
}

DEMCMCcrossover=function(n, N, n.params, K, Theta, log.dens.like, data, prior, CR){
  accept=0
  for(i in n:N){
    # Specify gamma
    # gamma=runif(1,.5,1)
    gamma=2.38/sqrt(2*n.params)
    if(i%%10==0){gamma=1}
    for(k in 1:K){
      # sample theta m and n
      mn=sample(setdiff(1:K,k),2,F)
      # sample Eps
      eps=runif(1, -.01, .01)
      # propose theta.star
      theta.star = Theta[k,,i-1]+ 
        gamma*(Theta[mn[1],,i-1]-
                 Theta[mn[2],,i-1])+
        eps
      # Crossover step
      mut=rbinom(5, size=1, (1-CR))
      theta.star[mut]=Theta[k,mut,i-1]
      # sample alpha
      alpha = runif(1)
      dens1 = log.dens.like(Theta[k,,i-1],  
                            data)+
        prior(Theta[k,,i-1])
      
      dens2 = log.dens.like(theta.star,  
                            data)+
        prior(theta.star)
      # MH step 
      if(alpha < exp(dens2 - dens1)){
        # if accept
        # store theta_k,i <- theta*
        Theta[k,,i] <- theta.star
        density[k,i] <- dens2
        accept=accept+1
      }else{
        Theta[k,,i] <- Theta[k,,i-1]      
        density[k,i] <- dens1
      }
      # density[k,i] <- log.dens.like(Theta[k,,i-1],  
      #                               data)+
      #   prior(Theta[k,,i-1])
    }
    print(i)
    # print(Theta[k,,i])
  }
  list("Theta"=Theta, "density"=density, "accept"=accept)
}