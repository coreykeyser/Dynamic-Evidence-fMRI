require(dqrng)
# Function to simulate RTs and responses for the DMC in C
rdmc <- function(N=1,A=3,tau=3.9,a=2,mu.c=.4,sigma=3,b.1=3.9,dt=.1,t.max=500,cong=1) {
  #  print(c(A,tau,mu.c,b.1))
  A.0 <- cong*exp(A)
  tau.0 <- exp(tau)
  b.1.0 <- exp(b.1)
  b.2.0 <- -exp(b.1)
  mu.c.0 <- exp(mu.c)
  
  # Define time points at which the process jumps  
  t <- seq(dt,t.max,by=dt)
  # Compute the drift rate of the superimposed processes
  
  # mu <- A.0*exp(-t/tau.0)*(exp(1)*t/(a-1)/tau.0)^(a-1)*((a-1)/t-1/tau.0) + mu.c
  X = expectation.X=numeric(length(t))
  expectation.X <- mu.c.0 * t + A.0 * exp(-t/tau.0) * (t * exp(1) / ((a-1)*tau.0))^(a-1)
  # rands=matrix(dqrnorm(length(t)*N,0,sigma), nrow = length(t))
  step=sqrt(dt)
  
  # Initialize RT and RE (response) vectors
  RT <- RE <- vector()
  # Simulate N trials
  for (i in 1:N){
    # Compute jumps
    # dX <- mu*dt + sigma*sqrt(dt)*rnorm(length(t))
    
    # X <- cumsum(rnorm(length(t),0,sigma) * sqrt(dt)) + expectation.X
    X <- cumsum(dqrnorm(length(t),0,sigma) * step) + expectation.X
    # rands=dqrnorm(length(t),0,sigma)
    # dX=cumsum(rands * step)
    # X = dX + expectation.X
    
    # Compute state vector
    # X <<- cumsum(dX)
    # Determine RT by which boundary was traversed first
    # steps.1 <- detect_index(X,function(x,b) (x>b) == T, b=b.1.0)
    # steps.2 <- detect_index(X,function(x,b) (x<b) == T, b=b.2.0)
    steps.1 <- which(X>b.1.0)[1]
    steps.2 <- which(X<b.2.0)[1]
    steps <-min(c(steps.1,steps.2),na.rm=TRUE)
    # print(steps)
    if (steps==Inf) steps <- t.max
    RT[i] <- steps*dt
    # Determine response
    RE[i] <- sign(X[steps])
  }
  
  out=cbind(RT=RT,RE=RE)
    # Return an Nx2 matrix of RTs and responses
return(out)
}
