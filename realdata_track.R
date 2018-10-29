
source("lca.R")

################################################################ import data

org.data = read.csv("evidenceAccumulationMaster.csv",header=T)

trials=unique(org.data[,1])
T=length(trials)

data=NULL
for(j in 1:T){
  temp=org.data[org.data[,1]==trials[j],]
  coherence = temp[,2]
  sub1 = temp[,3]
  sub2 = temp[,4]
  sub3 = temp[,5]
  sub4 = temp[,6]
  sub5 = temp[,7]
  data[[j]]=list('coh'=coherence, 'sub1'=sub1,'sub2'=sub2,'sub3'=sub3,'sub4'=sub4,'sub5'=sub5)
}

save(data,file="tracking_data.RData")

xs = seq(0,60,length=60)
plot(xs,data[[2]]$coh,lwd=4,type='l',ylim=c(-1,1))
lines(xs,data[[2]]$sub1,col='red',lty=2,lwd=4)
lines(xs,data[[2]]$sub2,col='blue',lty=2,lwd=4)
lines(xs,data[[2]]$sub3,col='green',lty=2,lwd=4)
lines(xs,data[[2]]$sub4,col='purple',lty=2,lwd=4)


################################################################ restructure data trial 1 

drift.v = data[[1]]$coh 
ldv = length(drift.v)

drift = matrix(NA,ldv,2)

for(i in 1:ldv){
  if(drift.v[i]>0){
    drift[i,1] = drift.v[i]
    drift[i,2] = 0
  } 
  if(drift.v[i]<0){
    drift[i,1] = 0
    drift[i,2] = drift.v[i]
  } 
  if(drift.v[i]==0){
    drift[i,1] = 0
    drift[i,2] = 0
  } 
  drift = abs(drift)
}

################################################################ pars

# To run FFI, you set beta=kappa=0 and free nu.
# To run LCA, you set nu=0 and free beta and kappa.

# dynamic models
alpha=100              # threshold
beta=.4   						# lateral inhibition
kappa=.1							# leakage
nu=0                  # feed forward inhibition
xi=1                  # this is xi/sqrt(dt/tau)
tau= .15               # nondecision time
drift = drift          # drift rate
startx=c(1,1)         # starting point
I0=.1 

delta_t=.1            # change in time constant
dt=1             # step size
maxtime=60           # total amount of time before terminating

pars=list(drift=drift, alpha=alpha, beta=beta, kappa=kappa, xi=xi, nu=nu,
          tau=tau, startx=startx, delta_t=delta_t, dt=dt, maxtime=maxtime, I0=.1)

temp1=rlca(n.items=2, startx=pars$startx, drift=pars$drift+pars$I0, max.time=pars$maxtime, K=pars$kappa, 
           L=pars$beta, nu=pars$nu, eta=pars$xi, thresh=pars$alpha, dt=pars$dt, tau=pars$delta_t, t0=pars$tau)

par(mfrow=c(2,1))
xs = seq(0,60,length=60)
plot(xs,drift[,1],type='l',lwd=4)
lines(xs,drift[,2],col='red',lty=2,lwd=4)
matplot(temp1$state[1:temp1$time,],lwd=c(2,2),type="l",lty=c(1,2),xlab="Time",ylab="Evidence",ylim=c(0,(max(temp1$state)+2)))
