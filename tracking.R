
source("lca.R")

################################################################

# To run FFI, you set beta=kappa=0 and free nu.
# To run LCA, you set nu=0 and free beta and kappa.

# dynamic models
alpha=100              # threshold
beta=.4   						# lateral inhibition
kappa=.1							# leakage
nu=0                  # feed forward inhibition
xi=1                  # this is xi/sqrt(dt/tau)
tau=.15               # nondecision time
drift = drift          # drift rate
startx=c(1,1)         # starting point
I0=.1 

delta_t=.1            # change in time constant
dt=.1             # step size
maxtime=5           # total amount of time before terminating

################################################################ generate coherency change 

coherencies = c(-.35,-.25,-.15,-.05,0,.05,.15,.25,.35)

drift.v = numeric(maxtime*100)
ldv = length(drift.v)
cohs = rep(coherencies,ldv)

out = NA

drift.v = rand.sample(cohs)

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

pars=list(drift=drift, alpha=alpha, beta=beta, kappa=kappa, xi=xi, nu=nu,
          tau=tau, startx=startx, delta_t=delta_t, dt=dt, maxtime=maxtime, I0=.1)


################################################################ run it

#### increasing lateral inhibiton
temp1=rlca(n.items=2, startx=pars$startx, drift=pars$drift+pars$I0, max.time=pars$maxtime*100, K=pars$kappa, 
      L=.2, nu=pars$nu, eta=pars$xi, thresh=pars$alpha, dt=pars$dt, tau=pars$delta_t, t0=pars$tau)

temp2=rlca(n.items=2, startx=pars$startx, drift=pars$drift+pars$I0, max.time=pars$maxtime*100, K=pars$kappa, 
           L=.4, nu=pars$nu, eta=pars$xi, thresh=pars$alpha, dt=pars$dt, tau=pars$delta_t, t0=pars$tau)

temp3=rlca(n.items=2, startx=pars$startx, drift=pars$drift+pars$I0, max.time=pars$maxtime*100, K=pars$kappa, 
           L=.6, nu=pars$nu, eta=pars$xi, thresh=pars$alpha, dt=pars$dt, tau=pars$delta_t, t0=pars$tau)

pdf("latInh.pdf",12,8)
xs = seq(0,maxtime*100,length=length(drift[,1]))
par(mfrow=c(4,1))
plot(xs,drift[,1],type = 'l',lwd=2,ylab="Drift",xlab='Time')
lines(xs,drift[,2],col="red",lwd=2,lty=2)
matplot(temp1$state[1:temp1$time,],lwd=c(2,2),type="l",lty=c(1,2),xlab="Time",ylab="Evidence",ylim=c(0,(max(temp1$state)+2)),
        main="Lateral inhibition = .2")
matplot(temp2$state[1:temp2$time,],lwd=c(2,2),type="l",lty=c(1,2),xlab="Time",ylab="Evidence",ylim=c(0,(max(temp2$state)+2)),
        main="Lateral inhibition = .4")
matplot(temp3$state[1:temp3$time,],lwd=c(2,2),type="l",lty=c(1,2),xlab="Time",ylab="Evidence",ylim=c(0,(max(temp3$state)+2)),
        main="Lateral inhibition = .6")
dev.off()
################################################################ plot it

# increasing leakage
temp1=rlca(n.items=2, startx=pars$startx, drift=pars$drift+pars$I0, max.time=pars$maxtime*100, K=.1, 
           L=pars$beta, nu=pars$nu, eta=pars$xi, thresh=pars$alpha, dt=pars$dt, tau=pars$delta_t, t0=pars$tau)
temp2=rlca(n.items=2, startx=pars$startx, drift=pars$drift+pars$I0, max.time=pars$maxtime*100, K=.3, 
           L=pars$beta, nu=pars$nu, eta=pars$xi, thresh=pars$alpha, dt=pars$dt, tau=pars$delta_t, t0=pars$tau)
temp3=rlca(n.items=2, startx=pars$startx, drift=pars$drift+pars$I0, max.time=pars$maxtime*100, K=.5, 
           L=pars$beta, nu=pars$nu, eta=pars$xi, thresh=pars$alpha, dt=pars$dt, tau=pars$delta_t, t0=pars$tau)

pdf("leak.pdf",12,8)
xs = seq(0,maxtime*100,length=length(drift[,1]))
par(mfrow=c(4,1))
plot(xs,drift[,1],type = 'l',lwd=2,ylab="Drift",xlab='Time')
lines(xs,drift[,2],col="red",lwd=2,lty=2)
matplot(temp1$state[1:temp1$time,],lwd=c(2,2),type="l",lty=c(1,2),xlab="Time",ylab="Evidence",ylim=c(0,(max(temp1$state)+2)),
        main = "Leak = 0.1")
matplot(temp2$state[1:temp2$time,],lwd=c(2,2),type="l",lty=c(1,2),xlab="Time",ylab="Evidence",ylim=c(0,(max(temp2$state)+2)),
        main = "Leak = 0.3")
matplot(temp3$state[1:temp3$time,],lwd=c(2,2),type="l",lty=c(1,2),xlab="Time",ylab="Evidence",ylim=c(0,(max(temp3$state)+2)),
        main = "Leak = 0.5")
dev.off()

################################################################ run it

#### increasing lateral inhibiton , no leakage
temp1=rlca(n.items=2, startx=pars$startx, drift=pars$drift+pars$I0, max.time=pars$maxtime*100, K=pars$kappa, 
           L=pars$beta, nu=pars$nu, eta=pars$xi, thresh=pars$alpha, dt=pars$dt, tau=pars$delta_t, t0=pars$tau)
temp2=rlca(n.items=2, startx=pars$startx, drift=pars$drift+pars$I0, max.time=pars$maxtime*100, K=pars$kappa, 
           L=pars$beta, nu=pars$nu, eta=pars$xi, thresh=pars$alpha, dt=pars$dt, tau=pars$delta_t, t0=pars$tau)
temp3=rlca(n.items=2, startx=pars$startx, drift=pars$drift+pars$I0, max.time=pars$maxtime*100, K=pars$kappa, 
           L=pars$beta, nu=pars$nu, eta=pars$xi, thresh=pars$alpha, dt=pars$dt, tau=pars$delta_t, t0=pars$tau)

pdf("LatInhNoleak.pdf",12,8)
xs = seq(0,maxtime*100,length=length(drift[,1]))
par(mfrow=c(4,1))
plot(xs,drift[,1],type = 'l',lwd=2,ylab="Drift",xlab='Time')
lines(xs,drift[,2],col="red",lwd=2,lty=2)
matplot(temp1$state[1:temp$time,],lwd=c(2,2),type="l",lty=c(1,2),xlab="Time",ylab="Evidence",ylim=c(0,(max(temp$state)+2)))
matplot(temp2$state[1:temp$time,],lwd=c(2,2),type="l",lty=c(1,2),xlab="Time",ylab="Evidence",ylim=c(0,(max(temp$state)+2)))
matplot(temp3$state[1:temp$time,],lwd=c(2,2),type="l",lty=c(1,2),xlab="Time",ylab="Evidence",ylim=c(0,(max(temp$state)+2)))
dev.off()
 


