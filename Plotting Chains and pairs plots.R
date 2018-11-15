# load("C:/Users/Chief/Dropbox/Ulrich/Noah Data Giv9_Ulrich2015_prior_sdx3_full_posterior_with_purify_in_samples_mu_exp.RData")
################
# Make a matrix of Theta
# I will use this later to make scatterplots
mTheta=matrix(,ncol = 4, nrow = 24)
mTheta=Theta[,,1]
colnames(mTheta)=colnames(Theta)
row.names(mTheta)=paste("Chain", 1:24)
row.names(Theta)=row.names(mTheta)
for(i in 1:dim(Theta)[3]){
  mTheta=rbind(mTheta,Theta[,,i])
}

################
# Overview of Length of all samples
# Lets plot the densities over all samples
# They are increasing, which is expected
ts.plot(t(density), col=1:24, ylab="Log Density")

# The shape parameter of the gamma seems ill behaved still
# One chain (Neon blue here) spikes ouside the rest of the 
# chains
ts.plot(t(exp(Theta[,"A",])), col=1:24, ylab="A")
# here is "zoom in" and plot the true value in thick black
ts.plot(t(exp(Theta[,"A",])), col=1:24, ylim=c(0, 250), ylab="A")
abline(h = (exp(true["A"])), lw=3)
# it was interesting to plot the parameter on the log scale
# I do that here
ts.plot(t((Theta[,"A",])), col=1:24, ylab="Log A")
abline(h = ((true["A"])), lw=3)

# Here is a histogram of the parameter
# Note we are back to real scale of "A"
# In red is the true parameter
hist(t(exp(Theta[,"A",])), breaks = 50, xlim = c(0, 250), xlab="A")
abline(v = (exp(true["A"])), col="red")
# Again, I plot "A" on log scale
hist(t((Theta[,"A",])), breaks = 50, xlab="Log A")
abline(v = ((true["A"])), col="red")

# We spike quite high at one point
# Good thing this is in the burnin though
# I would want an MLE algorithm to do better however
ts.plot(t(exp(Theta[,"tau",])), col=1:24, ylab="tau")
# Again, we'll zoom in
ts.plot(t(exp(Theta[,"tau",])), col=1:24, ylim=c(0,500), ylab="tau")
abline(h = (exp(true["tau"])), lw=3)
# histograms is also zoomed in on the same window as the last plot
hist(t(exp(Theta[,"tau",])), breaks = 1e6, xlim = c(0, 500), xlab="tau")
abline(v = (exp(true["tau"])), col="red")

# Mu.c seemed fine at first glance
# Note the histogram tells a different story
ts.plot(t(exp(Theta[,"mu.c",])), col=1:24, ylab="mu.c")
abline(h = (exp(true["mu.c"])), lw=3)
# Histogram of mu.c
# We can see here that most of the estimates are zero
# this happens aftter sampling posteriors
hist(t(exp(Theta[,"mu.c",])), breaks = 50, xlab="mu.c")
abline(v = (exp(true["mu.c"])), col="red")

# B.1 converges but does not seem to capture the truth well
ts.plot(t(exp(Theta[,"b.1",])), col=1:24, ylab="b.1")
abline(h = (exp(true["b.1"])), lw=3)

hist(t(exp(Theta[,"b.1",])), breaks = 50, xlab="b.1")
abline(v = (exp(true["b.1"])), col="red")

# I am trying to work on a way to plot the pairs
# This doesn't look that great 
pairs(exp(mTheta), col=1:24)

# Let's dissect the chains into their two sections
# First the Burnin
# Later the Posterior Samples
################
# Init to Burnin
# Density looks fine
# Pretty cool that the yellow and light blue chains are ill performing
# I guess that is expected
ts.plot(t(density[,1:burnin]), col=1:24, ylab="Log Density")
# Chains generally will look sticky since I am not "purifying"
# I think this is fine since we are just looking for starting values
ts.plot(t(exp(Theta[,"A",1:burnin])), col=1:24, ylab="A")
abline(h = (exp(true["A"])), lw=3)
hist(t(exp(Theta[,"A",1:burnin])), breaks=50, xlab="A")
abline(v = (exp(true["A"])), col="red")

# It's so odd that the yellow chain can even be accepted...
ts.plot(t(exp(Theta[,"tau",1:burnin])), col=1:24, ylab="tau")
# "Zoom In"
ts.plot(t(exp(Theta[,"tau",1:burnin])), col=1:24, ylim=c(0,250), ylab="tau")
abline(h = (exp(true["tau"])), lw=3)
# Note we are still on the same "zoom in" scale
hist(t(exp(Theta[,"tau",1:burnin])), breaks = 1e6, xlim=c(0,250), xlab="tau")
abline(v = (exp(true["tau"])), col="red")

ts.plot(t(exp(Theta[,"mu.c",1:burnin])), col=1:24, ylab="mu.c")
abline(h = (exp(true["mu.c"])), lw=3)
hist(t(exp(Theta[,"mu.c",1:burnin])), breaks=50,xlim=c(0,3), ylab="mu.c")
abline(v = (exp(true["mu.c"])), col="red")

ts.plot(t(exp(Theta[,"b.1",1:burnin])), col=1:24, ylab="b.1")
abline(h = (exp(true["b.1"])), lw=3)
hist(t(exp(Theta[,"b.1",1:burnin])), breaks=50, xlab="b.1")
abline(v = (exp(true["b.1"])), col="red")

# Pairs plot looks better
# I imagine taking on the ill performing chains may help
pairs(exp(mTheta[1:(24*burnin),]), col=1:24)

# Va bene
pairs(exp(mTheta[(50*24):(24*burnin),]), col=1:24)

# Ok let's look at the posteriors now
################
# Burnin to Posterior
# Note in this density plot by stroke of luch light blue gets accepted
ts.plot(t(density[,burnin:N]), col=1:24, ylab="Density")
# Very outside our range 
ts.plot(t(exp(Theta[,"A",burnin:N])), col=1:24, ylab="A")
# Zoom In
# Seems to do decent but I am not sure that the truth is captured
ts.plot(t(exp(Theta[,"A",burnin:N])), col=1:24, ylim=c(0,200), ylab="A")
abline(h = (exp(true["A"])), lw=3)
# I am untrusing here
hist(t(exp(Theta[,"A",burnin:N])), breaks=50, xlim=c(0,200), xlab="A")
abline(v = (exp(true["A"])), col="red")

# Seems like correlation between tau and A?
ts.plot(t(exp(Theta[,"tau",burnin:N])), col=1:24, ylab="tau")
# Zoom in
ts.plot(t(exp(Theta[,"tau",burnin:N])), col=1:24, ylim=c(0,250), ylab="tau")
abline(h = (exp(true["tau"])), lw=3)
# I believe this a little more
hist(t(exp(Theta[,"tau",burnin:N])),  breaks=50, xlim=c(0,250), xlab="tau")
abline(v = (exp(true["tau"])), col="red")

# I guess Mu.c looks fine, I would hope for something tighter
ts.plot(t(exp(Theta[,"mu.c",burnin:N])), col=1:24, ylab="mu.c")
abline(h = (exp(true["mu.c"])), lw=3)
# I don't like that most of the chains "floored out"
hist(t(exp(Theta[,"mu.c",burnin:N])), xlab="mu.c")
abline(v = (exp(true["mu.c"])), col="red")

# looks like a complete miss to me
ts.plot(t(exp(Theta[,"b.1",burnin:N])), col=1:24, ylab="b.1")
abline(h = (exp(true["b.1"])), lw=3)
hist(t(exp(Theta[,"b.1",burnin:N])), xlab="b.1")
abline(v = (exp(true["b.1"])), col="red")

# I could definitely use some help interpeting this
# I guess that Mu.c and b.1 estimates mix the best
# but b.1 is quite off and the other estimates have outliers
pairs(exp(mTheta[(24*burnin):dim(mTheta)[1],]), col=1:24)

# TL;DR
# It looks like there may be some issues on parameter recovery with the prior from Ulrich's data, 
# However, this is only with 400 samples
# I was able to capture one set of parameters decently well without a prior 
# but this sample with a prior, maybe not
# I am skeptical that all the parameters of the gamma distribution can be fit
# I only fit 4 of the parameters like you do
# There are still more parameters that can be fit (starting point may be the only one)
# I think that this model may benefit from a good prior 
# rather than the one I stole for Ulrich's data

# Going forward:
# I suspect there is something to be learned from the above correalation plot,
# but I do not know exactly how to interpret it or ameliorate issues
# I would like to talk to you about some priors for the model, maybe we should consult Ulrich
# Maybe we could just run these chains longer?
# We are fitting this data hierarchically

# I am going to need help from you to implement the hierarchical model
# I understand the princple, basically updating the prior distribution, 
# but I just haven't built one from scratch

# Plotting Observed data
plot.data=obs.data
plot.data=ifelse(plot.data[,2]==1, plot.data,-plot.data)
hist(plot.data)

sim.data=array(NA, dim=c(24,5000,2))
for(i in 1:24){
  sim.data[i,,]=rmdc(N = 5000, A = true$A, tau = Theta[i,"tau",N], 
       mu.c = Theta[i,"mu.c",N], b.1 = Theta[i,"b.1",N])
}

plot.data=sim.data[5,,]
plot.data=ifelse(plot.data[,2]==1, plot.data,-plot.data)
hist(plot.data)


# You can skip this if you'd like
# I am trying to learn some fancy joint posterior plots
################
library(MASS)
# Joint Density plots
# These are pretty, but please note these are on the log scale
z=kde2d(x=as.numeric(Theta[,"A",burnin:N]),h=10000, y=as.numeric(Theta[,"tau",burnin:N]))
image(z, xlab = "A", ylab = "tau")
contour(z, add = TRUE, drawlabels = FALSE)
points(x=true$A, y=true$tau, pch = "X")

z=kde2d(x=as.numeric(Theta[,"A",burnin:N]),h=10, y=as.numeric(Theta[,"mu.c",burnin:N]))
image(z, xlab = "A", ylab = "mu.c")
contour(z, add = TRUE, drawlabels = FALSE)
points(x=true$A, y=true$mu.c, pch = "X")

z=kde2d(x=as.numeric(Theta[,"A",burnin:N]),h=10, y=as.numeric(Theta[,"b.1",burnin:N]))
image(z, xlab = "A", ylab = "b.1")
contour(z, add = TRUE, drawlabels = FALSE)
points(x=true$A, y=true$b.1, pch = "X")

z=kde2d(x=as.numeric(Theta[,"tau",burnin:N]),h=10, y=as.numeric(Theta[,"mu.c",burnin:N]))
image(z, xlab = "tau", ylab = "mu.c")
contour(z, add = TRUE, drawlabels = FALSE)
points(x=true$tau, y=true$mu.c, pch = "X")

z=kde2d(x=as.numeric(Theta[,"tau",burnin:N]),h=10, y=as.numeric(Theta[,"b.1",burnin:N]))
image(z, xlab = "tau", ylab = "b.1")
contour(z, add = TRUE, drawlabels = FALSE)
points(x=true$tau, y=true$b.1, pch = "X")

z=kde2d(x=as.numeric(Theta[,"mu.c",burnin:N]),h=2, y=as.numeric(Theta[,"b.1",burnin:N]))
image(z, xlab = "mu.c", ylab = "b.1")
contour(z, add = TRUE, drawlabels = FALSE)
points(x=true$mu.c, y=true$b.1, pch = "X")
