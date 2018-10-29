# LCA Continuous Draft New 10/22
library(tidyverse)
library(here)
source(here("LCCA.R"))

# Real Data
#######
data=read_csv(here("evidenceAccumulationMaster.csv"))

colnames(data)

chrisTB=select(data, Trial, Coherence, `Chris Responses`)

chrisTB=chrisTB %>%
  mutate(drift1 = ifelse(Coherence>0, Coherence, 0),
         drift2 = ifelse(Coherence<0, -Coherence, 0),
         Iter = seq_along(Coherence))

# Play with rlcca
#######
source(here("LCCA.R"))

# To run FFI, you set beta=kappa=0 and free nu.
# To run LCA, you set nu=0 and free beta and kappa.

# dynamic models
thresh=100              # threshold
L=.4   						# lateral inhibition
K=.1							# leakage
nu=0                  # feed forward inhibition
xi=1                  # this is xi/sqrt(dt/tau)
tau=.15               # nondecision time
drift = drift          # drift rate
startx=c(1,1)         # starting point
I0=.1 

delta_t=.1            # change in time constant
dt=.1             # step size
maxtime=5           # total amount of time before terminating

lcaArr=array(NA, dim = c(10, 1500, 2), dimnames = list(paste("sim",1:10,sep=''),
                                                         paste("time",1:1500,sep=''),
                                                         paste("drift",1:2,sep='')))
for(i in 1:10){
  lcaArr[i,,]=rlcca(n.items = 2, max.time = length(chrisTB$drift1), startx = startx, 
            drift = bind_cols(list("drift1"=chrisTB$drift1, "drift2"=chrisTB$drift2)), 
           K = K, L = L, nu = nu, eta = .0000001, thresh = 1000, dt = dt, tau = tau, t0 = .2)$state
}

# dimension reduction
lcaTB=reshape2::melt(lcaArr)
TB=as.tibble(lcaTB)

stb=TB%>%
  filter(Var1 == "sim1", Var3=="drift1")

unique(stb$Var3)

ggplot(stb)+geom_line(aes(Var2, value))

colnames(out$state)=list("ACCUM_1", "ACCUM_2")
accumTB=as.tibble(out$state)
accumTB=accumTB %>%
  mutate(Iter = seq_along(ACCUM_1), 
         ACCUM_1_norm = normalit(ACCUM_1),
         ACCUM_2_norm = normalit(ACCUM_2))

masterTB=bind_cols(accumTB,chrisTB)
masterTB=masterTB %>%
  mutate(Response_norm = normalit(`Chris Responses`))


ggplot(masterTB)+geom_line(aes(Iter, ACCUM_1_norm), col="red")+
  geom_line(aes(Iter, -ACCUM_2_norm), col="blue")+ylab("Accumulators")#+
  geom_line(aes(Iter, 2*(Response_norm-.5)), col="black")
