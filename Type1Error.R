#####Type 1 error simulations####

source("BAN_run.R")
source("B_BAN_run.R")
source("FRN_run.R")

library(bayestestR)

#population parameters
set.seed(2356)
tru.b.mu = c(7,0) #lower y-value implies lower pain level so treatment effect is significant
tru.b.sig = c(1,0) #standard deviation suggesting between-patient variation
tru.y.sig = 2 #within-patient deviance

#indiv parameters
b0 = rnorm(J,tru.b.mu[1],tru.b.sig[1])
b1 = rnorm(J,tru.b.mu[2],tru.b.sig[2])
b = cbind(b0,b1)
trup = c(tru.b.mu,tru.b.sig,tru.y.sig,b0,b1)
set.seed(NULL) #undo seed setting

J = 30 #number of patients
nperiod = 6 #number of periods (two per cycle)
nobv = 1 #number of observations per patient

designs = c("FRN","BAN","BBAN")
sims = 1000
h0 = matrix(0,length(designs),sims)

for (j in 1:sims){
  
  for (i in 1:3){
    if (i==1){
      out <- FRN_run(J, nperiod, nobv, b, trup)} #Fixed randomization design
    if (i==2){ 
      out <- BAN_run(J, nperiod, nobv, b, trup)} #Bayesian bandit adaptive design
    if (i==3){
      out <- B_BAN_run(J, nperiod, nobv, b, trup)} #Batch updated Bayesian bandit adaptive design
    
    allocprob = out[1] #allocation probabilities through the adaptive periods
    post_samples = out[2] #final posterior draws
    hout = out[3] #final output from stan
    period = out[4] #study period identifier
    id = out[5] #patient identifier
    t = out[6] #treatment allocations 
    y = out[7] #observed outcomes
    save(out, file=paste("OP1_",designs[i],"_run",j,".RData",sep=""))
    
    q=post_samples$beta_mu[,2]
    ci_hdi <- ci(q, method = "HDI",ci=0.95) #posterior credible interval
    aa = ci_hdi$CI_low
    bb = ci_hdi$CI_high
    if (0>aa & 0<bb){h0[i,j]=0} #Decision rule based on CI
    else {h0[i,j]=1} #type I error
  }
}

type1error = rowSums(h0)/sims 
rbind(designs,type1error) 
