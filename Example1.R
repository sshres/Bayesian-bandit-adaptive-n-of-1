source("AGGBAN_run.R")
source("INDBAN_run.R")
source("B_BAN_run.R")
source("FRN_run.R")

#population parameters
set.seed(2356)
tru.b.mu = c(7,-1) #population mean
tru.b.sig = c(1,3) #standard deviation suggesting between-patient variation
tru.y.sig = 2 #within-patient deviance

J = 30 #number of patients
nperiod = 6 #number of periods (two per cycle)
nobv = 5 #number of observations per patient

#indiv parameters
b0 = rnorm(J,tru.b.mu[1],tru.b.sig[1])
b1 = rnorm(J,tru.b.mu[2],tru.b.sig[2])
b = cbind(b0,b1)
trup = c(tru.b.mu,tru.b.sig,tru.y.sig,b0,b1)
set.seed(NULL) #undo seed setting



designs = c("FRN","AGGBAN","INDBAN", "BBAN")
for (i in 1:4){
  if (i==1){
    out <- FRN_run(J, nperiod, nobv, b, trup)} #Fixed randomization design
  if (i==2){ 
    out <- AGGBAN_run(J, nperiod, nobv, b, trup)} #Aggregated Bayesian bandit adaptive design
  if (i==3){ 
    out <- INDBAN_run(J, nperiod, nobv, b, trup)} #Individual Bayesian bandit adaptive design
  if (i==4){
    out <- B_BAN_run(J, nperiod, nobv, b, trup)} #Batch updated Bayesian bandit adaptive design

  allocprob = out[1] #allocation probabilities through the adaptive periods
  post_samples = out[2] #final posterior draws
  hout = out[3] #final output from stan
  period = out[4] #study period identifier
  id = out[5] #patient identifier
  t = out[6] #treatment allocations 
  y = out[7] #observed outcomes
  save(out, file=paste("Example1_",designs[i],".RData",sep=""))
}





