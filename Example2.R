source("BAN_run.R")
source("B_BAN_run.R")
source("FRN_run.R")

#population parameters
set.seed(90)
tru.b.mu = c(7,0) #population mean
tru.b.sig = c(0.5,0.3) #standard deviation suggesting between-patient variation
tru.y.sig = 2 #within-patient deviance

#indiv parameters
b0 = rnorm(J,tru.b.mu[1],tru.b.sig[1])
b1 = rnorm(J,tru.b.mu[2],tru.b.sig[2])
b = cbind(b0,b1)
b[19,2]=-1.5
b[4,2]=-2
b[25,2]=-1.7
b[11,2]=-1.85
set.seed(NULL) #undo seed setting

J = 30 #number of patients
nperiod = 6 #number of periods (two per cycle)
nobv = 5 #number of observations per patient

designs = c("FRN","BAN","BBAN")
for (i in 1:3){
  if (i==1){
    out <- FRN_run(J, nperiod, nobv, b, truep)} #Fixed randomization design
  if (i==2){ 
    out <- BAN_run(J, nperiod, nobv, b, truep)} #Bayesian bandit adaptive design
  if (i==3){
    out <- B_BAN_run(J, nperiod, nobv, b, truep)} #Batch updated Bayesian bandit adaptive design
  
  allocprob = out[1] #allocation probabilities through the adaptive periods
  post_samples = out[2] #final posterior draws
  hout = out[3] #final output from stan
  period = out[4] #study period identifier
  id = out[5] #patient identifier
  t = out[6] #treatment allocations 
  y = out[7] #observed outcomes
  save(out, file=paste("Example2_",designs[i],".RData",sep=""))
}



