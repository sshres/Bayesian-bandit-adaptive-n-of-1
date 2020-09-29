
source("INDBAN_run.R")

#indiv parameters
J = 27 #number of patients
nperiod = 6 #number of periods (two per cycle)
nobv = 30 #number of observations per patient
#parameters from Stunnenburg et al 2018 supplementary web material
b0 = c(7.87,6.93,4.93,7.07,5.14,5.09,4.71,0.00,3.56,4.23,5.5,7.07,6.62,
       5.87,5.71,3.31,6.08,2.77,8,7.71,6.8,6.57,6.8,6.22,6.71,6.71,4.62)
b1 = c(-4.4,-3.1,-1.1,-4.7,-1.2,-3.2,-3.7,0,-2.66,0.05,-3.5,-2.7,-4.49,-3.22,
       -2.55,-2.63,-1.13,-4.41,-0.01,-5.49,-2.64,-1.65,-4.65,-3.14,-4.50,-6.99,-3.57)
b = cbind(b0,b1)
tru.y.sig = 2 #within-patient deviance

trup = c(tru.y.sig,b0,b1)

out <- INDBAN_run(J, nperiod, nobv, b, trup) #Individual Bayesian bandit adaptive design

allocprob = out[1] #allocation probabilities through the adaptive periods
post_samples = out[2] #final posterior draws
hout = out[3] #final output from stan
period = out[4] #study period identifier
id = out[5] #patient identifier
t = out[6] #treatment allocations 
y = out[7] #observed outcomes
save(out, file=paste("RealEg.RData",sep=""))





