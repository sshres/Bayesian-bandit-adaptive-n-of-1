INDBAN_run <- function(J, nperiod, nobv, b, truep){
  
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  
  hmc_samples <- function(current_data){
    hmodel <- stan_model(file='model_agg.stan') #stan model specified in an external file
    houtput <- sampling(hmodel,
                        data = current_data,
                        chains = 4,             # number of Markov chains
                        warmup = 4000,          # number of warmup iterations per chain
                        iter = 6000,           # total number of iterations per chain
                        refresh = 0,
                        verbose = FALSE,
                        control = list(adapt_delta = 0.99, max_treedepth = 15))
    return(houtput)
  }
  
  hmc_samples_ind <- function(current_data){
    imodel <- stan_model(file='model_ind.stan')
    ioutput <- sampling(imodel,
                        data = current_data,
                        chains = 4,             # number of Markov chains
                        warmup = 4000,          # number of warmup iterations per chain
                        iter = 6000,           # total number of iterations per chain
                        refresh = 0,
                        verbose = FALSE,
                        control = list(adapt_delta = 0.99, max_treedepth = 15))
    return(ioutput)
  }
  N = 2*J*nobv #number of total datapoints
  t = y = id = period = cyc = array(0,N)
  count = 0
  
  #first cycle of treatments 
  for (i in 1:J){
    tr1 = sample(c(0,1),1)  #decide which treatment to choose randomly
    tr2 = 1-tr1   #ensure the other treatment is chosen next
    Trt = c(tr1,tr2)
    for (j in 1:length(Trt)){
      for (z in 1:nobv){
        count = count + 1
        y[count] = b[i,1] + b[i,2]*Trt[j] + rnorm(1,0,tru.y.sig) #generate outcome
        id[count] = i #patient identifier
        period[count] = j
        cyc[count] = 1 #cycle
        t[count] = Trt[j]} #treatment
    }
  }
  
  #generate posterior samples after initial patient data from 2 periods
  current_data = list(J=J,N=N,y=y,t=t,id=id)
  hout = hmc_samples_ind(current_data)
  post_samples = extract(hout)
  
  Trt = c(0,1)
  allocprob = matrix(0,nperiod-2,J) #array to store treatment allocation probabilities
  #for each period
  for (numpd in 3:nperiod){
    for (pid in 1:J){
      nsamp = dim(post_samples$beta1)[1] #check how many MCMC samples
      indA = indB = 0
      for (n in 1:nsamp){
        y_expA = post_samples$beta1[n,pid] + post_samples$beta2[n,pid]*Trt[1]
        y_expB = post_samples$beta1[n,pid] + post_samples$beta2[n,pid]*Trt[2]
        if (y_expA < y_expB){indA = indA + 1} #indicator to check whether A or B is better
        else {indB = indB + 1}                #lower y values implies lower pain score
      }
      
      wA = indA/nsamp #allocation probability 
      wB = indB/nsamp
      
      c=0.5*(numpd/nperiod)   #tuning parameter
      wA = (wA)^c/(wA^c + wB^c) 
      wB = 1- wA
      allocprob[numpd-2,pid] = wA #keep track of allocation prob
      
      t.next = sample(c(0,1), size=1, prob=c(wA,wB))  #treatment chosen based on allocations probs
      
      #generate new observations for given patient based on chosen treatment
      for (z in 1:nobv){
        ynew = b[pid,1] + b[pid,2]*t.next + rnorm(1,0,tru.y.sig)
        y = c(y,ynew)
        id = c(id,pid)
        period = c(period, numpd)
        t = c(t, t.next)}
      
}#posterior updating
    nn = length(y)
    current_data = list(J=J,N=nn,y=y,t=t,id=id)
    hout = hmc_samples_ind(current_data)
    post_samples = extract(hout)
  }

  hout1 = hmc_samples(current_data)
  post_samples_aggr = extract(hout1)
  
  out <- list(allocprob,post_samples,hout,period,id,t,y)
  return(out)
}

