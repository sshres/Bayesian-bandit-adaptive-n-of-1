FRN_run <- function(J, nperiod, nobv, b, truep){
  
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)

  numt = 2
  N = nperiod*J*nobv #number of total datapoints
  t = y = id = period = cyc = array(0,N)
  count = 0
  
  for (k in 1:(nperiod/2)){
    for (i in 1:J){
      tr1 = sample(c(0,1),1)  #this part decides which treatment to choose randomly
      tr2 = 1-tr1   #this part ensures the other treatment is chosen next
      Trtreg = c(tr1,tr2)
      
      for (j in 1:numt){
        for (z in 1:nobv){
          count = count + 1
          y[count] = b[i,1] + b[i,2]*Trtreg[j] + rnorm(1,0,tru.y.sig) #observe outcome
          id[count] = i
          period[count] = j+2*(k-1)
          cyc[count] = k
          t[count] = Trtreg[j] #save treatment choice
        }
      }
    }
  }
  
  current_data = list(J=J,N=N,y=y,t=t,id=id)
  
  hout <- stan(file='model_agg.stan', 
                  data = current_data,
                  chains = 4,             # number of Markov chains
                  warmup = 1000,          # number of warmup iterations per chain
                  iter = 3000,           # total number of iterations per chain
                  refresh = 0,
                  verbose = FALSE,
                  control = list(adapt_delta = 0.99,max_treedepth = 15))
  
  post_samples=extract(hout)
  out <- list(post_samples,hout,period,id,t,y)
  return(out)
}





