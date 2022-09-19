//  Fit a Weibull model with interval censored data
//  


data {

  int<lower=0> T; // calibration period
  
  vector[T] Tstart;
  
  vector[T] Tend; 
  
  vector[T] Dropped; 
  
  int<lower=0> N; 
  
  int<lower=0> R; // remaining buyers

}


parameters {
  
  real<lower=0> lambda;
  real<lower=0> c;
  
}


model {

  lambda ~ normal(10,10);
  c ~ normal(1,10);
  
  //HANDLE THE FIRST INTERVAL INVOLVING 0 DIFFERENTLY
  target += Dropped[1] * weibull_lcdf(Tend[1] | c, lambda);
  
  for (i in 2:T) {
                          
   target += Dropped[i] * log_diff_exp(weibull_lcdf(Tend[i] | c, lambda),weibull_lcdf(Tstart[i] | c, lambda));
                           
  }
  target += R * weibull_lccdf(T | c, lambda) ;

}



