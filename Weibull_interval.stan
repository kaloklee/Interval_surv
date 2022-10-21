//  Fit a Weibull model with interval censored data
//  


data {

  int<lower=0> T; // calibration period
  
  vector[T] Tstart;
  
  vector[T] Tend; 
  
  vector[T] Dropped; 
  
  int<lower=0> N; //number of observations
  
  real<lower=0> R; // remaining after T

}


parameters {
  
  real<lower=0> lambda;
  real<lower=0> c;
  
}


model {

  //flat priors are used!!!
  
  //HANDLE THE FIRST INTERVAL INVOLVING 0 DIFFERENTLY
  target += Dropped[1] * weibull_lcdf(Tend[1] | c, lambda);
  
  for (i in 2:T) {
                          
   target += Dropped[i] * log_diff_exp(weibull_lccdf(Tend[i] | c, lambda),weibull_lccdf(Tstart[i] | c, lambda));
                           
  }
  target += R * weibull_lccdf(Tend[T] | c, lambda) ;

}


generated quantities{

//create expected historgram based on the model parameters (only consider parameter uncertainty)
//create predicted historgram based on the model parameters and randomness (consider parameter uncertainty and process uncertainty)

   vector[T+1] expected;
   vector[T+1] predicted=rep_vector(0, T+1);
   real chi_sq=0;
   
   for (i in 1:T) {

       expected[i]=N*exp(log_diff_exp( weibull_lccdf(Tstart[i] | c, lambda),
                                        weibull_lccdf(Tend[i] | c, lambda) 
                                      )
                         );
   }
   expected[T+1]=N*exp(weibull_lccdf(Tend[T] | c, lambda));
   
   for (j in 1:N) {
      vector[N] time;
      time[j] = weibull_rng(c, lambda);
      for (t in 1:T)  {
        predicted[t] += ( time[j]> t-1 && time[j]<= t );  
      }
      predicted[T+1] += ( time[j]>T );

  }
  
  
  for (t in 1:T) {
    
    chi_sq += square(Dropped[t]-expected[t])/expected[t];
    
  }
    
    chi_sq += square(R-expected[T+1])/expected[T+1];
  
}





