//  Fit a K-segment Weibull model with interval censored data
//  


data {

  int<lower=0> T; // calibration period
  
  vector[T] Tstart;
  
  vector[T] Tend; 
  
  vector[T] Dropped; 
  
  int<lower=0> N; //number of observations
  
  int<lower=0> R; // remaining after T
  
  int<lower=0> K; //number of segments

}


parameters {
  
  simplex[K] w;  //segment weights
  vector<lower=0>[K] c;
  positive_ordered[K] lambda;  //using an ordered vector for identification
  
}


model {

  vector[K] log_w = log(w); 
  vector[K] lps1;
  vector[K] lps3;

  lambda ~ gamma(1,3);
  c ~ gamma(1,3);

  //HANDLE THE FIRST INTERVAL INVOLVING 0 DIFFERENTLY
  for (k in 1:K) {
      lps1[k] = log_w[k] + weibull_lcdf(Tend[1] | c[k], lambda[k]) ;
  }
  target += Dropped[1] * log_sum_exp(lps1);
  
  for (i in 2:T) {
   vector[K] lps2;
   for (k in 1:K) {
      lps2[k] = log_w[k] + 
                log_diff_exp(weibull_lcdf(Tend[i] | c[k], lambda[k]),weibull_lcdf(Tstart[i] | c[k], lambda[k]));
   }  
   target += Dropped[i] * log_sum_exp(lps2);
  }
  
  for (k in 1:K) {
    lps3[k] = log_w[k] + weibull_lccdf(T | c[k], lambda[k]);
  }
  target += R * log_sum_exp(lps3) ;

}


generated quantities{

//create expected historgram based on the model parameters (only consider parameter uncertainty)
//create predicted historgram based on the model parameters and randomness (consider parameter uncertainty and process uncertainty)

    vector[T+1] expected;
    vector[T+1] predicted=rep_vector(0, T+1);
    real chi_sq=0;

    for (i in 1:T+1) {
      vector[K] lps2;
      for (k in 1:K) {
        if (i <= T) 
        
            lps2[k] = log( w[k] ) + log_diff_exp(weibull_lccdf(i-1 | c[k],lambda[k]),
                                               weibull_lccdf(i | c[k],lambda[k])
                                              );
        
        else lps2[k] = log( w[k] ) + weibull_lccdf( T+1 | c[k],lambda[k] ) ;
      }
    
      expected[i]=N*exp(log_sum_exp(lps2));
    }
    
    
    

    for (j in 1:N) {
       vector[N] time;
       int seg;
       seg = categorical_rng(w);
       time[j] = weibull_rng(c[seg],lambda[seg]);
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

