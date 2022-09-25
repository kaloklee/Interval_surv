//  Fit a Weibull-Gamma model with interval censored data
//  


functions {

  real gamma_Weibull_log_surv(real t , real alpha, real r, real c) {

    return   -r *  log1p( pow(t,c) / (alpha) ) ;
    
  }

}



data {

  int<lower=0> T; // calibration period
  
  vector[T] Tstart;
  
  vector[T] Tend; 
  
  vector[T] Dropped; 
  
  int<lower=0> N; 
  
  int<lower=0> R; // remaining people

}


parameters {
  
  real<lower=0> inv_a; // 1/alpha;
  real<lower=0> r_over_a; // r/alpha, mean of Gamma dist
  real<lower=0> c;

}

transformed parameters{
  
  real<lower=0> alpha = inv(inv_a);
  real<lower=0> r = alpha*r_over_a;
}

model {

  //weakly informative priors for regularization

  r_over_a ~ normal(1,10);
  inv_a ~ normal(1,10);
  
  for (i in 1:T) {

  target += Dropped[i]*log_diff_exp( gamma_Weibull_log_surv(Tstart[i],alpha, r, c) ,
                                     gamma_Weibull_log_surv(Tend[i],alpha, r, c) 
                                    ) ;              
  }
  
   target +=  R * (gamma_Weibull_log_surv(T,alpha, r, c));

}


generated quantities{

   vector[T+1] expected;
   vector[T+1] predicted=rep_vector(0, T+1);
   real chi_sq=0;
   
   for (i in 1:T) {

       expected[i]=N*exp(log_diff_exp( gamma_Weibull_log_surv(i-1,alpha, r, c),
                                       gamma_Weibull_log_surv(i  ,alpha, r, c) 
                                      )
                         );
   }
   expected[T+1]=N*exp(gamma_Weibull_log_surv(T+1,alpha, r, c));
   
   for (j in 1:N) {
      vector[N] time;
      vector[N] F;
      F[j] = uniform_rng(0,1);
      time[j] = pow(-alpha*(1 - pow(1 - F[j],-inv(r))),inv(c));
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


