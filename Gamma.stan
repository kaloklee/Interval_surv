// Gamma Distribution



data{
  
  int Nuc;
  real yuc[Nuc];
  int Nc;
  real yc[Nc];

}

parameters{
  
  real<lower=0> alpha;
  real<lower=0> r;
  
}
  
model{
  
  
  alpha ~ normal(1,10);
  r ~ normal(1,10);
  
  for (i in 1:Nuc) {
    target += gamma_lpdf(yuc[i] | alpha,r ) ;
  }
  
  for (i in 1:Nc){

      target += gamma_lccdf(yc[i] | alpha,r );
     
              
  }
  
}

generated quantities{
  
  real S[50];
  
  
  for (i in 1:50){
    
    S[i]=exp(gamma_lccdf(i | alpha,r ));
    
  }
  
}

