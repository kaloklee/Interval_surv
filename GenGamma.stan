// Generalized Gamma as implemented in flexsurv R package
// https://search.r-project.org/CRAN/refmans/flexsurv/html/GenGamma.html
// only require 'sigma' >0
  
functions{

real generalized_gamma_lpdf_fs(real x, real mu, real sigma, real Q) {

  real w = (log(x)-mu)/sigma;
  real Qis = 1/(Q * Q);
  real Qw = Q * w;
  real ll;
  if (Q < 0) ll=-log(sigma*x) + log(-Q)*(1 - 2*Qis) 
              + Qis*(Qw - exp(Qw)) - lgamma(Qis);
  if (Q > 0) ll=-log(sigma*x) + log(Q)*(1 - 2*Qis) 
              + Qis*(Qw - exp(Qw)) - lgamma(Qis);
  else ll=std_normal_lpdf(w); 
  return ll;
}

real generalized_gamma_cdf_fs(real x, real mu, real sigma, real Q) {
  real w = (log(x)-mu)/sigma;
  real Qis = 1/(Q * Q);
  real v = Qis*exp(Q * w);
  real lik;
  if (Q > 0) lik=gamma_cdf(v,Qis,1);
  else if (Q < 0) lik=1-gamma_cdf(v,Qis,1);
  else if (Q==0) lik=std_normal_cdf(w);  
  return lik;

 }

}

data{
  
  int Nuc;
  real yuc[Nuc];
  int Nc;
  real yc[Nc];
  
}

parameters{
  
  real mu;
  real<lower=0> sigma;
  real Q;
  
}
  
model{
  
  mu ~ normal(0,5);
  sigma ~ normal(0,1);
  Q ~ normal(0,5);
  
  for (i in 1:Nuc) {
    target += generalized_gamma_lpdf_fs(yuc[i],mu,sigma,Q);
  }
  
  for (i in 1:Nc){
    target += log1p(generalized_gamma_cdf_fs(yc[i],mu,sigma,Q));
  }
  
}


// generated quantities{
//   
//   real S[60];
// 
//     for (i in 1:60){
//       S[i] = 1-generalized_gamma_cdf_2S(i,k,mu,sigma);
//     }
//   
//   
//   
// }

