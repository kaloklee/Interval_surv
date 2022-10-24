library(tidyverse)


#War Data from Schmittlen and Morrison 1980

Dropped <- c(165, 46, 33, 22, 11, 11, 5, 4, 4, 3)
Tstart <- c(0:9);
Tend <- Tstart+1;
R = 315-sum(Dropped);

#Create a logsumexp function to avoid computation error
logsumexp <- function (a,b) {
  
  c = pmax(a,b);
  return (c + log(exp(a-c)+exp(b-c))) ;
  
}

logdiffexp <- function (a,b) {
  
  c = pmax(a,b);
  return (c + log(exp(a-c)-exp(b-c))) ;
  
}


Weibull_Exp_LL <- function(par,Tstart,Tend,Dropped,R) {
  
  r=par[1];
  shape=par[2];
  scale=par[3];
  p=par[4];
  t=length(Tend);
  
  #uncensored part
  #Exponential piece
  LL1=log(1-p)+logdiffexp( 
    pexp(Tstart, r, lower = FALSE, log = TRUE),
    pexp(Tend, r, lower = FALSE, log = TRUE));

  #Weibull piece
  LL2=log(p)+logdiffexp( 
    pweibull(Tstart, shape, scale, lower = FALSE, log = TRUE),
    pweibull(Tend, shape, scale, lower = FALSE, log = TRUE));
  
  #combined
  LLu = Dropped*logsumexp(LL1,LL2);
  
  #censored part
  
  #Exponential piece
  LL1=log(1-p)+pexp(Tend[t], r, lower = FALSE, log = TRUE);
    
  #Weibull piece
  LL2=log(p)+pweibull(Tend[t], shape, scale, lower = FALSE, log = TRUE);
  
  #combined
  LLc = R*logsumexp(LL1,LL2);
  
  return (-sum(LLu)-sum(LLc));
}


solution=optim(par=c(2.5,1,2.5,.5),fn=Weibull_Exp_LL,
               Tstart=Tstart,Tend=Tend,Dropped=Dropped,R=R,
               method="L-BFGS-B",
               lower=c(1e-5,1e-5,1e-5,1e-5), upper=c(Inf,Inf,Inf,.999))


###Survival function
WE_exp <- function(par,Tstart,Tend,N) {
  
  r=par[1];
  shape=par[2];
  scale=par[3];
  p=par[4];
  t=length(Tend);
  
  lp1 = log(1-p) + logdiffexp(pexp(Tstart, r, lower = FALSE, log = TRUE),
                                 pexp(Tend, r, lower = FALSE, log = TRUE));
    
  lp2 = log(p) + logdiffexp(pweibull(Tstart, shape, scale, lower = FALSE, log = TRUE),
                            pweibull(Tend, shape, scale, lower = FALSE, log = TRUE));
  
  Exp1 = N*exp(logsumexp(lp1,lp2));
  
  lp3 = log(1-p) + pexp(Tend[t], r, lower = FALSE, log = TRUE);
  lp4 = log(p) + pweibull(Tend[t], shape, scale, lower = FALSE, log = TRUE);

  Exp2 = N*exp(logsumexp(lp3,lp4));
    
  return (c(Exp1, Exp2));
  
}
  
Expected = WE_exp(solution$par,Tstart,Tend,315)

Expected

a1 <- data.frame( value = Expected,
                  type = c("Model"))
a2 <- data.frame( value = c(Dropped,315-sum(Dropped)),
                  type = c("Actual"))

result<-data.frame(cbind(
  x=c("0-1","1-2","2-3","3-4","4-5","5-6","6-7","7-8","8-9","9-10",">10"),
  rbind(a1,a2)))

result$x2 <- factor(result$x, levels=c("0-1","1-2","2-3","3-4","4-5","5-6","6-7","7-8","8-9","9-10",">10"))

ggplot() + 
  geom_bar(data = result, aes(x = x2, y = value, fill = type), position = "dodge", stat = "identity") +
  labs (title = "", x="x", y="y") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(fill=guide_legend(title="")) +
  theme_bw()
