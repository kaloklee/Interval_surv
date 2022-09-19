library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(posterior)

#User data path
setwd("~/R/Generalized Gamma")  

#Data from Schmittlen and Morrison 1980

Dropped <- c(165, 46, 33, 22, 11, 11, 5, 4, 4, 3)
Tstart <- c(0:9);
Tend <- Tstart+1;

data_list <- list(
  Dropped = Dropped, 
  Tstart = Tstart,
  Tend = Tend,
  T = 10,
  N = 315,
  R = 315-sum(Dropped)
)

#Pick a file
file <- file.path("Weibull_interval.stan")
modv7<- cmdstan_model(file,stanc_options = list("O1"), quiet=TRUE)

fit_v7 <- modv7$sample(
  data = data_list,
  seed = 120,
  chains = 3,
  parallel_chains = 3,
  refresh = 500,
  iter_warmup = 1000  ,
  iter_sampling = 1000  
  #  ,adapt_delta = 0.95,
  #  max_treedepth = 20
)

fit_v7$draws(c( "lambda","c")) %>% summarise_draws("mean", "quantile2", "rhat", .cores=getOption("mc.cores", 1))





#calculate average conditional expectation

WG_pred <- fit_wgt$summary(c("S"))  

a1 <- data.frame( value = WG_pred$mean,
                  type = rep("Predicted",11))

a2 <- data.frame( value = c(Dropped,315-sum(Dropped)),
                  type = c("Actual"))

result<-data.frame(cbind(
  x=rep(c("0-1","1-2","2-3","3-4","4-5","5-6","6-7","7-8","8-9","9-10",">10"),2),
  rbind(a1,a2)))

result$x2 <- factor(result$x, levels=c("0-1","1-2","2-3","3-4","4-5","5-6","6-7","7-8","8-9","9-10",">10"))

ggplot() + 
  geom_bar(data = result, aes(x = x2, y = value, fill = type), position = "dodge", stat = "identity") +
  labs (title = "Fit of Gamma", x="x", y="y") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(fill=guide_legend(title=""))
