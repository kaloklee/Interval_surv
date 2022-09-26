library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(posterior)

setwd("~/R/Survival_Models")

#War Data from Schmittlen and Morrison 1980

Dropped <- c(165, 46, 33, 22, 11, 11, 5, 4, 4, 3)
Tstart <- c(0:9);
Tend <- Tstart+1;

data_list <- list(
  Dropped = Dropped, 
  Tstart = Tstart,
  Tend = Tend,
  T = 10,
  N = 315,
  R = 315-sum(Dropped),
  K = 2
)


########################################################################

file <- file.path("Weibull_K_interval.stan")
mod_K<- cmdstan_model(file,stanc_options = list("O1"), quiet=TRUE)

fit_K <- mod_K$sample(
  data = data_list,
  seed = 121,
  chains = 4,
  parallel_chains = 4,
  refresh = 500,
  iter_warmup = 2000  ,
  iter_sampling = 2000  
  ,adapt_delta = 0.99
  ,  max_treedepth = 20
)

fit_K$summary()
fit_K$draws(c( "predicted")) %>% summarise_draws("mean", "quantile2", "rhat", .cores=getOption("mc.cores", 1))
fit_K$draws(c( "expected")) %>% summarise_draws("mean", "quantile2", "rhat", .cores=getOption("mc.cores", 1))

posteriorK <- fit_K$draws(format="df")

######################################################################

file <- file.path("Exp_K_interval.stan")
mod_E<- cmdstan_model(file,stanc_options = list("O1"), quiet=TRUE)

fit_E <- mod_E$sample(
  data = data_list,
  seed = 121,
  chains = 4,
  parallel_chains = 4,
  refresh = 500,
  iter_warmup = 2000  ,
  iter_sampling = 2000  
  ,adapt_delta = 0.99
  ,  max_treedepth = 20
)

fit_E$summary()
fit_E$draws(c( "predicted")) %>% summarise_draws("mean", "quantile2", "rhat", .cores=getOption("mc.cores", 1))
fit_E$draws(c( "expected")) %>% summarise_draws("mean", "quantile2", "rhat", .cores=getOption("mc.cores", 1))


posteriorE <- fit_E$draws(format="df")


######################################################################

file <- file.path("Weibull_interval.stan")
mod_W<- cmdstan_model(file,stanc_options = list("O1"), quiet=TRUE)

fit_W <- mod_W$sample(
  data = data_list,
  seed = 121,
  chains = 4,
  parallel_chains = 4,
  refresh = 500,
  iter_warmup = 1000  ,
  iter_sampling = 1000  
#  ,adapt_delta = 0.95
#  ,max_treedepth = 20
)

#check convergence
fit_W$summary()

#observe wider intervals for "predicted" since it accounts for process uncertainty
fit_W$draws(c( "predicted")) %>% summarise_draws("mean", "quantile2", "rhat", .cores=getOption("mc.cores", 1))
fit_W$draws(c( "expected")) %>% summarise_draws("mean", "quantile2", "rhat", .cores=getOption("mc.cores", 1))

posteriorW <- fit_W$draws(format="df")

#Show the average posterior histogram

E_pred <- fit_E$summary(c("predicted"))  
K_pred <- fit_K$summary(c("predicted"))
W_pred <- fit_W$summary(c("predicted"))

a1 <- data.frame( value = E_pred$mean,
                  type = rep("2seg Exponential",11))
a2 <- data.frame( value = K_pred$mean,
                  type = rep("2seg Weibull",11))
a3 <- data.frame( value = W_pred$mean,
                  type = rep("Weibull",11))
a4 <- data.frame( value = c(Dropped,315-sum(Dropped)),
                  type = c("Actual"))

result<-data.frame(cbind(
  x=rep(c("0-1","1-2","2-3","3-4","4-5","5-6","6-7","7-8","8-9","9-10",">10"),4),
  rbind(a1,a2,a3,a4)))

result$x2 <- factor(result$x, levels=c("0-1","1-2","2-3","3-4","4-5","5-6","6-7","7-8","8-9","9-10",">10"))

ggplot() + 
  geom_bar(data = result, aes(x = x2, y = value, fill = type), position = "dodge", stat = "identity") +
  labs (title = "", x="x", y="y") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(fill=guide_legend(title="")) +
  theme_bw()
#see the graph "histograms.png"

#Compare the chi-sq fit statistics
chi_sq <- posteriorW %>%
  mutate(type = "Weibull") %>%
  select(type, chi_sq) %>%
  bind_rows( posteriorE %>% 
             mutate(type = "2seg Exponential") %>%
             select(type, chi_sq) 
           )%>%
  bind_rows( posteriorK %>% 
             mutate(type = "2seg Weibull") %>%
             select(type, chi_sq)             
            )

ggplot(chi_sq, aes(x=chi_sq, fill=type)) +
  geom_density(alpha=.15) +
  guides(fill=guide_legend(title=""))+
  theme_bw()
#see the graph "chi_sq.png"
