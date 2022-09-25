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
  R = 315-sum(Dropped)
)

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
  #  ,adapt_delta = 0.95,
  #  max_treedepth = 20
)

#check convergence
fit_W$draws(c( "lambda","c","chi_sq")) %>% summarise_draws("mean", "quantile2", "rhat", .cores=getOption("mc.cores", 1))

#observe wider intervals for "predicted" since it accounts for process uncertainty
fit_W$draws(c( "predicted")) %>% summarise_draws("mean", "quantile2", "rhat", .cores=getOption("mc.cores", 1))
fit_W$draws(c( "expected")) %>% summarise_draws("mean", "quantile2", "rhat", .cores=getOption("mc.cores", 1))


posteriorW <- fit_W$draws(format="df")

mcmc_trace(posteriorW, pars = "c")
mcmc_trace(posteriorW, pars = "lambda")

mcmc_intervals(posteriorW, pars=c("lambda","c"))

mcmc_areas(posteriorW, prob = 0.90, prob_outer = 1, pars=c("lambda","c") )

color_scheme_set("red")
mcmc_pairs(posteriorW, pars=c("lambda","c"))

color_scheme_set("blue")
mcmc_scatter(posteriorW,pars=c("lambda","c")) 
mcmc_hex(posteriorW,pars=c("lambda","c"))


file <- file.path("WeibullGamma_interval.stan")
mod_WG<- cmdstan_model(file,stanc_options = list("O1"), quiet=TRUE)

fit_WG <- mod_WG$sample(
  data = data_list,
  seed = 122,
  chains = 4,
  parallel_chains = 4,
  refresh = 500,
  iter_warmup = 1000  ,
  iter_sampling = 1000  
  #  ,adapt_delta = 0.95,
  #  max_treedepth = 20
)

#Check convergence
fit_WG$draws(c( "alpha","c", "r")) %>% summarise_draws("mean", "quantile2", "rhat", .cores=getOption("mc.cores", 1))
fit_WG$draws(c( "r_over_a")) %>% summarise_draws("mean", "quantile2", "rhat", .cores=getOption("mc.cores", 1))
fit_WG$draws(c( "inv_a")) %>% summarise_draws("mean", "quantile2", "rhat", .cores=getOption("mc.cores", 1))
fit_WG$draws(c( "chi_sq")) %>% summarise_draws("mean", "quantile2", "rhat", .cores=getOption("mc.cores", 1))

#observe wider intervals for "predicted" since it accounts for process uncertainty
fit_WG$draws(c( "predicted")) %>% summarise_draws("mean", "quantile2", "rhat", .cores=getOption("mc.cores", 1))
fit_WG$draws(c( "expected")) %>% summarise_draws("mean", "quantile2", "rhat", .cores=getOption("mc.cores", 1))

posterior1 <- fit_WG$draws(format="df")
mcmc_trace(posterior1, pars = "alpha")
mcmc_trace(posterior1, pars = "c")
mcmc_trace(posterior1, pars = "r")
mcmc_trace(posterior1, pars = "r_over_a")
mcmc_trace(posterior1, pars = "inv_a")

mcmc_intervals(posterior1, pars="r")
mcmc_intervals(posterior1, pars="c")
mcmc_intervals(posterior1, pars="alpha")

mcmc_areas(posterior1, prob = 0.90, prob_outer = 1, pars = 'r' )
mcmc_areas(posterior1, prob = 0.90, prob_outer = 1, pars = 'c' )
mcmc_areas(posterior1, prob = 0.90, prob_outer = 1, pars = 'alpha' )

mcmc_dens_overlay(posterior1, pars=c("r", "alpha","c")) 

color_scheme_set("red")
mcmc_pairs(posterior1, pars=c("r", "alpha"))
mcmc_pairs(posterior1, pars=c("inv_a", "r_over_a"))

color_scheme_set("blue")
mcmc_scatter(
  posterior1,
  pars = c("inv_a", "r_over_a"),
  np = nuts_params(fit_WG), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)
mcmc_hex(posterior1, pars = c("inv_a", "r_over_a"))

 
#Show the average posterior histogram

WG_pred <- fit_WG$summary(c("predicted"))  
W_pred <- fit_W$summary(c("predicted"))


a1 <- data.frame( value = WG_pred$mean,
                  type = rep("WG",11))
a2 <- data.frame( value = W_pred$mean,
                  type = rep("Weibull",11))
a3 <- data.frame( value = c(Dropped,315-sum(Dropped)),
                  type = c("Actual"))

result<-data.frame(cbind(
  x=rep(c("0-1","1-2","2-3","3-4","4-5","5-6","6-7","7-8","8-9","9-10",">10"),3),
  rbind(a1,a2,a3)))

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
  bind_rows( posterior1 %>% 
               mutate(type = "WG") %>%
               select(type, chi_sq)
             )

ggplot(chi_sq, aes(x=chi_sq, fill=type)) +
  geom_density(alpha=.2) +
  guides(fill=guide_legend(title=""))+
  theme_bw()
#see the graph "chi_sq.png"