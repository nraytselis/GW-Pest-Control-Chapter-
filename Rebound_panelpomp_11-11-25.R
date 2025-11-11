library(pomp)
library(panelPomp)
library(tidyverse)
library(adaptMCMC)

samps = readRDS("Rebound_parameters_death.RDA")
params = samps$samples[which.max(samps$log.p),]
# names(params)[5] = "comp_d"
# params = c(params, "comp_b" = 0.001)
# Imports the data and experimental conditions from the experiment
Multiple_Time_Series <- read.csv("C:/RData/Data_Rebound_4-2-25.csv")
exp_df <- read.csv("C:/RData/Rebound_experimental_conditions.csv")

# Defines the days that harvesting actually occurred
Harvest_days = c(8, 36, 64)

# This converts target of proportion mortality imposed for a single day into an instantaneous death rate for use in model
exp_df[,"Harvest_rate"] = -log(1 - exp_df[,"Harvest"]/100)*(exp_df$Harvest_days %in% Harvest_days)

# Initial conditions function for pomp
GW_rinit <- Csnippet("
  N = rpois(N0);
  J = rpois(J0);
  A = rpois(A0);
")

# r process model
GW_step_process <- Csnippet("
  double VOL = 15;
  
  double Pred_A = Harvest_rate;
  double Pred_J = Harvest_rate;
  double Pred_N = Harvest_rate * Sieve;
  
  double birth_rate = b_M*exp(-comp_b / VOL * (c_N * N + c_J * J + A));

  double d_A_c = d_A*exp(comp_d / VOL * (c_N * N + c_J * J + A));
  double d_J_c = d_J*exp(comp_d / VOL * (c_N * N + c_J * J + A));
  double d_N_c = d_N*exp(comp_d / VOL * (c_N * N + c_J * J + A));



  /* Specify adult rates */
  const double A_rates[2] = {d_A_c, Pred_A};
  double A_trans[2];
  reulermultinom(2, A, &A_rates, dt, &A_trans);
  
  /* Specify juvenile rates */
  const double J_rates[3] = {m_J, d_J_c, Pred_J};
  double J_trans[3];
  reulermultinom(3, J, &J_rates, dt, &J_trans);
  
  /* Specify nauplii rates */
  const double N_rates[3] = {m_N, d_N_c, Pred_N};
  double N_trans[3];
  reulermultinom(3, N, &N_rates, dt, &N_trans);
  
  /* births */
  int births = rpois(birth_rate * A / 2 * dt);
  
  /* New totals */
  int A_out = A - A_trans[0] - A_trans[1] + J_trans[0];
  int J_out = J - J_trans[0] - J_trans[1] - J_trans[2] + N_trans[0];
  int N_out = N - N_trans[0] - N_trans[1] - N_trans[2] + births;
  N = N_out >= 0 ? N_out : 0;
  J = J_out >= 0 ? J_out : 0;
  A = A_out >= 0 ? A_out : 0;
")

# r measure model
GW_rmeasure <- function(N, J, A, counted_volume, VOL = 15, SAMPLEVOL = 0.5, aliquot_volume = 10, k, ...) {
  c(
    N1 = rnbinom(n = 1, mu = N * counted_volume / aliquot_volume * SAMPLEVOL / VOL, size = k),
    JOA1 = rnbinom(n = 1, mu = (J + 0.5 * (1 + 1 / 3) * A) * counted_volume / aliquot_volume * SAMPLEVOL / VOL, size = k),
    AF1 = rnbinom(n = 1, mu = 0.5 * 2 / 3 * A * counted_volume / aliquot_volume * SAMPLEVOL / VOL, size = k)
  )
}
# d measure model
GW_dmeasure <- Csnippet("
  double VOL = 15;
  double SAMPLEVOL = 0.5;
  double aliquot_volume = 10;
  lik = dnbinom_mu(N1, k, fmax(N, 1) * counted_volume / aliquot_volume * SAMPLEVOL / VOL, give_log) +
        dnbinom_mu(JOA1, k, fmax(J + 0.5 * (1 + 1.0 / 3.0) * A, 1) * counted_volume / aliquot_volume * SAMPLEVOL / VOL, give_log) +
        dnbinom_mu(AF1, k, fmax(0.5 * 2.0 / 3.0 * A, 1) * counted_volume / aliquot_volume * SAMPLEVOL / VOL, give_log);
")

# Template with updated parameters
template <- pomp(
  data = data.frame(
    Day = 0:17*7, #This needs to match the times in the actual dataset  (Generalize this)
    Tank = NA,
    N1 = NA, JOA1 = NA, AF1 = NA
  ),
  times = "Day", t0 = 0,
  rprocess = discrete_time(GW_step_process, delta.t = 1),
  rmeasure = GW_rmeasure,
  dmeasure = GW_dmeasure,
  rinit = GW_rinit,
  obsnames = c("N1", "JOA1", "AF1"),
  statenames = c("N", "J", "A"),
  covarnames = c("Sieve", "Harvest_rate", "counted_volume"),
  paramnames = names(params),
  params = params
)

# Checking how many tanks there are
U <- length(unique(exp_df$Tank))
# Empty list for the pomp templates
poList <- setNames(vector(mode = "list", length = U),
                   nm = paste0("unit", 1:U))
# Populates each element of the list with the model template, covariates, and the data
for (u in seq_len(U)) {
  covariates <- subset(exp_df, Tank == u, select = c("Harvest_days", "Sieve", "Harvest_rate", "counted_volume"))
  cov_table <- covariate_table(covariates, times = "Harvest_days")
  poList[[u]] <- pomp(template, covar = cov_table)
  data_u <- pomp(
    data = subset(Multiple_Time_Series, Tank == u, select = c("Day", "Tank", "N1", "JOA1", "AF1")),
    times="Day",
    t0=timezero(template)
  )
  poList[[u]]@data <- data_u@data
}

# Creates the panel Pomp object
GW_panel <- panelPomp(object = poList, shared = coef(template))

# Simulate and plot the model
#plot(simulate(GW_panel, nsim = 1))

prior.likelihood = function(x){
  prior.lik = with(as.list(x),
                   sum(dbeta(c(c_J, c_N), 1, 1, log=T)) + # These parameters are unknown between 0-1
                     dunif(k, min=0, max=100, log=T) + # aggregation parameter for observatiosn, k < 100
                     dnorm(b_M, mean=5.67, sd=0.19, log=T) + # Estimated from the lab (egg data)
                     sum(dunif(c(N0, J0, A0), min=0, max=1000000, log=T)) + # maximum search rate of type 2 fxn response and N_0s are unknown
                     #dunif(m_N, min=0.025, max=1, log=T) + # Maximum birth rate is less than 100/day per female, k is < 100
                     #dunif(m_J, min=0.025, max=0.5, log=T) + # Maximum birth rate is less than 100/day per female, k is < 100
                     dnorm(m_N, mean = 0.457, sd = 0.11, log=T) + # maturation rate estimated in lab
                     dnorm(m_J, mean=0.067, sd = 0.035, log=T) + # maturation rate estimated in lab
                     sum(dunif(c(d_N, d_J, d_A), min=0, max=1, log=T)) + # maturation rates are less than 0.33
                     sum(dunif(c(comp_b, comp_d), min=0, max=0.01, log=T)) # competition strength is less than 0.01
  )
  return(prior.lik)
}

prior.likelihood(params)

full_LL = function(x){
  prior_LL = prior.likelihood(x)
  if(!is.finite(prior_LL)){return(prior_LL)}else{
    prior_LL  + pfilter(GW_panel, shared = x, Np=500)@ploglik
  }
}

full_LL(params)

# # Read in the previous best-fit parameters and estimated var-covar matrix
# samps = readRDS("Rebound_parameters_death.RDA")
# pars = samps$samples[which.max(samps$log.p),]
# variances = samps$cov.jump
# 
# full_LL(pars)
# 
# model_fit = MCMC(full_LL, init=pars, scale=as.matrix(variances), adapt=50000, acc.rate = 0.3, n=250000)
# 
# # Plot trajectory of LL
# par(mfrow=c(1,1))
# plot(model_fit$log.p)
# 
# # Save result of adaptive MCMC
# saveRDS(model_fit, file = "Rebound_parameters_death.RDA")
# 
# # Return best fit parameters and best Full LL
# round(model_fit$samples[which.max(model_fit$log.p),], 5)
# max(model_fit$log.p)

# Read in the previous best-fit parameters and estimated var-covar matrix
samps = readRDS("Rebound_parameters_death.RDA")
pars = samps$samples[which.max(samps$log.p),]
variances = samps$cov.jump

full_LL(pars)
# variances = diag(length(params))*1e-2*params
# pars = params
model_fit = MCMC(full_LL, init=pars, scale=as.matrix(variances), adapt=50000, acc.rate = 0.3, n=250000)

# Plot trajectory of LL
par(mfrow=c(1,1))
plot(model_fit$log.p)

# Save result of adaptive MCMC
saveRDS(model_fit, file = "Rebound_parameters3.RDA")

# Return best fit parameters and best Full LL
signif(model_fit$samples[which.max(model_fit$log.p),], 3)
max(model_fit$log.p)


# Read in the previous best-fit parameters and estimated var-covar matrix
samps = readRDS("R:/CivitelloLab/Rebound_parameters3.RDA")
pars = samps$samples[which.max(samps$log.p),]
variances = samps$cov.jump

# full_LL(pars)
# 
# model_fit = MCMC(full_LL, init=pars, scale=as.matrix(variances), adapt=5000, acc.rate = 0.3, n=25000)
# 
# # Plot trajectory of LL
# par(mfrow=c(1,1))
# plot(model_fit$log.p)
# 
# # Save result of adaptive MCMC
# saveRDS(model_fit, file = "Rebound_parameters_death3.RDA")
# 
# # Return best fit parameters and best Full LL
# round(model_fit$samples[which.max(model_fit$log.p),], 5)
# max(model_fit$log.p)


# Do not use the simulate function to generate predictions, use the pfilter function (following Martinez-Bakker et al)
#sim.data = simulate(GW_panel, nsim=1, shared=pars)

# Using pfilter lets you average the prediction (pred.mean) over all of the particles. Martinez et al. used 2000
sim.data = pfilter(GW_panel, shared=pars, Np=2000, pred.mean=T)

# Function needed to extract pred.mean data from a panelpomp unit object
extract_predictions = function(x){
  x@pred.mean
}

# Lapply to extract simulated data from the panelpomp unit object list
data.list = lapply(FUN=extract_predictions, X=sim.data@unit.objects)

# Converts this list of simulated data to a dataframe
sim.df = data.frame("Day" = (0:17)*7, t(data.list[[1]]), "Tank" = 1)
for(i in 2:60){
  sim.df = rbind(sim.df, data.frame("Day" = (0:17)*7, t(data.list[[i]]), "Tank" = i))
}

write.csv(sim.df, file = "Harvest_sim_output.csv")

head(Multiple_Time_Series)

plotting_predictions = left_join(Multiple_Time_Series, sim.df, by=c("Day", "Tank"))
#colnames(plotting_predictions)[c(7:9)] = c("AF1", "JOA1", "N1")

exp_df_summary = exp_df %>% group_by(Tank) %>% summarize(Tank = mean(Tank), Harvest=max(Harvest), Sieve=mean(Sieve))

plotting_predictions = left_join(plotting_predictions, exp_df_summary, by="Tank")
plotting_predictions = plotting_predictions %>% group_by(Harvest, Sieve, Day) %>% summarise(AF1_p = (1/6)*mean(A), JOA1_p = (5/6)*mean(A) + mean(J), N1_p = mean(N),N1 = mean(Nauplii_per_Tank), AF1 = mean(Adults_per_Tank), JOA1 = mean(Juveniles_per_tank))

ggplot(data=plotting_predictions, aes(x=Day, group=interaction(Harvest, Sieve), color=interaction(Harvest, Sieve))) + geom_line(aes(y=AF1_p)) +
  theme(legend.position = "none")

ggplot(data=plotting_predictions, aes(x=Day, group=interaction(Harvest, Sieve), color=interaction(Harvest, Sieve))) + geom_line(aes(y=JOA1_p)) +
  theme(legend.position = "none")

ggplot(data=plotting_predictions, aes(x=Day, group=interaction(Harvest, Sieve), color=interaction(Harvest, Sieve))) + geom_line(aes(y=N1_p)) +
  theme(legend.position = "none")

ggplot(data=plotting_predictions, aes(x=N1+1, y=N1_p+1)) + geom_point() + geom_abline(slope=1, intercept=0) + scale_x_log10() + scale_y_log10()
ggplot(data=plotting_predictions, aes(x=JOA1+1, y=JOA1_p+1)) + geom_point() + geom_abline(slope=1, intercept=0) + scale_x_log10() + scale_y_log10()
ggplot(data=plotting_predictions, aes(x=AF1+1, y=AF1_p+1)) + geom_point() + geom_abline(slope=1, intercept=0) + scale_x_log10() + scale_y_log10()


##### Work in progress - R2 for panelpomp #####

R2_prediction_for_Rebound = function(parameters, panel, dataset= read.csv("C:/RData/Data_Rebound_4-2-25.csv")){
  # Using pfilter lets you average the prediction (pred.mean) over all of the particles. Martinez et al. used 2000
  sim.data = pfilter(panel, shared=parameters, Np=2000, pred.mean=T)
  
  # Lapply to extract simulated data from the panelpomp unit object list
  data.list = lapply(FUN=extract_predictions, X=sim.data@unit.objects)
  
  # Converts this list of simulated data to a dataframe
  sim.df = data.frame("Day" = (0:17)*7, t(data.list[[1]]), "Tank" = 1)
  for(i in 2:60){
    sim.df = rbind(sim.df, data.frame("Day" = (0:17)*7, t(data.list[[i]]), "Tank" = i))
  }
  
  obs_and_preds = left_join(dataset, sim.df, by=c("Day", "Tank"))
  obs_and_preds = obs_and_preds %>% mutate(AF_p = (1/6)*A, JOA_p = (5/6)*A + J, N_p = N,N_o = Nauplii_per_Tank, 
                                                       AF_o = Adults_per_Tank, JOA_o = Juveniles_per_tank)
  
  # Need to worry about complete cases because there are three missing observations (this is ok w/ na.rm=T)

  # These are means over the whole experiment, not tank-specific means, which maybe they should be (check paper)
  mean_AF = mean(obs_and_preds$AF_o, na.rm=T)
  mean_JOA = mean(obs_and_preds$JOA_o, na.rm=T)
  mean_N = mean(obs_and_preds$N_o, na.rm=T)

  # This is not log-transformed, though i think it should be (check paper)
  R2_A = 1 - sum((obs_and_preds$AF_p - obs_and_preds$AF_o)^2, na.rm=T)/sum((mean_AF - obs_and_preds$AF_o)^2, na.rm=T)
  R2_A
}

R2_prediction_for_Rebound(parameters=pars, panel=GW_panel)

print(o$AF_o, n.print=1e6)
#
#
#
# #### To do: Get tank-level abundance estimates to line up with these predictions ####
