library(panelPomp)
library(adaptMCMC)
library(pomp)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(stringr)
library(NatParksPalettes)

####Summary####
#This code plots the model fits from the mcmc chains against the experimental data from the Pest Control experiment. 

setwd("~/Desktop/Rscripts/Data")
###processing full chains

fullA = readRDS("Rebound_parameters_full_disp250kff_5.RDA")
fullB = readRDS("Rebound_parameters_full_disp250kff_3.RDA")
fullC = readRDS("Rebound_parameters_full_disp250kff_05.RDA")
fullD = readRDS("Rebound_parameters_full_disp250kff_01.RDA")
full = c(fullA, fullB, fullC, fullD)


get_best_fit = function(chain.list){
  L = length(chain.list)
  chain.scores = numeric()
  for(i in 1:L){
    chain.scores[i] = max(chain.list[[i]]$log.p)
  }
  list(chain.list[[which.max(chain.scores)]]$samples[which.max(chain.list[[which.max(chain.scores)]]$log.p),],
       chain.list[[which.max(chain.scores)]]$cov.jump,
       max(chain.list[[which.max(chain.scores)]]$log.p))
  
}

samps = get_best_fit(full)

pars = samps[[1]]
variances = samps[[2]]

# Imports the data and experimental conditions from the experiment
Multiple_Time_Series <- read.csv("Data_Rebound_4-2-25.csv")
exp_df <- read.csv("Rebound_experimental_conditions.csv")

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
    if(Sieve == 0){ Pred_N = Pred_N * naup_catch;}
    
    
  double birth_rate = b_M*exp(-comp_b / VOL * (c_N * N + c_J * J + A));

  double d_A_c = d_A*exp(comp_d / VOL * (c_N * N + c_J * J + A));
  double d_J_c = d_J*exp(comp_d / VOL * (c_N * N + c_J * J + A));
  double d_N_c = d_N*exp(comp_d / VOL * (c_N * N + c_J * J + A)) + cann*A;

  double m_N_c = m_N*exp(-comp_m / VOL * (c_N * N + c_J * J + A));
  double m_J_c = m_J*exp(-comp_m / VOL * (c_N * N + c_J * J + A));


  /* Specify adult rates */
  const double A_rates[2] = {d_A_c, Pred_A};
  double A_trans[2];
  reulermultinom(2, A, &A_rates, dt, &A_trans);
  
  /* Specify juvenile rates */
  const double J_rates[3] = {m_J_c, d_J_c, Pred_J};
  double J_trans[3];
  reulermultinom(3, J, &J_rates, dt, &J_trans);
  
  /* Specify nauplii rates */
  const double N_rates[3] = {m_N_c, d_N_c, Pred_N};
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
    Day = 0:119, #This needs to match the times in the actual dataset  (Generalize this)
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
  paramnames = names(pars),
  params = pars
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

# Simulate the model
params = list()
for(l in 1:32){
  for(i in (5:25)*10) #50:250)*1000
    params = c(params, list(full[[l]]$samples[i,]))
}

post_sim = function(x){simulate(GW_panel, shared=x, nsim=1)} 

sim = lapply(X=params,FUN=post_sim) 

# Create an empty list to hold data frames from each simulation
df_list <- vector("list", length(sim))

# Loop over each simulation (1-100)
for (s in seq_along(sim)) {
  sim_obj <- sim[[s]]@unit_objects
  
  # Loop over units in the simulation (1-60)
  unit_dfs <- lapply(seq_along(sim_obj), function(u) {
    unit_obj <- sim_obj[[u]]
    mat <- unit_obj@states  # 3 x 18 matrix (latent states - true abundance)
    
    # Create a data frame with 18 rows (one per time step aka week)
    data.frame(
      time = unit_obj@times,
      sim = s,
      unit = u,
      N = mat["N", ],
      J = mat["J", ],
      A = mat["A", ]
    )
  }) 
  
  # Combine all 60 units for this simulation into one data frame
  df_list[[s]] <- do.call(rbind, unit_dfs)
}


# Combine all simulations into one big data frame
full_df <- do.call(rbind, df_list)

#make new columns for AF1, JOA1,N1
full_df$N1 <- full_df$N
full_df$JOA1 <- round(full_df$J + 2/3*(full_df$A))
full_df$AF1 <- round(1/3*(full_df$A)) 

#add column for total copepods
full_df$total = full_df$AF1 + full_df$JOA1 + full_df$N1

#divide by 15 for per L
full_df$total = full_df$total

#identify covariates (treatments) for each unit
unit_covariates <- lapply(poList, function(unit) {
  covar_df <- unit@covar
})

covariates_to_extract <- c("Sieve", "Harvest_rate")
timestep <- 9

covariate_values <- lapply(unit_covariates, function(cov) {
  cov@table[covariates_to_extract, timestep]
})

cov_table_assignments <- do.call(rbind, covariate_values)
cov_table_assignments <- data.frame(unit = names(covariate_values), cov_table_assignments)

cov_table_assignments <- cov_table_assignments %>%
  mutate_at("unit", str_replace, "unit", "")

#bind the simulations and covariate (treatment) dfs
sims_and_covs = merge(full_df,cov_table_assignments,by="unit")
sims_and_covs$Sieve[sims_and_covs$Harvest_rate == "0"] <- "2" #controls 


sims_and_covs_summary = sims_and_covs %>% group_by(Sieve,Harvest_rate,time) %>% summarise(mean_A = mean(AF1), 
                                                                                     mean_J = mean(JOA1), 
                                                                                     mean_N = mean(N1),
                                                                                     mean_total = mean(total),
                                                                                     sd_A = sd(AF1),
                                                                                     sd_J = sd(JOA1),
                                                                                     sd_N = sd(N1),
                                                                                     sd_total = sd(total),
                                                                                     n_A = n(),
                                                                                     n_J = n(),
                                                                                     n_N = n(),
                                                                                     n_total = n(),
                                                                                     se_A = sd_A / sqrt(n_A),
                                                                                     se_J = sd_J / sqrt(n_J),
                                                                                     se_N = sd_N / sqrt(n_N),
                                                                                     se_total = sd_total / sqrt(n_total),
                                                                                     lower_ci_A = quantile(AF1, prob=0.025),
                                                                                     upper_ci_A = quantile(AF1, prob=0.975),
                                                                                     lower_ci_J = quantile(JOA1, prob=0.025),
                                                                                     upper_ci_J = quantile(JOA1, prob=0.975),
                                                                                     lower_ci_N =quantile(N1, prob=0.025),
                                                                                     upper_ci_N = quantile(N1, prob=0.975),
                                                                                     lower_ci_total =quantile(total, prob=0.025),
                                                                                     upper_ci_total = quantile(total, prob=0.975)) 
                                                                                     



#treatment groups: 1 - exp(-value) = the percent removed
#Sieve = 1 is 35
sims_and_covs_summary$Percent_Removed = 1-exp(-sims_and_covs_summary$Harvest_rate)
sims_and_covs_summary$Percent_Removed = as.factor(sims_and_covs_summary$Percent_Removed)


sims_and_covs_summary <- mutate(sims_and_covs_summary, Sieve = case_when(
  Sieve <= 0 ~ 100,
  Sieve == 1 ~ 35,
  Sieve == 2 ~ 0)) 

sims_and_covs_summary$treatment <- as.factor(paste(sims_and_covs_summary$Sieve,sims_and_covs_summary$Percent_Removed))


# #plot all 
ggplot(sims_and_covs_summary, aes(x = time, y = mean_A)) +
  geom_line() +
  facet_wrap(~ treatment) +
  theme_minimal() +
  labs(title = "Mean A over Time by Treatment Group") +
  geom_ribbon(aes(ymin=lower_ci_A, ymax=upper_ci_A, alpha =0.2),fill = "thistle2") + theme_classic()

ggplot(sims_and_covs_summary, aes(x = time, y = mean_J)) +
  geom_line() +
  facet_wrap(~ treatment) +
  theme_minimal() +
  labs(title = "Mean J over Time by Treatment Group") +
  geom_ribbon(aes(ymin=lower_ci_J, ymax=upper_ci_J, alpha =0.2),fill = "thistle2") + theme_classic()

ggplot(sims_and_covs_summary, aes(x = time, y = mean_N)) +
  geom_line() +
  facet_wrap(~ treatment) +
  theme_minimal() +
  labs(title = "Mean N over Time by Treatment Group") +
  geom_ribbon(aes(ymin=lower_ci_N, ymax=upper_ci_N, alpha =0.2),fill = "thistle2") + theme_classic()

ggplot(sims_and_covs_summary, aes(x = time, y = mean_total)) +
  geom_line() +
  facet_wrap(~ treatment) +
  theme_minimal() +
  labs(title = "Mean Total Copepods over Time by Treatment Group") +
  geom_ribbon(aes(ymin=lower_ci_total, ymax=upper_ci_total, alpha =0.2),fill = "thistle2") + theme_classic()



####Bring in file from super computer####
#The above code was crashing my computer, so we had to run on the cluster and save the file below for downstream analysis. 
sims_and_covs_summary = readRDS("sims_and_covs.RDA")

#bring in raw data (mesocosm data)
Multiple_Time_Series <- read.csv("Data_Rebound_4-2-25.csv")

real_df <- Multiple_Time_Series %>%
  rename(unit = Tank) %>%
  select(Day, unit, Sieve_Size,Percent_Removed, Nauplii_per_Tank,Juveniles_per_tank,Adults_per_Tank)

real_df$Sieve_Size[real_df$Percent_Removed == "0"] <- 0
real_df$Percent_Removed[real_df$Percent_Removed == "10"] <- 95
real_df$Percent_Removed[real_df$Percent_Removed == "20"] <- 98

real_df$treatment <- as.factor(paste(real_df$Sieve_Size,real_df$Percent_Removed/100))
real_df <- na.omit(real_df)
real_df$total <- real_df$Nauplii_per_Tank + real_df$Juveniles_per_tank + real_df$Adults_per_Tank
real_df$total <- real_df$total 

#calculate treatment means and SDs
real_df_summary = real_df %>% group_by(treatment,Day,Sieve_Size,Percent_Removed) %>% summarise(mean_A = mean(Adults_per_Tank), 
                                                                                          mean_J = mean(Juveniles_per_tank), 
                                                                                          mean_N = mean(Nauplii_per_Tank),
                                                                                          mean_total = mean(total),
                                                                                          sd_A = sd(Adults_per_Tank),
                                                                                          sd_J = sd(Juveniles_per_tank),
                                                                                          sd_N = sd(Nauplii_per_Tank),
                                                                                          sd_total = sd(total),
                                                                                          n_A = n(),
                                                                                          n_J = n(),
                                                                                          n_N = n(),
                                                                                          n_total = n(),
                                                                                          se_A = sd_A / sqrt(n_A),
                                                                                          se_J = sd_J / sqrt(n_J),
                                                                                          se_N = sd_N / sqrt(n_N),
                                                                                          se_total = sd_total / sqrt(n_total)) 


sims_and_covs_summary$treatment = gsub("100 0.1", "100 0.95", sims_and_covs_summary$treatment)
sims_and_covs_summary$treatment = gsub("100 0.2", "100 0.98", sims_and_covs_summary$treatment)
sims_and_covs_summary$treatment = gsub("35 0.1", "35 0.95", sims_and_covs_summary$treatment)
sims_and_covs_summary$treatment = gsub("35 0.2", "35 0.98", sims_and_covs_summary$treatment)

                                                                                        

#Edit labels for plotting 

treatment_labeller <- function(labels) {
  any_percent_values <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98)
  
  sapply(labels, function(x) {
    if (x == "0 0") {
      "control"
    } else {
      parts <- str_split_fixed(x, " ", 2)
      conc <- as.numeric(parts[1])
      perc <- as.numeric(parts[2])
      perc_rounded <- round(perc * 100)
      
      if (conc == 35 && perc %in% any_percent_values) {
        paste0("Size-unbiased ", perc_rounded, "%")
      } else if (conc == 100 && perc %in% any_percent_values) {
        paste0("Size-biased ", perc_rounded, "%")
      } else {
        NA  # Or handle differently if desired
      }
    }
  })
}

selected_treatments <- c(
  "0 0",
  "35 0.3", "35 0.5", "35 0.8", "35 0.98",
  "100 0.5", "100 0.8", "100 0.98"
)


sims_and_covs_summary_filtered = sims_and_covs_summary %>%  filter(treatment %in% selected_treatments)
real_df_summary_filtered = real_df_summary %>% filter(treatment %in% selected_treatments)

#Plot all for total copepods 

my_colors <- c(
  "0"   = "yellowgreen",
  "35" = "darkslateblue",
  "100"  = "cadetblue3"
)

#Plot all treatments 
ggplot(sims_and_covs_summary, aes(x = time, y = mean_total)) +    
  theme_minimal() +   
  geom_vline(xintercept = c(8, 36, 64), linetype = "dashed",
             color = "black", linewidth = 0.5, alpha = 0.8) +   
  geom_ribbon(aes(ymin=lower_ci_total, ymax=upper_ci_total), 
              fill = "grey", alpha = 0.8) +   
  theme_classic() +   
  geom_point(data = real_df_summary,
             aes(x = Day, y = mean_total, group = as.factor(Sieve_Size), color = as.factor(Sieve_Size))) +   
  geom_line(data = real_df_summary,
            aes(x = Day, y = mean_total, group = as.factor(Sieve_Size), color = as.factor(Sieve_Size)))  +    
  geom_errorbar(data = real_df_summary,
                aes(x = Day, ymin = mean_total - se_total, ymax = mean_total + se_total,
                    color = as.factor(Sieve_Size))) +   # make sure errorbars follow the same palette
  geom_line(color = "black") +   
  facet_wrap(~ treatment, labeller = labeller(treatment = treatment_labeller),
             nrow = 4, ncol = 5) +   
  labs(x = "Time (days)", y = "Total copepod abundance") +   
  theme_classic(base_size = 15) +   
  theme(legend.position = "none") + 
  scale_color_manual(values = my_colors) 


####Plot for main text####
ggplot(sims_and_covs_summary_filtered, aes(x = time, y = mean_total)) +    
  theme_minimal() +   
  geom_vline(xintercept = c(8, 36, 64), linetype = "dashed",
             color = "black", linewidth = 0.5, alpha = 0.8) +   
  geom_ribbon(aes(ymin=lower_ci_total, ymax=upper_ci_total), 
              fill = "grey", alpha = 0.8) +   
  theme_classic() +   
  geom_point(data = real_df_summary_filtered, size = 3,
             aes(x = Day, y = mean_total, group = as.factor(Sieve_Size), color = as.factor(Sieve_Size))) +   
  geom_line(data = real_df_summary_filtered, linewidth = 1.5,
            aes(x = Day, y = mean_total, group = as.factor(Sieve_Size), color = as.factor(Sieve_Size)))  +    
  geom_errorbar(data = real_df_summary_filtered, size = 1,
                aes(x = Day, ymin = mean_total - se_total, ymax = mean_total + se_total,
                    color = as.factor(Sieve_Size))) +   # make sure errorbars follow the same palette
  geom_line(color = "black", size = 1) +   
  facet_wrap(~ treatment, labeller = labeller(treatment = treatment_labeller),
             nrow = 2, ncol = 4) +   
  labs(x = "Time (days)", y = "Total copepod abundance") +   
  theme_classic(base_size = 20) +   
  theme(legend.position = "none") + 
  scale_color_manual(values = my_colors) 

 
####Plot for supplemental#### 

p1 = natparks.pals("BryceCanyon")
p2  = natparks.pals("Everglades")
p3 = natparks.pals("Cuyahoga")
natparks = c(p1,p2,p3)


#Need to do correction for juveniles and adults (AF, JOA) before plotting separately. 
adults = ggplot(sims_and_covs_summary, aes(x = time, y = mean_A)) +
  theme_minimal() +
  geom_ribbon(aes(ymin=lower_ci_A, ymax=upper_ci_A, alpha =0.2),fill = "olivedrab4") + theme_classic() +
  geom_point(data = real_df_summary, aes(x=Day,y=mean_A, color = treatment))+ 
  geom_errorbar(data=real_df_summary, aes(x=Day, ymin=mean_A-se_A, ymax=mean_A+se_A, color = treatment)) +
  geom_line() + labs(x = "Time (days)",y = "AF Abundance") + 
  facet_wrap(~ treatment,labeller = labeller(treatment = treatment_labeller)) + scale_color_manual(values = natparks)

juv = ggplot(sims_and_covs_summary, aes(x = time, y = mean_J)) +
  theme_minimal() +
  geom_ribbon(aes(ymin=lower_ci_J, ymax=upper_ci_J, alpha =0.2),fill = "olivedrab4") + theme_classic() +
  geom_point(data = real_df_summary, aes(x=Day,y=mean_J, color = treatment))+ 
  geom_errorbar(data=real_df_summary, aes(x=Day, ymin=mean_J-se_J, ymax=mean_J+se_J,color = treatment)) +
  geom_line() + labs(x = "Time (days)",y = "JOA Abundance") + 
  facet_wrap(~ treatment,labeller = labeller(treatment = treatment_labeller)) + scale_color_manual(values = natparks)


naup = ggplot(sims_and_covs_summary, aes(x = time, y = mean_N)) +
  theme_minimal() +
  geom_ribbon(aes(ymin=lower_ci_N, ymax=upper_ci_N, alpha =0.2),fill = "olivedrab4") + theme_classic() +
  geom_point(data = real_df_summary, aes(x=Day,y=mean_N, color = treatment))+ 
  geom_errorbar(data=real_df_summary, aes(x=Day, ymin=mean_N-se_N, ymax=mean_N+se_N, color = treatment)) +
  geom_line() + labs(x = "Time (days)",y = "N Abundance") + 
  facet_wrap(~ treatment,labeller = labeller(treatment = treatment_labeller)) + scale_color_manual(values = natparks)

library(ggpubr)
ggarrange(adults,juv,naup,
          nrow = 1,
          ncol =3,
          labels = c("A","B","C"), 
          legend = "none"
)







