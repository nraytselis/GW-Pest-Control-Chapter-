library(deSolve)
library(ggplot2)
library(tidyr)
library(ggsci)
library(deSolve)
library(mgcv)
library(itsadug)
library(cowplot)
library(stringr)
library(tidyverse)

setwd("~/Desktop/Rscripts/Data")

####Summary####
#This code runs the ODE for the the Pest Control Experiment with harvests/culling of copepod populations and plots how densities of infectious adults,
#non-infectious adults, and total copepods change across a removal gradient. It is parameterized with mcmc chains from our state space modeling 
#fits to mesocosm data. The mesocosm experiment involved culling a specific % of the copepod population every 28 days to 
#mimic the application of Abate pesticides. Control was either size-biased towards adults or equally applied to all stages. 

Harvest_ODE =function(t, y, parameters) {

  with(as.list(parameters),{
    N=y[1]; J=y[2]; A=y[3]; Es = y[4:(4+latent_stages - 1)]; I = y[4+latent_stages]
    VOL = 1

    d_A_c = d_A*exp(comp_d / VOL * (c_N * N + c_J * J + A + sum(Es) + I)) #density dependence in deaths

    d_J_c = d_J*exp(comp_d / VOL * (c_N * N + c_J * J + A + sum(Es) + I)) #density dependence in deaths

    d_N_c = d_N*exp(comp_d / VOL * (c_N * N + c_J * J + A + sum(Es) + I)) #density dependence in deaths
    
    m_N_c = m_N*exp(-comp_m/VOL*(c_N*N + c_J*J + A + sum(Es) + I)) #density dependence in maturation
    
    m_J_c = m_J*exp(-comp_m/VOL*(c_N*N + c_J*J + A + sum(Es) +I)) #density dependence in maturation

    dNdt = b_M*(A + sum(Es))/2*exp(-comp_b/VOL*(c_N*N + c_J*J + A + sum(Es) + I)) - (m_N_c+d_N_c)*N - cann*(A + I + sum(Es))*N

    dJdt = m_N_c*N - (m_J_c+d_J_c)*J

    dAdt = m_J_c*J - d_A_c*A - lambda*A

    # development of all stages
    latent_progression = latent_rate*Es
          # lost to next stage   #death      #gained from last stage
    dEsdt = -latent_progression - d_A_c*Es + c(lambda*A, latent_progression[1:(latent_stages - 1)])

    dIdt = as.numeric(latent_progression[latent_stages]) - d_A_c*I

    result = c(dNdt,dJdt,dAdt, dEsdt, dIdt)

    return(list(result))
  }
  )
}



Initial_conditions = c(N = 7500, J = 6000, A = 700)/15 #abundances divided by the volume of the tanks
timespan = 98 #don't simulate last 4 weeks

#ReboundParams <- readRDS(file = "Rebound_parameters_full_100000p8.RDA") #bring in best parameters

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

full = get_best_fit(full)

pars = full[[1]]

#trying to achieve 5% prevalence: Garrett, K.B., Box, E.K., Cleveland, C.A. et al. 
#Dogs and the classic route of Guinea Worm transmission: an evaluation of copepod ingestion. 
#Sci Rep 10, 1430 (2020). https://doi.org/10.1038/s41598-020-58191-4

pars["lambda"] = 0.0034 #5% infection among adults
pars["latent_stages"] = 60
pars["latent_rate"] = 4.3
# parameters["comp_M"] = 0
# parameters["comp_b"] = 0
# parameters["comp_d"] = 0
# parameters["can"] = 0.0001


parameters = signif(pars,3)

Exposed_names = paste0("E", 1:parameters["latent_stages"])
Exposed_values = rep(0, times=parameters["latent_stages"])
names(Exposed_values) = Exposed_names
Exposed_values
#
Initial_conditions = c(N = 7500, J = 6000, A = 700, Exposed_values, I = 0)/15
#
ReboundSim = data.frame(ode(y = Initial_conditions, times=1:timespan, parms=parameters, hmax=1,
                method="lsoda", func=Harvest_ODE))
ReboundSim[,"Es"] = rowSums(ReboundSim) - ReboundSim[,"N"] - ReboundSim[,"J"]- ReboundSim[,"A"] - ReboundSim[,"I"] - ReboundSim[,"time"]
# plot(ReboundSim)
# par(mfrow=c(3, 2))
# plot(N ~ time, data=ReboundSim, typ="l")
# plot(J ~time, data=ReboundSim, typ="l")
# plot(A ~time, data=ReboundSim, typ="l")
# plot(Es ~ time, data=ReboundSim, typ="l")
# plot(I ~ time, data=ReboundSim, typ="l")
#

#generalize event data function
Harvest_scheme <- function(states = names(Initial_conditions), harvest_days=c(8,36,64),harvest_prop=0,harvest_nauplii=TRUE){
  var = rep(states,times = length(harvest_days))
  time = rep(harvest_days,each=length(states))
  value = 1 - harvest_prop
  event.data = data.frame(var,time,value,"method"="multiply")
  if(harvest_nauplii == FALSE){
    event.data = subset(event.data,var != "N")
  }
  event.data
}

event_data = Harvest_scheme(harvest_nauplii = TRUE)

ReboundSim = data.frame(ode(y = Initial_conditions, times=1:timespan, parms=parameters, hmax=1,
                            method="lsoda", func=Harvest_ODE,events=list(data=event_data)))

ReboundSim[,"Es"] = rowSums(ReboundSim) - ReboundSim[,"N"] - ReboundSim[,"J"]- ReboundSim[,"A"] - ReboundSim[,"I"] - ReboundSim[,"time"]
#plot(ReboundSim)
#par(mfrow=c(3, 2))
#plot(N ~ time, data=ReboundSim, typ="l")
#plot(J ~time, data=ReboundSim, typ="l")
#plot(A ~time, data=ReboundSim, typ="l")
#plot(Es ~ time, data=ReboundSim, typ="l")
#plot(I ~ time, data=ReboundSim, typ="l")



#Infecteds_over_time = ggplot(data=ReboundSim,aes(x=time,y=I)) + geom_line() + theme_minimal()
#Infecteds_over_time
#Adults_over_time = ggplot(data=ReboundSim,aes(x=time,y=A)) + geom_line() + theme_minimal()
#plot_grid(Adults_over_time,Infecteds_over_time)


ReboundSim = data.frame(ReboundSim)
ReboundSim[,"Exposed"] = apply(X=ReboundSim[,which(str_detect(colnames(ReboundSim),"E"))],MARGIN=1,FUN = sum)

ReboundSim = ReboundSim %>% select(time,N,J,A,Exposed,I)
ReboundSim = ReboundSim %>% pivot_longer(cols = c(N,J,A,Exposed,I))

#ggplot(ReboundSim, aes(x=time, y = value + 0.01 , group = name, color = name)) + geom_line() + scale_y_log10()



#####
Initial_conditions = c(N = 7500, J = 6000, A = 700, Exposed_values, I = 0)/15
timespan = 98 #do not simulate last 4 weeks
harvest_rates = seq(0,0.9,0.1)
harvest_nauplii = c(TRUE,FALSE)
harvest_params = expand.grid(hr = harvest_rates,hn = harvest_nauplii)
colnames(harvest_params) = c("harvest_rate","harvest_nauplii")
mean_N = numeric()
mean_J = numeric()
mean_A = numeric()
mean_I = numeric()
mean_E = numeric()
mean_total = numeric()
mean_NIA = numeric()
mean_AJ = numeric()

for(i in 1:nrow(harvest_params)){
  event_data = Harvest_scheme(harvest_prop = harvest_params[i,"harvest_rate"],harvest_nauplii = harvest_params[i,"harvest_nauplii"])
  ReboundSim = ode(y = Initial_conditions, times=1:timespan, parms=parameters,
                   method="lsoda", func=Harvest_ODE,events=list(data=event_data))
  mean_N[i] = mean(ReboundSim[,"N"])
  mean_J[i] = mean(ReboundSim[,"J"])
  mean_A[i] = mean(ReboundSim[,"A"])
  mean_I[i] = mean(ReboundSim[,"I"])
  mean_E[i] = mean(apply(X=ReboundSim[,which(str_detect(colnames(ReboundSim),"E"))],MARGIN=1,FUN = sum))
 
  E_time = apply(ReboundSim[, str_detect(colnames(ReboundSim), "E")], 1, sum)
  total_time = ReboundSim[,"N"] + ReboundSim[,"J"] + ReboundSim[,"A"] + ReboundSim[,"I"] + E_time
  NIA_time = ReboundSim[,"A"] + E_time
  AJ_time = ReboundSim[,"J"] + ReboundSim[,"A"] + ReboundSim[,"I"] + E_time
  
  #mean total and se total all stages 
  mean_total[i] = mean(total_time)
  mean_NIA[i] = mean(NIA_time)
  mean_AJ[i] = mean(AJ_time)
}
sim_output = cbind(harvest_params,mean_N,mean_J,mean_A, mean_I,mean_E,mean_total, mean_NIA, mean_AJ)

sim_output = data.frame(sim_output)


# Plot all on one graph 
sim_long_mean <- sim_output %>%
  select(harvest_rate, harvest_nauplii,
         mean_total, mean_I, mean_NIA) %>%
  pivot_longer(
    cols = c(mean_total,mean_I, mean_NIA),
    names_to = "group",
    values_to = "mean"
  ) %>%
  mutate(group_label = dplyr::recode(group,
                                     "mean_total" = "Total Copepods",
                                     "mean_I" = "Infecteds",
                                     "mean_NIA" = "NIA"))


sim_long_mean$harvest_nauplii = as.factor(sim_long_mean$harvest_nauplii) 

#5.5% infection prevelance among adults and copepodites 
#> 50/910 = 0.05494505, lambda = 0.032

ggplot(sim_long_mean, aes(x = as.factor(harvest_rate), y = mean, 
                     color = group, shape = as.factor(harvest_nauplii), group = interaction(group, harvest_nauplii),linetype=harvest_nauplii)) +
  geom_line(position = position_dodge(width = 0.5), linewidth = 2) +
  theme_classic() +
  labs(
    x = "Removal Rate",
    y = expression('Mean Density, L' ^ -1),
    color = "Category",
    shape = "Harvest Nauplii",
    color = "Stage",
    linetype = "Mortality Type"
  ) + theme(text = element_text(size = 20)) +
  scale_color_manual(
    values = c(mean_NIA = "yellowgreen", mean_I = "darkslateblue", mean_total = "cadetblue3"),
    labels = c(
      mean_NIA = "Non-infectious Adults",
      mean_I = "Infectious Adults",
      mean_total = "Total (all stages)"
    )
  ) + scale_linetype_discrete(name = "Mortality Type", labels = c("Size-biased", "Size-unbiased")) +
  theme(legend.key.size = unit(0.5, "cm")) + scale_y_log10() +
  theme(legend.position = c(0.8, 0.4))

