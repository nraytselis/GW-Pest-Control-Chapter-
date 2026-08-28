library(deSolve) 
library(ggplot2) 
library(tidyr)
library(ggsci)
library(deSolve)
library(mgcv)
library(itsadug)
library(cowplot)
library(wesanderson)
library(stringr)
library(dplyr)
library(ggpubr)
setwd("~/Desktop/Rscripts/Data")


####Summary####
#Plotting density of different copepod stages over time based on Pest Control ODE using parameters from best fit of mcmc chains. 
#These plots show that although total density and density of other size classes show signs of compensation/overcompensation, infectious copepods 
#do not recover. 



Harvest_ODE =function(t, y, parameters) {
  
  with(as.list(parameters),{
    N=y[1]; J=y[2]; A=y[3]; Es = y[4:(4+latent_stages - 1)]; I = y[4+latent_stages]; EJs = y[(5+latent_stages):((5+2*latent_stages - 1))]; IJ = y[5+2*latent_stages]
    VOL = 1
    
    d_A_c = d_A*exp(comp_d / VOL * (c_N * N + c_J * (J + sum(EJs)+IJ) + A + sum(Es)  + I)) #density dependence in deaths
    
    d_J_c = d_J*exp(comp_d / VOL * (c_N * N + c_J * (J + sum(EJs)+IJ) + A + sum(Es) + I)) #density dependence in deaths
    
    d_N_c = d_N*exp(comp_d / VOL * (c_N * N + c_J * (J + sum(EJs)+IJ) + A + sum(Es) +  I)) #density dependence in deaths
    
    m_N_c = m_N*exp(-comp_m/VOL*(c_N*N + c_J*(J + sum(EJs)+IJ) + A + sum(Es) + I)) #density dependence in maturation
    
    m_J_c = m_J*exp(-comp_m/VOL*(c_N*N + c_J*(J + sum(EJs)+IJ) + A + sum(Es) +  I)) #density dependence in maturation
    
    dNdt = b_M*(A + sum(Es))/2*exp(-comp_b/VOL*(c_N*N + c_J*(J + sum(EJs)+IJ) + A + sum(Es) + I)) - (m_N_c+d_N_c)*N - cann*(A + I + sum(Es))*N
    
    dJdt = m_N_c*N - (m_J_c+d_J_c)*J - lambda*J
    
    dAdt = m_J_c*J - d_A_c*A - lambda*A
    
    # development of all stages
    latent_progression = latent_rate*Es # This one is for adults
    latent_progressionJs = latent_rate*EJs # This one is for juveniles
    # development of exposed juveniles into exposed adults
    exposed_development = m_J_c*EJs
    
              # lost to next stage  #gained from juveniles   #death       #gained from last stage
    dEsdt = -latent_progression + exposed_development - d_A_c*Es + c(lambda*A, latent_progression[1:(latent_stages - 1)])
    # lost to next stage   #developed to adult            #death      #gained from last stage
    dEJsdt = -latent_progressionJs - exposed_development - d_J_c*EJs + c(lambda*J, latent_progressionJs[1:(latent_stages - 1)])
    
  
    
    dIdt = as.numeric(latent_progression[latent_stages]) - d_A_c*I
    dIJdt = as.numeric(latent_progressionJs[latent_stages]) - d_J_c*IJ
    
    result = c(dNdt,dJdt,dAdt, dEsdt, dIdt, dEJsdt, dIJdt)
    
    return(list(result))
  }
  )
}



Initial_conditions = c(N = 7500, J = 6000, A = 700)/15 #abundances divided by the volume of the tanks
timespan = 100 #don't simulate last 4 weeks

# Dave is just bringing in best fit parameter file, rather than reading all of the mcmc chains
#parameters = readRDS("C:/RData/bestfit.rds")

#bring in mcmc chains
fullA = readRDS("Rebound_parameters_full_disp250kff_5.RDA")
fullB = readRDS("Rebound_parameters_full_disp250kff_3.RDA")
fullC = readRDS("Rebound_parameters_full_disp250kff_05.RDA")
fullD = readRDS("Rebound_parameters_full_disp250kff_01.RDA")
ReboundParams = c(fullA, fullB, fullC, fullD)


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

 ReboundParams = get_best_fit(ReboundParams)
 parameters = as.list(ReboundParams[[1]])

#define other parameter values 
parameters["latent_stages"] = 60
parameters["latent_rate"] = 4.3
parameters["lambda"] = 0.032
parameters = unlist(parameters)
parameters = signif(parameters,3) 

#define exposed classes (seperate stages because of linear chain trick)
Exposed_names = paste0("E", 1:parameters["latent_stages"])
Exposed_values = rep(0, times=parameters["latent_stages"])
names(Exposed_values) = Exposed_names
Exposed_values

#define exposed classes (seperate stages because of linear chain trick)
ExposedJ_names = paste0("EJ", 1:parameters["latent_stages"])
ExposedJ_values = rep(0, times=parameters["latent_stages"])
names(ExposedJ_values) = ExposedJ_names
ExposedJ_values

Initial_conditions = c(N = 7500, J = 6000, A = 700, Exposed_values, I = 0, ExposedJ_values, IJ = 0)/15

ReboundSim = data.frame(ode(y = Initial_conditions, times=1:timespan, parms=parameters, hmax=1,
                            method="lsoda", func=Harvest_ODE))

# Hard coded this to work for the 60 latent stages
ReboundSim[,"Es"] = rowSums(ReboundSim[,4:63])#rowSums(ReboundSim) - ReboundSim[,"N"] - ReboundSim[,"J"]- ReboundSim[,"A"] - ReboundSim[,"I"] - ReboundSim[,"time"]
ReboundSim[,"EJs"] = rowSums(ReboundSim[,65:124])#rowSums(ReboundSim) - ReboundSim[,"N"] - ReboundSim[,"J"]- ReboundSim[,"A"] - ReboundSim[,"I"] - ReboundSim[,"time"]

ggplot(ReboundSim) + geom_line(aes(x=time, y=N)) + 
  geom_line(aes(x=time, y=J), color="blue") +
  geom_line(aes(x=time, y=A), color="red") +
  geom_line(aes(x=time, y=Es), color="red", linetype = 2) +
  geom_line(aes(x=time, y=EJs), color="blue", linetype = 2) +
  geom_line(aes(x=time, y=I), color="red", linetype = 3) +
  geom_line(aes(x=time, y=IJ), color="blue", linetype = 3)

##### Dave stopped making edits here#####







#generalize event data function, change harvest prop depending on % harvested 
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
ReboundSim = data.frame(ReboundSim)
ReboundSim[,"Exposed"] = apply(X=ReboundSim[,which(str_detect(colnames(ReboundSim),"E"))],MARGIN=1,FUN = sum) 
ReboundSim = ReboundSim %>% select(time,N,J,A,Exposed,I,IJ)
ReboundSim = ReboundSim %>% pivot_longer(cols = c(N,J,A,Exposed,I,IJ)) 


#90%
event_data_90 = Harvest_scheme(harvest_nauplii = TRUE,harvest_prop = 0.9)
ReboundSim_90 = data.frame(ode(y = Initial_conditions, times=1:timespan, parms=parameters, hmax=1,
                            method="lsoda", func=Harvest_ODE,events=list(data=event_data_90)))
ReboundSim_90 = data.frame(ReboundSim_90)
ReboundSim_90[,"Exposed"] = apply(X=ReboundSim_90[,which(str_detect(colnames(ReboundSim_90),"E"))],MARGIN=1,FUN = sum) 
ReboundSim_90 = ReboundSim_90 %>% select(time,N,J,A,Exposed,I,IJ)
ReboundSim_90 = ReboundSim_90 %>% pivot_longer(cols = c(N,J,A,Exposed,I,IJ)) 


#50%
event_data_50 = Harvest_scheme(harvest_nauplii = TRUE,harvest_prop = 0.5)
ReboundSim_50 = data.frame(ode(y = Initial_conditions, times=1:timespan, parms=parameters, hmax=1,
                               method="lsoda", func=Harvest_ODE,events=list(data=event_data_50)))
ReboundSim_50 = data.frame(ReboundSim_50)
ReboundSim_50[,"Exposed"] = apply(X=ReboundSim_50[,which(str_detect(colnames(ReboundSim_50),"E"))],MARGIN=1,FUN = sum) 
ReboundSim_50 = ReboundSim_50 %>% select(time,N,J,A,Exposed,I,IJ)
ReboundSim_50 = ReboundSim_50 %>% pivot_longer(cols = c(N,J,A,Exposed,I,IJ)) 


#control
event_data = Harvest_scheme(harvest_nauplii = FALSE)

ReboundSim = data.frame(ode(y = Initial_conditions, times=1:timespan, parms=parameters, hmax=1,
                            method="lsoda", func=Harvest_ODE,events=list(data=event_data)))
ReboundSim = data.frame(ReboundSim)
ReboundSim[,"Exposed"] = apply(X=ReboundSim[,which(str_detect(colnames(ReboundSim),"E"))],MARGIN=1,FUN = sum) 
ReboundSim = ReboundSim %>% select(time,N,J,A,Exposed,I,IJ)
ReboundSim = ReboundSim %>% pivot_longer(cols = c(N,J,A,Exposed,I,IJ)) 

#90%
event_data_90 = Harvest_scheme(harvest_nauplii = FALSE,harvest_prop = 0.9)
ReboundSim_90 = data.frame(ode(y = Initial_conditions, times=1:timespan, parms=parameters, hmax=1,
                               method="lsoda", func=Harvest_ODE,events=list(data=event_data_90)))
ReboundSim_90 = data.frame(ReboundSim_90)
ReboundSim_90[,"Exposed"] = apply(X=ReboundSim_90[,which(str_detect(colnames(ReboundSim_90),"E"))],MARGIN=1,FUN = sum) 
ReboundSim_90 = ReboundSim_90 %>% select(time,N,J,A,Exposed,I,IJ)
ReboundSim_90 = ReboundSim_90 %>% pivot_longer(cols = c(N,J,A,Exposed,I,IJ)) 

#50%
event_data_50 = Harvest_scheme(harvest_nauplii = FALSE,harvest_prop = 0.5)
ReboundSim_50 = data.frame(ode(y = Initial_conditions, times=1:timespan, parms=parameters, hmax=1,
                               method="lsoda", func=Harvest_ODE,events=list(data=event_data_50)))
ReboundSim_50 = data.frame(ReboundSim_50)
ReboundSim_50[,"Exposed"] = apply(X=ReboundSim_50[,which(str_detect(colnames(ReboundSim_50),"E"))],MARGIN=1,FUN = sum) 
ReboundSim_50 = ReboundSim_50 %>% select(time,N,J,A,Exposed,I,IJ)
ReboundSim_50 = ReboundSim_50 %>% pivot_longer(cols = c(N,J,A,Exposed,I,IJ)) 

#combined plot 
sim_90_data <- as.data.frame(sim_90[["data"]]) 
sim_50_data <- as.data.frame(sim_50[["data"]]) 
control_data <- as.data.frame(control[["data"]]) 

control_data$treatment <- "control"
sim_50_data$treatment <- "50% with nauplii harvest"
sim_90_data$treatment <- "90% with nauplii harvest"

combined_df <- bind_rows(control_data, sim_90_data,sim_50_data)

sim_90_noNaup_data <- as.data.frame(sim_90_noNaup[["data"]]) 
sim_50_noNaup_data <- as.data.frame(sim_50_noNaup[["data"]]) 
control_noNaup_data <- as.data.frame(control_noNaup[["data"]]) 

control_noNaup_data$treatment <- "control"
sim_50_noNaup_data$treatment <- "50% no nauplii harvest"
sim_90_noNaup_data$treatment <- "90% no nauplii harvest"

combined_noNaup_df <- bind_rows(control_noNaup_data, sim_90_noNaup_data,sim_50_noNaup_data)

combined_noNaup_df$name <- factor(combined_noNaup_df$name)

combined_noNaup_df$treatment <- gsub(" no nauplii harvest", "", combined_noNaup_df$treatment)

combined_df$treatment <- gsub(" with nauplii harvest", "", combined_df$treatment)

# Filter out 50% treatment from both datasets
combined_df <- combined_df %>% filter(treatment != "50%")
combined_noNaup_df <- combined_noNaup_df %>% filter(treatment != "50%")

#combine the adults and exposed
#Pivot wider to add columns
combined_df_wide = combined_df %>% pivot_wider(names_from = name, values_from = value)
combined_df_wide$NIA <- combined_df_wide$A + combined_df_wide$Exposed
#remove A and E columns
combined_df_wide$A <- NULL
combined_df_wide$Exposed <- NULL 
#pivot longer   
combined_df = combined_df_wide %>% pivot_longer(cols=c("N","J","I","NIA","IJ"),names_to = "name",values_to = "value")
  
combined_noNaup_df_wide = combined_noNaup_df %>% pivot_wider(names_from = name, values_from = value)
combined_noNaup_df_wide$NIA <- combined_noNaup_df_wide$A + combined_noNaup_df_wide$Exposed
#remove A and E columns
combined_noNaup_df_wide$A <- NULL
combined_noNaup_df_wide$Exposed <- NULL 
#pivot longer   
combined_noNaup_df = combined_noNaup_df_wide %>% pivot_longer(cols=c("N","J","I","NIA","IJ"),names_to = "name",values_to = "value")

# Clean treatment labels (remove suffixes)
combined_df$treatment <- gsub(" with nauplii harvest", "", combined_df$treatment)
combined_noNaup_df$treatment <- gsub(" no nauplii harvest", "", combined_noNaup_df$treatment)

#add additional selective/no selective column
combined_df = combined_df %>% mutate(mortalitytype = "Complete")
combined_noNaup_df = combined_noNaup_df %>% mutate(mortalitytype = "Selective")

#combine two dataframes
combined_df = bind_rows(combined_df,combined_noNaup_df) 

#combine mortalitytype and treatment columns
combined_df$treatment <- as.factor(paste(combined_df$treatment,combined_df$mortalitytype))
combined_df$treatment = gsub("control Complete", "control", combined_df$treatment )
combined_df$treatment = gsub("control Selective", "control", combined_df$treatment )

combined_df_I = combined_df %>% filter(name == "I")
combined_df_J = combined_df %>% filter(name == "J")
combined_df_N = combined_df %>% filter(name == "N")
combined_df_NIA = combined_df %>% filter(name == "NIA")
combined_df_IJ = combined_df %>% filter(name == "IJ")

combined_df_I <- combined_df_I %>%
  mutate(treat_group = case_when(
    str_detect(treatment, regex("control", ignore_case = TRUE))   ~ "control",
    str_detect(treatment, regex("complete", ignore_case = TRUE))  ~ "complete",
    str_detect(treatment, regex("selective", ignore_case = TRUE)) ~ "selective",
    TRUE ~ NA_character_   # fallback if nothing matches
  ))

combined_df_J <- combined_df_J %>%
  mutate(treat_group = case_when(
    str_detect(treatment, regex("control", ignore_case = TRUE))   ~ "control",
    str_detect(treatment, regex("complete", ignore_case = TRUE))  ~ "complete",
    str_detect(treatment, regex("selective", ignore_case = TRUE)) ~ "selective",
    TRUE ~ NA_character_
  ))

combined_df_N <- combined_df_N %>%
  mutate(treat_group = case_when(
    str_detect(treatment, regex("control", ignore_case = TRUE))   ~ "control",
    str_detect(treatment, regex("complete", ignore_case = TRUE))  ~ "complete",
    str_detect(treatment, regex("selective", ignore_case = TRUE)) ~ "selective",
    TRUE ~ NA_character_
  ))

combined_df_NIA <- combined_df_NIA %>%
  mutate(treat_group = case_when(
    str_detect(treatment, regex("control", ignore_case = TRUE))   ~ "control",
    str_detect(treatment, regex("complete", ignore_case = TRUE))  ~ "complete",
    str_detect(treatment, regex("selective", ignore_case = TRUE)) ~ "selective",
    TRUE ~ NA_character_
  ))

combined_df_IJ <- combined_df_IJ %>%
  mutate(treat_group = case_when(
    str_detect(treatment, regex("control", ignore_case = TRUE))   ~ "control",
    str_detect(treatment, regex("complete", ignore_case = TRUE))  ~ "complete",
    str_detect(treatment, regex("selective", ignore_case = TRUE)) ~ "selective",
    TRUE ~ NA_character_
  ))

combined_df_I <- combined_df_I %>%
  mutate(treat_group = factor(treat_group, 
                              levels = c("control", "complete","selective")))

combined_df_J <- combined_df_J %>%
  mutate(treat_group = factor(treat_group, 
                              levels = c("control", "complete","selective")))

combined_df_N <- combined_df_N %>%
  mutate(treat_group = factor(treat_group, 
                              levels = c("control", "complete","selective")))

combined_df_NIA <- combined_df_NIA %>%
  mutate(treat_group = factor(treat_group, 
                              levels = c("control", "complete","selective")))

combined_df_IJ <- combined_df_IJ %>%
  mutate(treat_group = factor(treat_group, 
                              levels = c("control", "complete","selective")))

####Publication Plots####

my_colors <- c(
  control   = "yellowgreen",
  selective = "darkslateblue",
  complete  ="cadetblue3"
)

I = ggplot(combined_df_I, aes(x = time, y = value + 0.01, color = treat_group)) + geom_vline(xintercept = c(8, 36, 64), linetype = "dashed", color = "black", linewidth = 1, alpha = 0.8) +
  geom_line(linewidth=1.5) +
  labs(x = "Time (days)", y = "Infectious Adults \nper L", color = "Removal Type") +
  theme_classic(base_size = 15) +
  scale_color_manual(values = my_colors, labels = c("Control", "Size-unbiased","Size-biased")) 

J = ggplot(combined_df_J, aes(x = time, y = value + 0.01, color = treat_group)) + geom_vline(xintercept = c(8, 36, 64), linetype = "dashed", color = "black", linewidth = 1, alpha = 0.8) +
  geom_line(linewidth=1.5) +
  labs(x = "Time (days)", y = "Juveniles \nper L", color = "Removal Type") +
  theme_classic(base_size = 15) +
  scale_color_manual(values = my_colors,labels = c("Control", "Size-unbiased","Size-biased")) + theme(axis.title.x = element_blank())

N = ggplot(combined_df_N, aes(x = time, y = value + 0.01, color = treat_group)) + geom_vline(xintercept = c(8, 36, 64), linetype = "dashed", color = "black", linewidth = 1, alpha = 0.8) +
  geom_line(linewidth=1.5) +
  labs(x = "Time (days)", y = "Nauplii \nper L", color = "Removal Type") +
  theme_classic(base_size = 15) +
  scale_color_manual(values = my_colors,labels = c("Control", "Size-unbiased","Size-biased")) + theme(axis.title.x = element_blank())

NIA = ggplot(combined_df_NIA, aes(x = time, y = value + 0.01, color = treat_group)) + geom_vline(xintercept = c(8, 36, 64), linetype = "dashed", color = "black", linewidth = 1, alpha = 0.8) +
  geom_line(linewidth=1.5) +
  labs(x = "Time (days)", y = "Non-Infected Adults \nper L", color = "Removal Type") +
  theme_classic(base_size = 15) +
  scale_color_manual(values = my_colors,labels = c("Control", "Size-unbiased","Size-biased")) + theme(axis.title.x = element_blank())

IJ = ggplot(combined_df_IJ, aes(x = time, y = value + 0.01, color = treat_group)) + geom_vline(xintercept = c(8, 36, 64), linetype = "dashed", color = "black", linewidth = 1, alpha = 0.8) +
  geom_line(linewidth=1.5) +
  labs(x = "Time (days)", y = "Infectious Juveniles \nper L", color = "Removal Type") +
  theme_classic(base_size = 15) +
  scale_color_manual(values = my_colors,labels = c("Control", "Size-unbiased","Size-biased")) + theme(axis.title.x = element_blank())


library(ggpubr)

ggarrange(
  I,
  IJ,
  J,
  N,
  NIA,        
  common.legend = TRUE,    # Share a common legend
  legend = "right",       # Put the legend at the bottom
  align = "h",             # Align plots horizontally (x-axis alignment)
  nrow = 5                 # Arrange in 2 rows
)


