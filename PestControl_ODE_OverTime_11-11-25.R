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

setwd("~/Desktop/Rscripts/Data")


####Summary####
#Plotting density of different copepod stages over time based on Pest Control ODE using parameters from best fit of mcmc chains. 
#These plots show that although total density and density of other size classes show signs of compensation/overcompensation, infectious copepods 
#do not recover. 



Harvest_ODE =function(t, y, parameters) {
  
  with(as.list(parameters),{
    N=y[1]; J=y[2]; A=y[3]; Es = y[4:(4+latent_stages - 1)]; I = y[4+latent_stages]
    VOL = 1
    
    d_A_c = d_A*exp(comp_d / VOL * (c_N * N + c_J * J + A + sum(Es))) #density dependence in deaths
    
    d_J_c = d_J*exp(comp_d / VOL * (c_N * N + c_J * J + A + sum(Es))) #density dependence in deaths
    
    d_N_c = d_N*exp(comp_d / VOL * (c_N * N + c_J * J + A + sum(Es))) #density dependence in deaths
    
    m_N_c = m_N*exp(-comp_m/VOL*(c_N*N + c_J*J + A + sum(Es))) #density dependence in maturation
    
    m_J_c = m_J*exp(-comp_m/VOL*(c_N*N + c_J*J + A + sum(Es))) #density dependence in maturation
    
    dNdt = b_M*(A + sum(Es))/2*exp(-comp_b/VOL*(c_N*N + c_J*J + A + sum(Es))) - (m_N_c+d_N_c)*N - cann*(A + I + sum(Es))*N
    
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
timespan = 100 #don't simulate last 4 weeks

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

Initial_conditions = c(N = 7500, J = 6000, A = 700, Exposed_values, I = 0)/15

ReboundSim = data.frame(ode(y = Initial_conditions, times=1:timespan, parms=parameters, hmax=1,
                            method="lsoda", func=Harvest_ODE))

ReboundSim[,"Es"] = rowSums(ReboundSim) - ReboundSim[,"N"] - ReboundSim[,"J"]- ReboundSim[,"A"] - ReboundSim[,"I"] - ReboundSim[,"time"]

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
ReboundSim = ReboundSim %>% select(time,N,J,A,Exposed,I)
ReboundSim = ReboundSim %>% pivot_longer(cols = c(N,J,A,Exposed,I)) 
control = ggplot(ReboundSim, aes(x=time, y = value + 0.01 , group = name, color = name)) + 
  geom_line() + scale_y_log10() + ylab("density per L") + theme_minimal()
control

#90%
event_data_90 = Harvest_scheme(harvest_nauplii = TRUE,harvest_prop = 0.9)
ReboundSim_90 = data.frame(ode(y = Initial_conditions, times=1:timespan, parms=parameters, hmax=1,
                            method="lsoda", func=Harvest_ODE,events=list(data=event_data_90)))
ReboundSim_90 = data.frame(ReboundSim_90)
ReboundSim_90[,"Exposed"] = apply(X=ReboundSim_90[,which(str_detect(colnames(ReboundSim_90),"E"))],MARGIN=1,FUN = sum) 
ReboundSim_90 = ReboundSim_90 %>% select(time,N,J,A,Exposed,I)
ReboundSim_90 = ReboundSim_90 %>% pivot_longer(cols = c(N,J,A,Exposed,I)) 
sim_90 = ggplot(ReboundSim_90, aes(x=time, y = value + 0.01 , group = name, color = name)) + geom_line() + scale_y_log10()
#sim_90

#50%
event_data_50 = Harvest_scheme(harvest_nauplii = TRUE,harvest_prop = 0.5)
ReboundSim_50 = data.frame(ode(y = Initial_conditions, times=1:timespan, parms=parameters, hmax=1,
                               method="lsoda", func=Harvest_ODE,events=list(data=event_data_50)))
ReboundSim_50 = data.frame(ReboundSim_50)
ReboundSim_50[,"Exposed"] = apply(X=ReboundSim_50[,which(str_detect(colnames(ReboundSim_50),"E"))],MARGIN=1,FUN = sum) 
ReboundSim_50 = ReboundSim_50 %>% select(time,N,J,A,Exposed,I)
ReboundSim_50 = ReboundSim_50 %>% pivot_longer(cols = c(N,J,A,Exposed,I)) 
sim_50 = ggplot(ReboundSim_50, aes(x=time, y = value + 0.01 , group = name, color = name)) + geom_line() + scale_y_log10()
#sim_50


#control, 50, 90 for harvest_nauplii_true and lambda = 1.3
plot_grid(control,sim_50,sim_90,labels=c('control', '50%','90%')) + draw_label("Log densities/L when harvest_nauplii = TRUE and lambda =1.3", fontface='bold')


#same thing but with no nauplii harvest 
#(override names so need to re run sim above to get graphs with nauplii harvest again)

#control
event_data = Harvest_scheme(harvest_nauplii = FALSE)

ReboundSim = data.frame(ode(y = Initial_conditions, times=1:timespan, parms=parameters, hmax=1,
                            method="lsoda", func=Harvest_ODE,events=list(data=event_data)))
ReboundSim = data.frame(ReboundSim)
ReboundSim[,"Exposed"] = apply(X=ReboundSim[,which(str_detect(colnames(ReboundSim),"E"))],MARGIN=1,FUN = sum) 
ReboundSim = ReboundSim %>% select(time,N,J,A,Exposed,I)
ReboundSim = ReboundSim %>% pivot_longer(cols = c(N,J,A,Exposed,I)) 
control_noNaup = ggplot(ReboundSim, aes(x=time, y = value + 0.01 , group = name, color = name)) + geom_line() + scale_y_log10()
#control_noNaup
#90%
event_data_90 = Harvest_scheme(harvest_nauplii = FALSE,harvest_prop = 0.9)
ReboundSim_90 = data.frame(ode(y = Initial_conditions, times=1:timespan, parms=parameters, hmax=1,
                               method="lsoda", func=Harvest_ODE,events=list(data=event_data_90)))
ReboundSim_90 = data.frame(ReboundSim_90)
ReboundSim_90[,"Exposed"] = apply(X=ReboundSim_90[,which(str_detect(colnames(ReboundSim_90),"E"))],MARGIN=1,FUN = sum) 
ReboundSim_90 = ReboundSim_90 %>% select(time,N,J,A,Exposed,I)
ReboundSim_90 = ReboundSim_90 %>% pivot_longer(cols = c(N,J,A,Exposed,I)) 
sim_90_noNaup = ggplot(ReboundSim_90, aes(x=time, y = value + 0.01 , group = name, color = name)) + geom_line() + scale_y_log10()
#sim_90_noNaup
#50%
event_data_50 = Harvest_scheme(harvest_nauplii = FALSE,harvest_prop = 0.5)
ReboundSim_50 = data.frame(ode(y = Initial_conditions, times=1:timespan, parms=parameters, hmax=1,
                               method="lsoda", func=Harvest_ODE,events=list(data=event_data_50)))
ReboundSim_50 = data.frame(ReboundSim_50)
ReboundSim_50[,"Exposed"] = apply(X=ReboundSim_50[,which(str_detect(colnames(ReboundSim_50),"E"))],MARGIN=1,FUN = sum) 
ReboundSim_50 = ReboundSim_50 %>% select(time,N,J,A,Exposed,I)
ReboundSim_50 = ReboundSim_50 %>% pivot_longer(cols = c(N,J,A,Exposed,I)) 
sim_50_noNaup = ggplot(ReboundSim_50, aes(x=time, y = value + 0.01 , group = name, color = name)) + geom_line() + scale_y_log10()
#sim_50_noNaup

#control, 50, 90 for harvest_nauplii_true and lambda = 1.3
plot_grid(control_noNaup,sim_50_noNaup,sim_90_noNaup)  #,labels=c('control', '50%','90%')) + draw_label("Log densities/L when harvest_nauplii = TRUE and lambda =1.3", fontface='bold')
plot_grid(control,sim_50,sim_90,control_noNaup,sim_50_noNaup,sim_90_noNaup)


#combined plot 
sim_90_data <- as.data.frame(sim_90[["data"]]) 
sim_50_data <- as.data.frame(sim_50[["data"]]) 
control_data <- as.data.frame(control[["data"]]) 

control_data$treatment <- "control"
sim_50_data$treatment <- "50% with nauplii harvest"
sim_90_data$treatment <- "90% with nauplii harvest"

combined_df <- bind_rows(control_data, sim_90_data,sim_50_data)

PlotwithNaupliiHarvest = ggplot(combined_df, aes(x = time, group = name, color = name, y = value + 0.01)) +
  geom_line(size=0.75) +
  facet_wrap(~factor(treatment, c("control", "50% with nauplii harvest", "90% with nauplii harvest")),ncol =1) + 
  labs(x = "Time (days)", y = "Density per L", color = "Stage") +
  theme_minimal() +
  theme(
    strip.text = element_text(hjust = 0, face = "bold"),
    strip.background = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  ) + theme_classic() + 
  theme( strip.text = element_text(face = "bold", size = 9),     
    strip.background = element_blank()) + scale_color_manual(values = wes_palette("Darjeeling1", n = 5)) +scale_y_log10() 

#combined plot for no harvest
sim_90_noNaup_data <- as.data.frame(sim_90_noNaup[["data"]]) 
sim_50_noNaup_data <- as.data.frame(sim_50_noNaup[["data"]]) 
control_noNaup_data <- as.data.frame(control_noNaup[["data"]]) 

control_noNaup_data$treatment <- "control"
sim_50_noNaup_data$treatment <- "50% no nauplii harvest"
sim_90_noNaup_data$treatment <- "90% no nauplii harvest"

combined_noNaup_df <- bind_rows(control_noNaup_data, sim_90_noNaup_data,sim_50_noNaup_data)

combined_noNaup_df$name <- factor(combined_noNaup_df$name)

PlotNoNaupliiHarvest = ggplot(combined_noNaup_df, aes(x = time, group = name, color = name, y = value + 0.01)) +
  geom_line(size=0.75) +
  facet_wrap(~factor(treatment, c("control", "50% no nauplii harvest", "90% no nauplii harvest")),ncol =1) + 
  labs(x = "Time (days)", y = "Density per L", color = "Stage") +
  theme_minimal() +
  theme(
    strip.text = element_text(hjust = 0, face = "bold"),
    strip.background = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  ) + theme_classic() + 
  theme( strip.text = element_text(face = "bold", size = 9),     
         strip.background = element_blank()) + scale_color_manual(values = wes_palette("Darjeeling1", n = 5)) + scale_y_log10() 

plot_grid(PlotwithNaupliiHarvest,PlotNoNaupliiHarvest)

combined_noNaup_df$treatment <- gsub(" no nauplii harvest", "", combined_noNaup_df$treatment)

combinedplot = ggplot(combined_noNaup_df, aes(x = time, group = interaction(name,treatment), linetype = treatment, color = name, y = value + 0.01)) +
  geom_line(size=0.75) +
  labs(x = "Time (days)", y = "Density per L", color = "Stage") +
  theme_minimal() +
  theme(
    strip.text = element_text(hjust = 0, face = "bold"),
    strip.background = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  ) + theme_classic() + 
  theme( strip.text = element_text(face = "bold", size = 9),     
         strip.background = element_blank()) + scale_y_log10() + scale_color_manual(values = wes_palette("Darjeeling1", n = 5)) +
  ggtitle("Selective Harvest") + theme(text = element_text(size = 15)) + 
  scale_linetype_manual(
    name = "treatment",
    values = c("50%" = "solid", 
               "90%" = "dotted", 
               "control" = "dashed"),
    labels = c("50%", "90%", "control")
  )
  

combinedplot

combined_df$treatment <- gsub(" with nauplii harvest", "", combined_df$treatment)

combinedplot2 = ggplot(combined_df, aes(x = time, group = interaction(name,treatment), linetype = treatment, color = name, y = value + 0.01)) +
  geom_line(size=0.75) +
  labs(x = "Time (days)", y = "Density per L", color = "Stage") +
  theme_minimal() +
  theme(
    strip.text = element_text(hjust = 0, face = "bold"),
    strip.background = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  ) + theme_classic() + 
  theme( strip.text = element_text(face = "bold", size = 9),     
         strip.background = element_blank()) + scale_y_log10()+ scale_color_manual(values = wes_palette("Darjeeling1", n = 5)) +
  ggtitle("Complete Harvest") + theme(text = element_text(size = 15)) + 
  scale_linetype_manual(
    name = "treatment",
    values = c("50%" = "solid", 
               "90%" = "dotted", 
               "control" = "dashed"),
    labels = c("50%", "90%", "control")
  )


combinedplot2

#plot_grid(combinedplot,combinedplot2)

library(ggpubr)

ggarrange(
  combinedplot2,
  combinedplot,
  labels = "AUTO",         # Adds automatic labels A, B, ...
  common.legend = TRUE,    # Share a common legend
  legend = "right",       # Put the legend at the bottom
  align = "h",             # Align plots horizontally (x-axis alignment)
  nrow = 2                 # Arrange in 2 rows
)







# Filter out 50% treatment from both datasets and combine E and I
library(ggpubr)  # Required for ggarrange()

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
combined_df = combined_df_wide %>% pivot_longer(cols=c("N","J","I","NIA"),names_to = "name",values_to = "value")
  
combined_noNaup_df_wide = combined_noNaup_df %>% pivot_wider(names_from = name, values_from = value)
combined_noNaup_df_wide$NIA <- combined_noNaup_df_wide$A + combined_noNaup_df_wide$Exposed
#remove A and E columns
combined_noNaup_df_wide$A <- NULL
combined_noNaup_df_wide$Exposed <- NULL 
#pivot longer   
combined_noNaup_df = combined_noNaup_df_wide %>% pivot_longer(cols=c("N","J","I","NIA"),names_to = "name",values_to = "value")


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

#THIS IS THE GRAPH!!!
ggplot(combined_df,aes(x=time,y=value + 0.01, color = name,linetype=treatment)) + geom_line(size = 0.75) +
  labs(x = "Time (days)", y = "Density per L", color = "Stage") +
  theme_classic(base_size = 15) +
  scale_y_log10() + scale_color_manual(values = wes_palette("Darjeeling1")) 

combined_df_I = combined_df %>% filter(name == "I")
combined_df_J = combined_df %>% filter(name == "J")
combined_df_N = combined_df %>% filter(name == "N")
combined_df_NIA = combined_df %>% filter(name == "NIA")


library(RColorBrewer)
I = ggplot(combined_df_I,aes(x=time,y=value + 0.01, color=treatment)) + geom_line(size = 1) +
  labs(x = "Time (days)", y = "I per L", color = "Stage") +
  theme_classic(base_size = 15) +scale_color_brewer(palette = "Dark2") 

J = ggplot(combined_df_J,aes(x=time,y=value + 0.01, color=treatment)) + geom_line(size = 1) +
  labs(x = "Time (days)", y = "J per L", color = "Stage") +
  theme_classic(base_size = 15) +scale_color_brewer(palette = "Dark2") 

N = ggplot(combined_df_N,aes(x=time,y=value + 0.01, color=treatment)) + geom_line(size = 1) +
  labs(x = "Time (days)", y = "N per L", color = "Stage") +
  theme_classic(base_size = 15) + scale_color_brewer(palette = "Dark2")  

NIA = ggplot(combined_df_NIA,aes(x=time,y=value + 0.01, color=treatment)) + geom_line(size = 1) +
  labs(x = "Time (days)", y = "NIA per L", color = "Stage") +
  theme_classic(base_size = 15) + scale_color_brewer(palette = "Dark2") 

library(RColorBrewer)

# Extract the Dark2 palette
dark2_cols <- brewer.pal(8, "Dark2")

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


####Publication Plots####

my_colors <- c(
  control   = "yellowgreen",
  selective = "darkslateblue",
  complete  ="cadetblue3"
)

I = ggplot(combined_df_I, aes(x = time, y = value + 0.01, color = treat_group)) + geom_vline(xintercept = c(8, 36, 64), linetype = "dashed", color = "black", linewidth = 1, alpha = 0.8) +
  geom_line(linewidth=1.5) +
  labs(x = "Time (days)", y = "Infectious Adults \nper L", color = "Removal Type") +
  theme_classic(base_size = 20) +
  scale_color_manual(values = my_colors, labels = c("Control", "Size-unbiased","Size-biased")) 

J = ggplot(combined_df_J, aes(x = time, y = value + 0.01, color = treat_group)) + geom_vline(xintercept = c(8, 36, 64), linetype = "dashed", color = "black", linewidth = 1, alpha = 0.8) +
  geom_line(linewidth=1.5) +
  labs(x = "Time (days)", y = "Juveniles \nper L", color = "Removal Type") +
  theme_classic(base_size = 20) +
  scale_color_manual(values = my_colors,labels = c("Control", "Size-unbiased","Size-biased")) + theme(axis.title.x = element_blank())

N = ggplot(combined_df_N, aes(x = time, y = value + 0.01, color = treat_group)) + geom_vline(xintercept = c(8, 36, 64), linetype = "dashed", color = "black", linewidth = 1, alpha = 0.8) +
  geom_line(linewidth=1.5) +
  labs(x = "Time (days)", y = "Nauplii \nper L", color = "Removal Type") +
  theme_classic(base_size = 20) +
  scale_color_manual(values = my_colors,labels = c("Control", "Size-unbiased","Size-biased")) + theme(axis.title.x = element_blank())

NIA = ggplot(combined_df_NIA, aes(x = time, y = value + 0.01, color = treat_group)) + geom_vline(xintercept = c(8, 36, 64), linetype = "dashed", color = "black", linewidth = 1, alpha = 0.8) +
  geom_line(linewidth=1.5) +
  labs(x = "Time (days)", y = "Non-Infected Adults \nper L", color = "Removal Type") +
  theme_classic(base_size = 20) +
  scale_color_manual(values = my_colors,labels = c("Control", "Size-unbiased","Size-biased")) + theme(axis.title.x = element_blank())


library(ggpubr)

ggarrange(
  J,
  N,
  NIA,
  I,        
  common.legend = TRUE,    # Share a common legend
  legend = "right",       # Put the legend at the bottom
  align = "h",             # Align plots horizontally (x-axis alignment)
  nrow = 4                 # Arrange in 2 rows
)

