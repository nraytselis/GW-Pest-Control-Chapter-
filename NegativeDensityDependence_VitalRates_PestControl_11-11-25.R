library(coda)
library(adaptMCMC)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(ggpubr)


setwd("~/Desktop/Rscripts/Data")


####Summary####
#This code uses the 32 mcmc chains to determine which vital rates are under the strongest negative density dependence in the pest control mesocosm experiment.


###processing full chains
  fullA = readRDS("Rebound_parameters_full_disp250kff_5.RDA")
  fullB = readRDS("Rebound_parameters_full_disp250kff_3.RDA")
  fullC = readRDS("Rebound_parameters_full_disp250kff_05.RDA")
  fullD = readRDS("Rebound_parameters_full_disp250kff_01.RDA")
  full = c(fullA, fullB, fullC, fullD)
  
  
  #extract samples from all chains in the list 
  get_samples <- function(chain.list) {
    L <- length(chain.list)
    samples_list <- vector("list", L)  # empty list to store samples
    
    for (i in 1:L) {
      samples <- chain.list[[i]]$samples  # extract 'samples' from each chain
      samples_list[[i]] <- coda::mcmc(samples)  # convert to mcmc object
    }
    
    return(samples_list)
  }
  
  all_chains <- get_samples(full)
  
  
  #Each stage is multiplied by the proportion of the population within that stage (averaged across the entire experiment)
  
  #births
  Per_capita_births = function(b_M, comp_b, c_N, c_J,Total_Abundance, VOL=15){ #these are the parameters that contribute to the birth rate 
    N = 10/18*Total_Abundance
    J = 7/18*Total_Abundance
    A = 1/18*Total_Abundance
    b_M*exp(-comp_b/VOL*(c_N*N + c_J*J + A))
  } 
  
  #maturation 
  Per_capita_maturation = function(m_N,comp_m,c_N,c_J,Total_Abundance, VOL=15){
    N = 10/18*Total_Abundance
    J = 7/18*Total_Abundance
    A = 1/18*Total_Abundance
    m_N*exp(-comp_m/VOL*(c_N*N + c_J*J + A)) 
  }
  
  
  #cannibalism
  Per_capita_cann = function(cann,I = 0, A, N = 0,Total_Abundance, VOL=15){
    A = 1/18*Total_Abundance
    cann*A
  }
  
  #deaths
  Per_capita_deaths = function(d_N,comp_d,c_N,c_J,Total_Abundance,cann,VOL=15){
    N = 10/18*Total_Abundance
    J = 7/18*Total_Abundance
    A = 1/18*Total_Abundance
    d_N*exp(comp_d / VOL * (c_N * N + c_J * J + A)) + cann*A
  }
  
  
  
  #thin all chains
  all_chains_coda <- convert.to.coda(all_chains)
  all_chains_coda_thin <- window(all_chains_coda,start=60000,end=100000,thin=1000)
  
  #make a combined list with all parameter sets from each chain
  pars <- as.mcmc.list(all_chains_coda_thin)
  
  #loop through all parameter sets
  pars_matrix <- as.matrix(pars)  
  n_iter <- nrow(pars_matrix)
  
  
  #adult range up to 5,000/tank
  #juvenile range up to 35,000/tank
  #nauplii range up to 50,000/tank
  #total abundance up to 90k/tank
  
  #set up ranges
  Total_Abundance_range <- 10 * (0:9000) #calculate from raw data 
  n_Total_Abundance<- length(Total_Abundance_range)
  
  #adult_range <- 10 * (0:500)  # 501 values
  #n_adults <- length(adult_range)
  
  #results matrices
  matrix_Total_Abundance_births <- matrix(NA, nrow = n_iter, ncol = n_Total_Abundance)
  matrix_Total_Abundance_maturation <- matrix(NA, nrow = n_iter, ncol = n_Total_Abundance)
  matrix_Total_Abundance_deaths <- matrix(NA, nrow = n_iter, ncol = n_Total_Abundance)
  
  #matrix cann
  matrix_adults_cann <- matrix(NA, nrow = n_iter, ncol = n_Total_Abundance)
  
  
  #births
  for (i in 1:n_iter) {
    p <- pars_matrix[i, ]
    
    matrix_Total_Abundance_births[i, ] <- Per_capita_births(
      b_M = p["b_M"],
      comp_b = p["comp_b"],
      c_N = p["c_N"],
      c_J = p["c_J"],
      Total_Abundance = Total_Abundance_range
    )
  }
  
  Hi = apply(matrix_Total_Abundance_births, 2, quantile, 0.975)
  Per_capita_birth_rate = apply(matrix_Total_Abundance_births, 2, quantile, probs=0.5)
  Lo = apply(matrix_Total_Abundance_births, 2, quantile, probs=0.025)
  plotting_data_TotalAbundance = data.frame("Total_Abundance" = Total_Abundance_range, Hi, Per_capita_birth_rate, Lo)
  births_Total = ggplot(data=plotting_data_TotalAbundance, aes(x=Total_Abundance, y=Per_capita_birth_rate)) + 
    geom_line() + geom_ribbon(aes(ymin=Lo, ymax=Hi, alpha=0.2),fill = "#476F84FF") + ylim(0,6.5) + labs(x="Total Abundance", y="Per Capita Birth Rate") + theme_classic(base_size = 13) + theme(legend.position="none")
    
  births_Total
  
  max = max(plotting_data_TotalAbundance$Per_capita_birth_rate)
  min = min(plotting_data_TotalAbundance$Per_capita_birth_rate)
  
  max/min
foldchangeB = (min-max)/min
  
  #maturation
  for (i in 1:n_iter) {
    p <- pars_matrix[i, ]
    
    matrix_Total_Abundance_maturation[i, ] <- Per_capita_maturation(
      m_N = p["m_N"],
      comp_m = p["comp_m"],
      c_N = p["c_N"],
      c_J = p["c_J"],
      Total_Abundance = Total_Abundance_range
    )
  }
  
  Hi = apply(matrix_Total_Abundance_maturation, 2, quantile, 0.975)
  Per_capita_maturation_rate = apply(matrix_Total_Abundance_maturation, 2, quantile, probs=0.5)
  Lo = apply(matrix_Total_Abundance_maturation, 2, quantile, probs=0.025)
  plotting_data_TotalAbundance = data.frame("Total_Abundance" = Total_Abundance_range, Hi, Per_capita_maturation_rate, Lo)
  maturation_Total = ggplot(data=plotting_data_TotalAbundance, aes(x=Total_Abundance, y=Per_capita_maturation_rate)) + 
    geom_line() + geom_ribbon(aes(ymin=Lo, ymax=Hi, alpha=0.2),fill ="olivedrab") + ylim(0,0.6) + labs(x="Total Abundance", y="Per Capita Maturation Rate") + theme_classic(base_size = 13) + theme(legend.position="none") 
  maturation_Total
  
  max = max(plotting_data_TotalAbundance$Per_capita_maturation_rate)
  min = min(plotting_data_TotalAbundance$Per_capita_maturation_rate)
  
foldchangeM = (min-max)/min
  
  #deaths
  for (i in 1:n_iter) {
    p <- pars_matrix[i, ]
    
    matrix_Total_Abundance_deaths[i, ] <- Per_capita_deaths(
      d_N    = p["d_N"],
      comp_d = p["comp_d"],
      c_N    = p["c_N"],
      c_J    = p["c_J"],
      cann   = p["cann"],
      Total_Abundance = Total_Abundance_range
    )
  }
  
  Hi = apply(matrix_Total_Abundance_deaths, 2, quantile, 0.975)
  Per_capita_death_rate = apply(matrix_Total_Abundance_deaths, 2, quantile, probs=0.5)
  Lo = apply(matrix_Total_Abundance_deaths, 2, quantile, probs=0.025)
  plotting_data_TotalAbundance = data.frame("Total_Abundance" = Total_Abundance_range, Hi, Per_capita_death_rate, Lo)
  deaths_Total = ggplot(data=plotting_data_TotalAbundance, aes(x=Total_Abundance, y=Per_capita_death_rate)) + 
    geom_line() + geom_ribbon(aes(ymin=Lo, ymax=Hi, alpha=0.2),fill = "gold2") + ylim(0,0.15) + labs(x="Total Abundance", y="Per Capita Death Rate") + theme_classic(base_size = 13) + theme(legend.position="none") 

  
  deaths_Total
  
  max = max(plotting_data_TotalAbundance$Per_capita_death_rate)
  min = min(plotting_data_TotalAbundance$Per_capita_death_rate)
  
 foldchangeD = (max-min)/min
  
  
  #cann
  for (i in 1:n_iter) {
    p <- pars_matrix[i, ]
    
    matrix_adults_cann[i, ] <- Per_capita_cann(
      cann = p["cann"],
      Total_Abundance = Total_Abundance_range
    )
  }
  
  Hi = apply(matrix_adults_cann, 2, quantile, 0.975)
  Median = apply(matrix_adults_cann, 2, quantile, probs=0.5)
  Lo = apply(matrix_adults_cann, 2, quantile, probs=0.025)
  plotting_data_A = data.frame("Total_Abundance" = Total_Abundance_range, Hi, Median, Lo)
  cann = ggplot(data=plotting_data_A, aes(x=Total_Abundance, y=Median)) + geom_line() + geom_ribbon(aes(ymin=Lo, ymax=Hi, alpha=0.2),fill = "lightblue") + 
    ylim(0,0.0001) + labs(x = "Total Adults", y = "Adult Cannibalism Rate on Nauplii") + theme_classic(base_size = 13) 
  cann
  
  
  
  max = max(plotting_data_TotalAbundance$Per_capita_death_rate)
  min = min(plotting_data_TotalAbundance$Per_capita_death_rate)
foldchangecann =(max-min/min)

  ggarrange(births_Total,maturation_Total,deaths_Total,cann,
            nrow = 2,
            ncol =2,
            labels = c("A","B","C","D"), 
            legend = "none"
  )


 foldchange = data.frame(foldchangeB,foldchangeM,foldchangeD)
 

 foldchange = pivot_longer(foldchange,cols= 1:3, names_to = "VitalRate",values_to = "Value")

 foldchange$Value = abs( foldchange$Value)
 foldchange$VitalRate <- c("Birth","Maturation","Death")
 
 #plot the fold changes 
 ggplot(data=foldchange,aes(x=VitalRate,y=Value)) + geom_point(size = 4, color = "#476F84FF") + theme_classic() + labs(x = "Vital Rate", y = "Fold Change") + theme_classic(base_size = 15)

 

 #plot vital rates together

 HiB = apply(matrix_Total_Abundance_births, 2, quantile, 0.975)
 Per_capita_birth_rateB = apply(matrix_Total_Abundance_births, 2, quantile, probs=0.5)
 LoB = apply(matrix_Total_Abundance_births, 2, quantile, probs=0.025)
 plotting_data_TotalAbundanceB = data.frame("Total_Abundance" = Total_Abundance_range, HiB, Per_capita_birth_rateB, LoB)
 
 HiM = apply(matrix_Total_Abundance_maturation, 2, quantile, 0.975)
 Per_capita_maturation_rateM = apply(matrix_Total_Abundance_maturation, 2, quantile, probs=0.5)
 LoM = apply(matrix_Total_Abundance_maturation, 2, quantile, probs=0.025)
 plotting_data_TotalAbundanceM = data.frame("Total_Abundance" = Total_Abundance_range, HiM, Per_capita_maturation_rate, LoM)
 
 HiD = apply(matrix_Total_Abundance_deaths, 2, quantile, 0.975)
 Per_capita_death_rateD = apply(matrix_Total_Abundance_deaths, 2, quantile, probs=0.5)
 LoD = apply(matrix_Total_Abundance_deaths, 2, quantile, probs=0.025)
 plotting_data_TotalAbundanceD = data.frame("Total_Abundance" = Total_Abundance_range, HiD, Per_capita_death_rate, LoD)
 
 #cann
 HiCann = apply(matrix_adults_cann, 2, quantile, 0.975)
 Per_capita_death_rateCann = apply(matrix_adults_cann, 2, quantile, probs=0.5)
 LoCann = apply(matrix_adults_cann, 2, quantile, probs=0.025)
 plotting_data_TotalAbundanceCann = data.frame("Total_Abundance" = Total_Abundance_range, HiCann, Per_capita_death_rateCann, LoCann)
 
 
 HiB <- as.numeric(HiB)
 Per_capita_birth_rateB <- as.numeric(Per_capita_birth_rateB)
 LoB <- as.numeric(LoB)
 
 HiM <- as.numeric(HiM)
 Per_capita_maturation_rateM <- as.numeric(Per_capita_maturation_rateM)
 LoM <- as.numeric(LoM)
 
 HiD <- as.numeric(HiD)
 Per_capita_death_rateD <- as.numeric(Per_capita_death_rateD)
 LoD <- as.numeric(LoD)
 
 HiCann <- as.numeric(HiCann)
 Per_capita_death_rateCann <- as.numeric(Per_capita_death_rateCann)
 LoCann <- as.numeric(LoCann)
 
 
 # Births
 plotting_data_TotalAbundanceB <- data.frame(
   Total_Abundance = Total_Abundance_range,
   Hi = HiB,
   Median = Per_capita_birth_rateB,
   Lo = LoB,
   VitalRateType = "Birth"
 )
 
 plotting_data_TotalAbundanceM <- data.frame(
   Total_Abundance = Total_Abundance_range,
   Hi = HiM,
   Median = Per_capita_maturation_rateM,
   Lo = LoM,
   VitalRateType = "Maturation"
 )
 
 plotting_data_TotalAbundanceD <- data.frame(
   Total_Abundance = Total_Abundance_range,
   Hi = HiD,
   Median = Per_capita_death_rateD,
   Lo = LoD,
   VitalRateType = "Death"
 )
 
 plotting_data_TotalAbundanceCann <- data.frame(
   Total_Abundance = Total_Abundance_range,
   Hi = HiCann,
   Median = Per_capita_death_rateCann,
   Lo = LoCann,
   VitalRateType = "Cann"
 )
 
 
 plotting_data_TotalAbundanceCann = plotting_data_TotalAbundanceCann %>% select(c(1:4))
 
 
 vitalrates_long <- bind_rows(
   plotting_data_TotalAbundanceB,
   plotting_data_TotalAbundanceM,
   plotting_data_TotalAbundanceD
 )
 
palette = c(Birth = "#EEDD88", Maturation = "#77AADD",Death = "#AAAA00")
 
 ggplot(vitalrates_long,
        aes(x = Total_Abundance/15, y = Median,
            color = VitalRateType, fill = VitalRateType)) +
   geom_line(size = 1) +
   geom_ribbon(aes(ymin = Lo, ymax = Hi), alpha = 0.2, color = NA) +
   labs(x = expression('Density, L' ^ -1), y = "Per Capita Vital Rate",) +
   theme_classic(base_size = 13) + scale_color_manual(values = palette,name = "Vital Rate") +
   scale_fill_manual(values = palette,name = "Vital Rate") +
   theme(legend.position = c(0.2, 0.2)) + scale_y_log10()
 
 

 #bring in time series data to add a rug
 
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
 real_df$total <- real_df$total/15 #per L
 


colors <- c(
  Birth   = "yellowgreen",
  Death = "darkslateblue",
  Maturation  ="cadetblue3"
)
 
 ggplot(vitalrates_long,
        aes(x = Total_Abundance/15, y = Median,
            color = VitalRateType, fill = VitalRateType)) +
   geom_line(size = 1.5) +
   geom_ribbon(aes(ymin = Lo, ymax = Hi), alpha = 0.2, color = NA) +
   labs(x = expression('Density, L' ^ -1), y = "Per Capita Vital Rate",) +
   theme_classic(base_size = 30) + scale_color_manual(values = colors,name = "Vital Rate") +
   scale_fill_manual(values = colors,name = "Vital Rate") +
   theme(legend.position = c(0.2, 0.2)) + scale_y_log10() + geom_rug(
     data = dplyr::filter(real_df, total <= 6000),
     aes(x = total),
     inherit.aes = FALSE,
     sides = "b"
   )
 
 

 
 

 
 
 
 
 

 