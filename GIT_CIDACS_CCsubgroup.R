############################################################

# Analysis from:
# 'hort-term ambient heat exposure and low APGAR score in newborns: A time-stratified case-crossover analysis 
# in São Paulo state, Brazil'

# Date of code creation:
# 12/2/2025

# Code author:
# Michelle Del Carretto

# Purpose:
# 1) Conduct subgroup analysis

############################################################
# Load datasets #
############################################################

setwd()
data <- read.csv("GIT_CCdataset.csv")

library(dlnm); library(survival);library(splines)

############################################################
# Dataset preparation #
############################################################

# set percentiles of interest
pop_quantiles <-c(14.8, 20.9, 26.1)

# set average humidity 
data$RH01 <- rowMeans(data[, c("Hmed", "Hmed1")]) 

# convert variables to factors
cols_to_convert <- c("bin_mage", "bin_edu","ethnicity","bin_parity_2","bin_prenatal_care", "sex", "q_IBP", "Koppen" )
data[cols_to_convert] <- lapply(data[cols_to_convert], factor) # Convert the specified columns to factors
str(data)

# stratify into dataframes
df_bel_20 <- subset(data, bin_mage == "<20")
df_20_34 <- subset(data, bin_mage == "20-34")
df_ab_35 <- subset(data, bin_mage == ">=35")

df_furthered <- subset(data, bin_edu == "Further education")
df_highschool <- subset(data, bin_edu == "Highschool or below")

df_asian <- subset(data, ethnicity == "Asian")
df_black <- subset(data, ethnicity == "Black")
df_indig <- subset(data, ethnicity == "Indigenous")
df_mixed <- subset(data, ethnicity == "Mixed")
df_white <- subset(data, ethnicity == "White")

df_nulli <- subset(data, bin_parity_2 == "Nulliparous")
df_primi <- subset(data, bin_parity_2 == "Primiparous")
df_multi <- subset(data, bin_parity_2 == "Multiparous")

df_delay <- subset(data, bin_prenatal_care == "Access after 1st trimester")
df_timely <- subset(data, bin_prenatal_care == "Access in 1st trimester")

df_mal <- subset(data, sex == "Male")
df_fem <- subset(data, sex == "Female")

df_one <- subset(data, q_IBP == "1")
df_two <- subset(data, q_IBP == "2")
df_three <- subset(data, q_IBP == "3")
df_four <- subset(data, q_IBP == "4")

df_tropical <- subset(data, Koppen %in% c("Af", "Aw")) # broad definition
df_temperate <- subset(data, Koppen %in% c("Cfa", "Cfb", "Cwa", "Cwb")) # broad definition
df_af <- subset(data, Koppen == "Af")
df_aw <- subset(data, Koppen == "Aw")
df_cfa <- subset(data, Koppen == "Cfa")
df_cfb <- subset(data, Koppen == "Cfb")
df_cwa <- subset(data, Koppen == "Cwa")
df_cwb <- subset(data, Koppen == "Cwb")

############################################################
# Subgroup analyses #
############################################################

# List of datasets to iterate over
datasets <-   list(df_bel_20, df_20_34, df_ab_35, df_furthered,        df_highschool,         df_asian, df_black, df_indig,    df_mixed, df_white, df_nulli,      df_primi,     df_multi,    df_delay,                     df_timely,       df_mal,  df_fem, df_one,  df_two, df_three, df_four )
dataset_names <- c("<20",     "20-34",  ">=35",   "Further education", "Highschool or below", "Asian",  "Black", "Indigenous", "Mixed",  "White",  "Nulliparous","Primiparous","Multiparous","Access after 1st trimester","Access in 1st trimester", "Male", "Female","IBP 1", "IBP 2", "IBP 3", "IBP 4")

# Initialize an empty list to store ORs
or_results <- list()  

# Loop through the datasets
for (i in 1:length(datasets)) {
  
  # Extract the dataset
  data <- datasets[[i]]
  
  # Create temperature matrix
  Tmed01 <- as.matrix(data[, c("Tmed", "Tmed1")])
  
  # Define crossbasis
  cb <- crossbasis(Tmed01, lag = 1, 
                   argvar = list(fun = "ns", df = 2), 
                   arglag = list(fun = "lin", int = T))  
  
  # Calculate RH over 0-1 day lag
  data$RH01 <- rowMeans(data[, c("Hmed", "Hmed1")]) 
  
  # Run model
  mod <- clogit(bin_apgar_numeric ~ cb + ns(RH01, df = 3) + strata(unique_id), data) 
  
  # summary
  print(paste("Model summary for:", dataset_names[i]))
  print(summary(mod))  
  
  # Predict ORs centered at the 50th percentile
  cp <- crosspred(cb, mod, cen = pop_quantiles[2], by = 0.1)
  
  # Get predictions for ORs
  or_results[[dataset_names[i]]] <- c(
    OR_Lag01_heat = round(cp$allRRfit["26.1"],2),
    CI_Low_Lag01_heat = round(cp$allRRlow["26.1"],2),
    CI_High_Lag01_heat = round(cp$allRRhigh["26.1"],2),
    OR_Lag0_heat = round(cp$matRRfit["26.1", "lag0"],2),
    CI_Low_Lag0_heat = round(cp$matRRlow["26.1", "lag0"],2),
    CI_High_Lag0_heat = round(cp$matRRhigh["26.1", "lag0"],2),
    OR_Lag1_heat = round(cp$matRRfit["26.1", "lag1"],2),
    CI_Low_Lag1_heat = round(cp$matRRlow["26.1", "lag1"],2),
    CI_High_Lag1_heat = round(cp$matRRhigh["26.1", "lag1"],2),
    OR_Lag01_cool = round(cp$allRRfit["14.8"],2),
    CI_Low_Lag01_cool = round(cp$allRRlow["14.8"],2),
    CI_High_Lag01_cool = round(cp$allRRhigh["14.8"],2),
    OR_Lag0_cool = round(cp$matRRfit["14.8", "lag0"],2),
    CI_Low_Lag0_cool = round(cp$matRRlow["14.8", "lag0"],2),
    CI_High_Lag0_cool = round(cp$matRRhigh["14.8", "lag0"],2),
    OR_Lag1_cool = round(cp$matRRfit["14.8", "lag1"],2),
    CI_Low_Lag1_cool = round(cp$matRRlow["14.8", "lag1"],2),
    CI_High_Lag1_cool = round(cp$matRRhigh["14.8", "lag1"],2)
  )
}

# Convert the results list to a data frame 
or_results_df <- as.data.frame(do.call(cbind, or_results))

# Print the OR results data frame
print(or_results_df)

############################################################
# Subgroup analyses (Koppen) #
############################################################

# set percentiles of interest
pop_quantiles_tropical <-c(17.4, 23.6,  28)
pop_quantiles_temperate <-c(14.6, 20.6, 25.6)

# extract temperature matrix for each
Tmed_tropical <- as.matrix(df_tropical[, c("Tmed", "Tmed1")])
Tmed_temperate <- as.matrix(df_temperate[, c("Tmed", "Tmed1")])
Tmed_af <- as.matrix(df_af[, c("Tmed", "Tmed1")])
Tmed_aw <- as.matrix(df_aw[, c("Tmed", "Tmed1")])
Tmed_cfa <- as.matrix(df_cfa[, c("Tmed", "Tmed1")])
Tmed_cfb <- as.matrix(df_cfb[, c("Tmed", "Tmed1")])
Tmed_cwa <- as.matrix(df_cwa[, c("Tmed", "Tmed1")])
Tmed_cwb <- as.matrix(df_cwb[, c("Tmed", "Tmed1")])

# define crossbases
cb_tropical <- crossbasis(Tmed_tropical,lag=1,argvar=list(fun = "ns", df = 2) ,arglag=list(fun = "lin",int=T)) 
cb_temperate <- crossbasis(Tmed_temperate,lag=1,argvar=list(fun = "ns", df = 2) ,arglag=list(fun = "lin",int=T)) 
cb_af <- crossbasis(Tmed_af,lag=1,argvar=list(fun = "ns", df = 2) ,arglag=list(fun = "lin",int=T)) 
cb_aw <- crossbasis(Tmed_aw,lag=1,argvar=list(fun = "ns", df = 2) ,arglag=list(fun = "lin",int=T))
cb_cfa <- crossbasis(Tmed_cfa,lag=1,argvar=list(fun = "ns", df = 2) ,arglag=list(fun = "lin",int=T)) 
cb_cfb <- crossbasis(Tmed_cfb,lag=1,argvar=list(fun = "ns", df = 2) ,arglag=list(fun = "lin",int=T))
cb_cwa <- crossbasis(Tmed_cwa,lag=1,argvar=list(fun = "ns", df = 2) ,arglag=list(fun = "lin",int=T)) 
cb_cwb <- crossbasis(Tmed_cwb,lag=1,argvar=list(fun = "ns", df = 2) ,arglag=list(fun = "lin",int=T))

# run models
mod_tropical <- clogit(bin_apgar_numeric ~ cb_tropical + ns(RH01, df=3) + strata(unique_id), df_tropical)
mod_temperate <- clogit(bin_apgar_numeric ~ cb_temperate + ns(RH01, df=3) + strata(unique_id), df_temperate)
mod_af <- clogit(bin_apgar_numeric ~ cb_af + ns(RH01, df=3) + strata(unique_id), df_af)
mod_aw <- clogit(bin_apgar_numeric ~ cb_aw + ns(RH01, df=3) + strata(unique_id), df_aw)
mod_cfa <- clogit(bin_apgar_numeric ~ cb_cfa + ns(RH01, df=3) + strata(unique_id), df_cfa)
mod_cfb <- clogit(bin_apgar_numeric ~ cb_cfb + ns(RH01, df=3) + strata(unique_id), df_cfb)
mod_cwa <- clogit(bin_apgar_numeric ~ cb_cwa + ns(RH01, df=3) + strata(unique_id), df_cwa)
mod_cwb <- clogit(bin_apgar_numeric ~ cb_cwb + ns(RH01, df=3) + strata(unique_id), df_cwb)

# get summaries
summary(mod_tropical)
summary(mod_temperate)

# model predictions
cp_tropical  <- crosspred(cb_tropical, mod_tropical, cen = pop_quantiles_tropical[2] , by=0.1)
cp_temperate <- crosspred(cb_temperate, mod_temperate, cen = pop_quantiles_temperate[2] , by=0.1)
cp_af  <- crosspred(cb_af, mod_af, cen = pop_quantiles_tropical[2] , by=0.1)
cp_aw <- crosspred(cb_aw, mod_aw, cen = pop_quantiles_tropical[2] , by=0.1)
cp_cfa  <- crosspred(cb_cfa, mod_cfa, cen = pop_quantiles_temperate[2] , by=0.1)
cp_cfb <- crosspred(cb_cfb, mod_cfb, cen = pop_quantiles_temperate[2] , by=0.1)
cp_cwa  <- crosspred(cb_cwa, mod_cwa, cen = pop_quantiles_temperate[2] , by=0.1)
cp_cwb <- crosspred(cb_cwb, mod_cwb, cen = pop_quantiles_temperate[2] , by=0.1)

# extract predictions to dataframe
table_koppen <- data.frame(
  Object = c("A-Tropcial","C-Temperate","Af","Aw","Cfa","Cfb","Cwa","Cwb"),
  OR_Lag01_heat = c(cp_tropical$allRRfit["28"],cp_temperate$allRRfit["25.6"] , cp_af$allRRfit["28"], cp_aw$allRRfit["28"], cp_cfa$allRRfit["25.6"], cp_cfb$allRRfit["25.6"], cp_cwa$allRRfit["25.6"], cp_cwb$allRRfit["25.6"]),
  CI_Low_Lag01_heat = c(cp_tropical$allRRlow["28"],cp_temperate$allRRlow["25.6"] , cp_af$allRRlow["28"], cp_aw$allRRlow["28"], cp_cfa$allRRlow["25.6"], cp_cfb$allRRlow["25.6"], cp_cwa$allRRlow["25.6"], cp_cwb$allRRlow["25.6"]),
  CI_High_Lag01_heat = c(cp_tropical$allRRhigh["28"], cp_temperate$allRRhigh["25.6"], cp_af$allRRhigh["28"], cp_aw$allRRhigh["28"], cp_cfa$allRRhigh["25.6"], cp_cfb$allRRhigh["25.6"], cp_cwa$allRRhigh["25.6"], cp_cwb$allRRhigh["25.6"]),
  OR_Lag0_heat = c(cp_tropical$matRRfit["28", "lag0"], cp_temperate$matRRfit["25.6", "lag0"],cp_af$matRRfit["28", "lag0"], cp_aw$matRRfit["28", "lag0"], cp_cfa$matRRfit["25.6", "lag0"], cp_cfb$matRRfit["25.6", "lag0"], cp_cwa$matRRfit["25.6", "lag0"], cp_cwb$matRRfit["25.6", "lag0"]),
  CI_Low_Lag0_heat = c(cp_tropical$matRRlow["28", "lag0"], cp_temperate$matRRlow["25.6", "lag0"], cp_af$matRRlow["28", "lag0"], cp_aw$matRRlow["28", "lag0"], cp_cfa$matRRlow["25.6", "lag0"], cp_cfb$matRRlow["25.6", "lag0"], cp_cwa$matRRlow["25.6", "lag0"], cp_cwb$matRRlow["25.6", "lag0"]),
  CI_High_Lag0_heat = c(cp_tropical$matRRhigh["28", "lag0"], cp_temperate$matRRhigh["25.6", "lag0"],cp_af$matRRhigh["28", "lag0"], cp_aw$matRRhigh["28", "lag0"], cp_cfa$matRRhigh["25.6", "lag0"], cp_cfb$matRRhigh["25.6", "lag0"], cp_cwa$matRRhigh["25.6", "lag0"], cp_cwb$matRRhigh["25.6", "lag0"]),
  OR_Lag1_heat = c(cp_tropical$matRRfit["28", "lag1"], cp_temperate$matRRfit["25.6", "lag1"], cp_af$matRRfit["28", "lag1"], cp_aw$matRRfit["28", "lag1"], cp_cfa$matRRfit["25.6", "lag1"], cp_cfb$matRRfit["25.6", "lag1"], cp_cwa$matRRfit["25.6", "lag1"], cp_cwb$matRRfit["25.6", "lag1"]),
  CI_Low_Lag1_heat = c(cp_tropical$matRRlow["28", "lag1"], cp_temperate$matRRlow["25.6", "lag1"], cp_af$matRRlow["28", "lag1"], cp_aw$matRRlow["28", "lag1"], cp_cfa$matRRlow["25.6", "lag1"], cp_cfb$matRRlow["25.6", "lag1"], cp_cwa$matRRlow["25.6", "lag1"], cp_cwb$matRRlow["25.6", "lag1"]),
  CI_High_Lag1_heat = c(cp_tropical$matRRhigh["28", "lag1"], cp_temperate$matRRhigh["25.6", "lag1"],cp_af$matRRhigh["28", "lag1"], cp_aw$matRRhigh["28", "lag1"], cp_cfa$matRRhigh["25.6", "lag1"], cp_cfb$matRRhigh["25.6", "lag1"], cp_cwa$matRRhigh["25.6", "lag1"], cp_cwb$matRRhigh["25.6", "lag1"]),
  OR_Lag01_cool = c(cp_tropical$allRRfit["17.4"],cp_temperate$allRRfit["14.6"] , cp_af$allRRfit["17.4"], cp_aw$allRRfit["17.4"], cp_cfa$allRRfit["14.6"], cp_cfb$allRRfit["14.6"], cp_cwa$allRRfit["14.6"], cp_cwb$allRRfit["14.6"]),
  CI_Low_Lag01_cool = c(cp_tropical$allRRlow["17.4"],cp_temperate$allRRlow["14.6"] , cp_af$allRRlow["17.4"], cp_aw$allRRlow["17.4"], cp_cfa$allRRlow["14.6"], cp_cfb$allRRlow["14.6"], cp_cwa$allRRlow["14.6"], cp_cwb$allRRlow["14.6"]),
  CI_High_Lag01_cool = c(cp_tropical$allRRhigh["17.4"], cp_temperate$allRRhigh["14.6"], cp_af$allRRhigh["17.4"], cp_aw$allRRhigh["17.4"], cp_cfa$allRRhigh["14.6"], cp_cfb$allRRhigh["14.6"], cp_cwa$allRRhigh["14.6"], cp_cwb$allRRhigh["14.6"]),
  OR_Lag0_cool = c(cp_tropical$matRRfit["17.4", "lag0"], cp_temperate$matRRfit["14.6", "lag0"],cp_af$matRRfit["17.4", "lag0"], cp_aw$matRRfit["17.4", "lag0"], cp_cfa$matRRfit["14.6", "lag0"], cp_cfb$matRRfit["14.6", "lag0"], cp_cwa$matRRfit["14.6", "lag0"], cp_cwb$matRRfit["14.6", "lag0"]),
  CI_Low_Lag0_cool = c(cp_tropical$matRRlow["17.4", "lag0"], cp_temperate$matRRlow["14.6", "lag0"], cp_af$matRRlow["17.4", "lag0"], cp_aw$matRRlow["17.4", "lag0"], cp_cfa$matRRlow["14.6", "lag0"], cp_cfb$matRRlow["14.6", "lag0"], cp_cwa$matRRlow["14.6", "lag0"], cp_cwb$matRRlow["14.6", "lag0"]),
  CI_High_Lag0_cool = c(cp_tropical$matRRhigh["17.4", "lag0"], cp_temperate$matRRhigh["14.6", "lag0"],cp_af$matRRhigh["17.4", "lag0"], cp_aw$matRRhigh["17.4", "lag0"], cp_cfa$matRRhigh["14.6", "lag0"], cp_cfb$matRRhigh["14.6", "lag0"], cp_cwa$matRRhigh["14.6", "lag0"], cp_cwb$matRRhigh["14.6", "lag0"]),
  OR_Lag1_cool = c(cp_tropical$matRRfit["17.4", "lag1"], cp_temperate$matRRfit["14.6", "lag1"], cp_af$matRRfit["17.4", "lag1"], cp_aw$matRRfit["17.4", "lag1"], cp_cfa$matRRfit["14.6", "lag1"], cp_cfb$matRRfit["14.6", "lag1"], cp_cwa$matRRfit["14.6", "lag1"], cp_cwb$matRRfit["14.6", "lag1"]),
  CI_Low_Lag1_cool = c(cp_tropical$matRRlow["17.4", "lag1"], cp_temperate$matRRlow["14.6", "lag1"], cp_af$matRRlow["17.4", "lag1"], cp_aw$matRRlow["17.4", "lag1"], cp_cfa$matRRlow["14.6", "lag1"], cp_cfb$matRRlow["14.6", "lag1"], cp_cwa$matRRlow["14.6", "lag1"], cp_cwb$matRRlow["14.6", "lag1"]),
  CI_High_Lag1_cool = c(cp_tropical$matRRhigh["17.4", "lag1"], cp_temperate$matRRhigh["14.6", "lag1"],cp_af$matRRhigh["17.4", "lag1"], cp_aw$matRRhigh["17.4", "lag1"], cp_cfa$matRRhigh["14.6", "lag1"], cp_cfb$matRRhigh["14.6", "lag1"], cp_cwa$matRRhigh["14.6", "lag1"], cp_cwb$matRRhigh["14.6", "lag1"]))
  
# Round all numeric values to 2 decimal places
table_koppen[,-1] <- round(table_koppen[,-1], 2)

# Transpose the data frame
table_koppen <- t(table_koppen)

# Convert it back to a data frame for better readability
table_koppen <- as.data.frame(table_koppen)

# Set column names for clarity (optional, if you want better column names)
colnames(table_koppen) <- table_koppen$Object

# View the transposed and rounded table
print(table_koppen)

