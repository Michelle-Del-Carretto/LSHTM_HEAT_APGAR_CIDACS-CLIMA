############################################################

# Analysis from:
# 'Acute heat exposure and low APGAR score in newborns: A time-stratified case-crossover analysis 
# in São Paulo state, Brazil'

# Date of code creation:
# 11/2/2025

# Code author:
# Michelle Del Carretto

# Purpose:
# 1) Conduct main analysis
# 2) Conduct analysis on categories of APGAR
# 3) Sensitivity analysis

############################################################
# Load datasets #
############################################################

setwd("~/Documents/LSHTM project/R code")
data <- read.csv("GIT_CCdataset.csv")

library(dlnm); library(survival);library(splines)

############################################################
# Main analysis (step-by-step APGAR 0-7) #
############################################################

# create temperature matrix
Tmed01 <- as.matrix(data[, c("Tmed", "Tmed1")])

# define crossbasis
cb <- crossbasis(Tmed01,lag=1,  # feed temperature matrix into crossbasis function (this will model the 2 dimensional relationship between tempeature and time/lags)
                 argvar=list(fun = "ns", df = 2), # how relationship between APGAR and temperature is modelled - flexible spline term with 2degrees of freedom or 1 knot
                 arglag=list(fun = "lin",int=T))  # how relationship between APGAR and lags is modelled - we are restricting the relationship to linear

# calculate RH over 0-1 day lag
data$RH01 <- rowMeans(data[, c("Hmed", "Hmed1")]) 

# run model
mod <- clogit(bin_apgar_numeric ~ cb + ns(RH01, df=3) + strata(unique_id), data) # mean humidity is modelled as a spline with 3 degrees of freedom

# summary
summary(mod) 

# set prediction values (using population-weighted temperature quintiles)
pop_quantiles <-c(14.8, 20.9, 26.1)

# predict ORs centered at the 5th percentile
cp <- crosspred(cb, mod, cen = pop_quantiles[2] , by=0.1)

# get predictions (heat)

cp$allRRfit["26.1"] # cumulative lag 0-1 OR
cp$allRRlow["26.1"] # lower CI
cp$allRRhigh["26.1"] # upper CI

cp$matRRfit["26.1", "lag0"] # lag 0 OR
cp$matRRlow["26.1", "lag0"] 
cp$matRRhigh["26.1", "lag0"] 

cp$matRRfit["26.1", "lag1"] # lag 1 OR
cp$matRRlow["26.1", "lag1"]
cp$matRRhigh["26.1", "lag1"] 

# plot

# overall (cumulative effect over 0-1 day lag) association at different temperatures
plot(cp, "overall", ylab="OR",xlab=expression(paste("Daily mean temperature (",degree,"C)")),cex.axis=1.5, cex.lab=1.75, lwd=2.5,ci.arg=list(density=100,col=grey(0.7)))
abline(v = c(14.8, 20.9, 26.1), col = c("blue", "black", "red"), lty = "dashed", lwd = 2)

# association at (95th percentile) vs (median) on each lag day
plot(cp , var = 26.1, ylab = "OR", xlab = "Lag (days)",cex.axis=1.5, cex.lab=1.75, lwd=2.5, ci.arg = list(density = 100, col = grey(0.7)),xaxt='n')
axis(1, at=c(0:1), cex.axis=1.5)

# association at different temperatures on lag 0
plot(cp, "slices", lag=0, ylab="OR",xlab=expression(paste("Daily mean temperature (",degree,"C)")),cex.axis=1.5, cex.lab=1.75, lwd=2.5, ci.arg=list(density=100,col=grey(0.7)))
abline(v = c(14.8, 20.9, 26.1), col = c("blue", "black", "red"), lty = "dashed", lwd = 2)

# association at different temperatures on lag 1
plot(cp, "slices", lag=1, ylab="OR", xlab=expression(paste("Daily mean temperature (",degree,"C)")),cex.axis=1.5, cex.lab=1.75, lwd=2.5, ci.arg=list(density=100,col=grey(0.7)))
abline(v = c(14.8, 20.9, 26.1), col = c("blue", "black", "red"), lty = "dashed", lwd = 2)

############################################################
# Categories of APGAR (loop) #
############################################################

# filter datasets for categories
data_0_2 <- subset(data, clean_apgar == '0-2')
data_3_5 <- subset(data, clean_apgar == '3-5')
data_6_7 <- subset(data, clean_apgar == '6-7')

# check using table
table(data_0_2$bin_apgar)
table(data_3_5$bin_apgar)
table(data_6_7$bin_apgar)
table(data$bin_apgar)

# List of datasets to iterate over
datasets <- list(data_0_2, data_3_5, data_6_7, data)
dataset_names <- c("APGAR 0-2", "APGAR 3-5", "APGAR 6-7", "APGAR 0-7")

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
  
  # Plotting (only overall plot)
  if (i == 1) {
    par(mfrow = c(2, 2))  
  }
  
  # Overall (cumulative effect over 0-1 day lag)
  plot(cp, "overall", 
       ylab = "OR", 
       xlab = expression(paste("Daily mean temperature (", degree, "C)")), 
       cex.axis = 1.5, cex.lab = 1.5, lwd = 2.5, 
       ci.arg = list(density = 100, col = grey(0.7)))
  abline(v = c(14.8, 20.9, 26.1), col = c("blue", "black", "red"), lty = "dashed", lwd = 2)
  title(main = paste("Overall Plot -", dataset_names[i]))
}


# Convert the results list to a data frame 
or_results_df <- as.data.frame(do.call(cbind, or_results))

# Print the OR results data frame
print(or_results_df)

############################################################
# Sensitivity analysis #
############################################################

### 0 lag model ###

# crossbasis
cb0 <- onebasis(data$Tmed, df = 2) # one dimensional as we are only examining day of delivery here

# model
model0 <- clogit(bin_apgar_numeric ~ cb0 + ns(Hmed, df=3) + strata(unique_id), data)

# summary
summary(model0)

# predictions
cp0 <- crosspred(cb0, model0, cen = pop_quantiles[2], by = 0.1)


### 01 lag models (temperature matrix and humidity average pre-created) ###

# crossbasis
cb2 <- crossbasis(Tmed01,lag=1,
                  argvar=list(fun = "ns", df = 3) , # spline with 2 knots (3 df) in the temperature dimension
                  arglag=list(fun = "lin",int=T)) # linear lag dimension

# models
mod01a <- clogit(bin_apgar_numeric ~ cb2 + ns(RH01, df=3) + strata(unique_id), data)
mod01b <- clogit(bin_apgar_numeric ~ cb + RH01 + strata(unique_id), data) # humidity is linear

# summary
summary(mod01a)
summary(mod01b)

# predictions
cp01a <- crosspred(cb2, mod01a, cen = pop_quantiles[2] , by=0.1)
cp01b <- crosspred(cb, mod01b, cen = pop_quantiles[2] , by=0.1)


### 06 lag model ###

# temp matrix 
Tmed06 <- as.matrix(data[, c("Tmed", "Tmed1","Tmed2", "Tmed3", "Tmed4", "Tmed5", "Tmed6")])

# hum average 
data$RH06 <- rowMeans(data[, c("Hmed", "Hmed1", "Hmed2", "Hmed3", "Hmed4", "Hmed5", "Hmed6")])

# crossbasis
cb06 <- crossbasis(Tmed06,lag=6,
                   argvar=list(fun = "ns", df = 2), # spline with 1 knot (2 df) in the temperature dimension
                   arglag=list(fun="ns", knots=logknots(6, 1),int=T)) # spline with 1 knot (2 df) in the lag dimension

# model
mod06 <- clogit(bin_apgar_numeric ~ cb06 + ns(RH06, df=3) + strata(unique_id), data)

# summary
summary(mod06)

# predictions
cp06 <- crosspred(cb06 , mod06, cen = pop_quantiles[2] , by=0.1)


# table results together
table_sensitivity <- data.frame(
  Object = c("0",   "01", "01", "06"),
  OR_Lag01 = c(0, cp01a$allRRfit["26.1"], cp01b$allRRfit["26.1"], cp06$allRRfit["26.1"]),
  CI_Low_Lag01 = c(0, cp01a$allRRlow["26.1"], cp01b$allRRlow["26.1"],cp06$allRRlow["26.1"]),
  CI_High_Lag01 = c(0, cp01a$allRRhigh["26.1"], cp01b$allRRhigh["26.1"],  cp06$allRRhigh["26.1"]),
  OR_Lag0 = c(cp0$allRRfit["26.1"], cp01a$matRRfit["26.1", "lag0"], cp01b$matRRfit["26.1", "lag0"], cp06$matRRfit["26.1", "lag0"]),
  CI_Low_Lag0 = c(cp0$allRRlow["26.1"], cp01a$matRRlow["26.1", "lag0"], cp01b$matRRlow["26.1", "lag0"], cp06$matRRlow["26.1", "lag0"]),
  CI_High_Lag0 = c(cp0$allRRhigh["26.1"], cp01a$matRRhigh["26.1", "lag0"], cp01b$matRRhigh["26.1", "lag0"],cp06$matRRhigh["26.1", "lag0"]),
  OR_Lag1 = c(0, cp01a$matRRfit["26.1", "lag1"], cp01b$matRRfit["26.1", "lag1"],cp06$matRRfit["26.1", "lag1"]),
  CI_Low_Lag1 = c(0, cp01a$matRRlow["26.1", "lag1"], cp01b$matRRlow["26.1", "lag1"],cp06$matRRlow["26.1", "lag1"]),
  CI_High_Lag1 = c(0, cp01a$matRRhigh["26.1", "lag1"], cp01b$matRRhigh["26.1", "lag1"],cp06$matRRhigh["26.1", "lag1"]))

table_sensitivity$estimate_ci_01 <- sprintf("%.2f (%.2f, %.2f)", table_sensitivity$OR_Lag01, table_sensitivity$CI_Low_Lag01, table_sensitivity$CI_High_Lag01)
table_sensitivity$estimate_ci_0 <- sprintf("%.2f (%.2f, %.2f)", table_sensitivity$OR_Lag0, table_sensitivity$CI_Low_Lag0, table_sensitivity$CI_High_Lag0)
table_sensitivity$estimate_ci_1 <- sprintf("%.2f (%.2f, %.2f)", table_sensitivity$OR_Lag1, table_sensitivity$CI_Low_Lag1, table_sensitivity$CI_High_Lag1)

print(table_sensitivity)






