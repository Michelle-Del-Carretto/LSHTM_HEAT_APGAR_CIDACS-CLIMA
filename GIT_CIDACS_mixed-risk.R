
############################################################

# Analysis from:
# 'Acute heat exposure and low APGAR score in newborns: A time-stratified case-crossover analysis 
# in São Paulo state, Brazil'

# Date of code creation:
# 23/3/2025

# Code author:
# Michelle Del Carretto

# Purpose:
# 1) Clean and setup health dataset (less restricted - mixed-risk births) 
# 4) Format dataset for case-crossover analysis
# 5) Attach exposure data
# 6) Case-crossover analyisis

############################################################
# Load datasets #
############################################################

# set working directory
setwd("~/Documents/LSHTM project/R code")

#load dataset 
load("/Users/michelledelcarretto/Documents/LSHTM project/R code/DNSP_2010_2019 copy.RData")


############################################################
# Clean and setup health dataset #
############################################################

###### Date of birth ######
# Create new date of birth column
combined_dataframe$dob <- combined_dataframe$DTNASC  

# Check 8 characters long & correct format
dob_char <- as.character(combined_dataframe$dob)
all(nchar(dob_char) == 8) # TRUE  
unique(combined_dataframe$dob[is.na(as.Date(combined_dataframe$dob, format="%d%m%Y")) & !is.na(combined_dataframe$dob)]) # factor(0)
rm(dob_char)

# specify date class
combined_dataframe$dob <- as.Date(combined_dataframe$dob, format="%d%m%Y")

# restrict: 2013-2019
combined_dataframe <- combined_dataframe[combined_dataframe$dob >= as.Date("2013-01-01"), ]


###### Singleton ######
summary(combined_dataframe$GRAVIDEZ)

# keep only level "1"
combined_dataframe <- combined_dataframe[combined_dataframe$GRAVIDEZ == "1" & !is.na(combined_dataframe$GRAVIDEZ), ]


###### Birth municipality ######
summary(combined_dataframe$CODMUNNASC)

# Convert to character
combined_dataframe$birth_municip <- as.character(combined_dataframe$CODMUNNASC)

# Keep only rows starting with "35" (i.e. in Sao Paulo state)
combined_dataframe <- combined_dataframe[grepl("^35", combined_dataframe$birth_municip), ]

# Convert back to factor
combined_dataframe$birth_municip <- as.factor(combined_dataframe$birth_municip)

# Check 
nlevels(combined_dataframe$birth_municip)  # 524 municipalities left 


###### APGAR 5 ######
summary(combined_dataframe$APGAR5)

# convert to numeric & remove NAs and 99
combined_dataframe$clean_apgar <- as.numeric(as.character(as.character(combined_dataframe$APGAR5)))
combined_dataframe <- combined_dataframe[!is.na(combined_dataframe$clean_apgar) & combined_dataframe$clean_apgar != "99", ]

# create variable with APGAR categories
combined_dataframe$clean_apgar <- factor(ifelse(combined_dataframe$clean_apgar %in% 0:2, 1,
                                                ifelse(combined_dataframe$clean_apgar %in% 3:5, 2,
                                                       ifelse(combined_dataframe$clean_apgar %in% 6:7, 3, 4))), 
                                         levels = c(1, 2, 3, 4), labels = c("0-2", "3-5", "6-7", "8-10"))
# Check 
summary(combined_dataframe$clean_apgar)

# Create binary APGAR variable
combined_dataframe$bin_apgar <- factor(ifelse(combined_dataframe$clean_apgar %in% c("0-2", "3-5", "6-7"), 1, 0), 
                                       levels = c(1, 0), 
                                       labels = c("Low Apgar", "Normal Apgar"))
# Check 
summary(combined_dataframe$bin_apgar)

###########################################################
# Save dataset (APGAR ≤7) #
############################################################

# filter for low Apgar births
low_apgar_df <- combined_dataframe[combined_dataframe$bin_apgar == "Low Apgar", ]

# Select specific columns
low_apgar_df <- low_apgar_df[, c("birth_municip", "dob", "clean_apgar", "bin_apgar" )]

# reoder by date of birth
low_apgar_df <- low_apgar_df[order(low_apgar_df$dob), ]

rm(combined_dataframe)

############################################################
# Load datasets #
############################################################

setwd("~/Documents/LSHTM project/Updated exposures")
load("Updated_exposures.RData")

setwd("~/Documents/LSHTM project/R code")

library(data.table)

############################################################
# Format health dataset #
############################################################


# rename date & municip columns
names(low_apgar_df)[names(low_apgar_df) == 'dob'] <- 'date'
names(low_apgar_df)[names(low_apgar_df) == 'birth_municip'] <- 'municipality_code'

# add ID
low_apgar_df$unique_id <- 1:nrow(low_apgar_df)

# first, create dataframe of all days within date range
range(low_apgar_df$date) # date range
alldays <- data.frame(date = seq(as.Date("2013-01-01"), # min date
                                 as.Date("2019-12-31"), # max date
                                 by = 'days'))

# Add year, month, and day of the week (dow) to dataframe of all days
alldays$year <- as.integer(format(alldays$date, "%Y"))
alldays$month <- as.integer(format(alldays$date, "%m"))
alldays$dow <- factor(weekdays(alldays$date))

# Add year, month, and dow to the low_apgar_df
low_apgar_df$date <- as.Date(low_apgar_df$date)
low_apgar_df$year <- as.integer(format(low_apgar_df$date, "%Y"))
low_apgar_df$month <- as.integer(format(low_apgar_df$date, "%m"))
low_apgar_df$dow <- factor(weekdays(low_apgar_df$date))


# convert alldays and cases dataframes to data tables 
setDT(alldays)
low_apgar_dt <- setDT(low_apgar_df)

# insert new rows through a Cartesian join 
low_apgar_dt <- alldays[
  low_apgar_dt[, .(unique_id, date, municipality_code,clean_apgar, 
                    bin_apgar = "Normal Apgar",  # Assign default value of "Normal Apgar"
                   year, month, dow)],
  on = .(year, month, dow), allow.cartesian = TRUE # match control rows on year, month, dow
][
  low_apgar_dt, on = .(unique_id, date), `:=`(bin_apgar = i.bin_apgar) #join with the original dataset (If i.bin_apgar -the bin_apgar value from the original low_apgar_dt- exists, it replaces "Normal Apgar" with that value)
][]

# convert to dataframe and rename 
data <- as.data.frame(low_apgar_dt)

# Recoding factor levels to 0 and 1 
data$bin_apgar_numeric <- ifelse(data$bin_apgar == "Low Apgar", 1, 
                                 ifelse(data$bin_apgar == "Normal Apgar", 0, NA))
data$bin_apgar_numeric <- as.numeric(data$bin_apgar_numeric)

# On average, how many controls per ID
summary_data <- aggregate(cbind(total_cases = bin_apgar_numeric, total_controls = 1 - bin_apgar_numeric) ~ unique_id, 
                          data = data, FUN = sum) # Group by unique_id and calculate sum of cases and controls
mean_data <- data.frame(                          # Calculate the mean number of cases and controls
  mean_cases = mean(summary_data$total_cases),
  mean_controls = mean(summary_data$total_controls))
mean_data

rm(alldays)
rm(low_apgar_dt)
rm(low_apgar_df)
rm(summary_data)
rm(mean_data)

############################################################
# Format & link exposure data #
############################################################

# keep observations between 1/1/2013 (and the week before) and 31/12/2019 
start_date <- as.Date("24-12-2012", format="%d-%m-%Y") #2 extra days to be safe as.Date("24-12-2012", format="%d-%m-%Y")
end_date <- as.Date("31-12-2019", format="%d-%m-%Y")
df_exp <- combined[combined$date >= start_date & combined$date <= end_date, ]
summary(df_exp$date)

# remove last character from municipality in exposure dataset
df_exp$municipality_code <- substr(df_exp$municipality_code, 1, nchar(df_exp$municipality_code) - 1)
municipality_code_char <- as.character(df_exp$municipality_code)

# set date class
df_exp$date = as.Date(df_exp$date, format = "%Y-%m-%d")
df_exp$municipality_code <- as.character(df_exp$municipality_code)

# Sort the data by municipality_code and date
df_exp <- df_exp[order(df_exp$municipality_code, df_exp$date), ]

# Use ave() to calculate lagged values of Tmed
df_exp$Tmed1 <- ave(df_exp$Tmed, df_exp$municipality_code, FUN = function(x) c(NA, head(x, -1)))
df_exp$Tmed2 <- ave(df_exp$Tmed, df_exp$municipality_code, FUN = function(x) c(NA, NA, head(x, -2)))
df_exp$Tmed3 <- ave(df_exp$Tmed, df_exp$municipality_code, FUN = function(x) c(NA, NA, NA, head(x, -3)))
df_exp$Tmed4 <- ave(df_exp$Tmed, df_exp$municipality_code, FUN = function(x) c(NA, NA, NA, NA, head(x, -4)))
df_exp$Tmed5 <- ave(df_exp$Tmed, df_exp$municipality_code, FUN = function(x) c(NA, NA, NA, NA, NA, head(x, -5)))
df_exp$Tmed6 <- ave(df_exp$Tmed, df_exp$municipality_code, FUN = function(x) c(NA, NA, NA, NA, NA, NA, head(x, -6)))


# Use ave() to calculate lagged values of Hmed
df_exp$Hmed1 <- ave(df_exp$Hmed, df_exp$municipality_code, FUN = function(x) c(NA, head(x, -1)))
df_exp$Hmed2 <- ave(df_exp$Hmed, df_exp$municipality_code, FUN = function(x) c(NA, NA, head(x, -2)))
df_exp$Hmed3 <- ave(df_exp$Hmed, df_exp$municipality_code, FUN = function(x) c(NA, NA, NA, head(x, -3)))
df_exp$Hmed4 <- ave(df_exp$Hmed, df_exp$municipality_code, FUN = function(x) c(NA, NA, NA, NA, head(x, -4)))
df_exp$Hmed5 <- ave(df_exp$Hmed, df_exp$municipality_code, FUN = function(x) c(NA, NA, NA, NA, NA, head(x, -5)))
df_exp$Hmed6 <- ave(df_exp$Hmed, df_exp$municipality_code, FUN = function(x) c(NA, NA, NA, NA, NA, NA, head(x, -6)))

# join datasets by date and region
data$municipality_code <- as.character(data$municipality_code)
data <- merge(data, df_exp, by = c("date", "municipality_code"), all.x = TRUE)

# Sort the data by unique_id
data <- data[order(data$unique_id), ]

rm(combined)
rm(df_exp)
rm(end_date)
rm(start_date)
rm(municipality_code_char)

############################################################
# Load datasets #
############################################################

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
