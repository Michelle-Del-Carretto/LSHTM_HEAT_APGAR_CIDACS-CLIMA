############################################################

# Analysis from:
# 'Acute heat exposure and low APGAR score in newborns: A time-stratified case-crossover analysis 
# in São Paulo state, Brazil'

# Date of code creation:
# 10/2/2025

# Code author:
# Michelle Del Carretto

# Purpose:
# 1) Clean and setup dataset
# 2) Calculate population-weighted percentiles (whole dataset)
# 3) Calculate population-weighted percentiles (Koppen zones)

############################################################
# Load datasets #
############################################################

setwd("~/Documents/LSHTM project/Updated exposures")
load("Updated_exposures.RData")
setwd("~/Documents/LSHTM project/R code")
birth_count <- read.csv("GIT_ts_dataset.csv")

############################################################
# Format datasets #
############################################################

# check all are in correct format
unique(combined$date[is.na(as.Date(combined$date, format="%d/%m/%Y")) & !is.na(combined$date)]) #character(0)
summary(combined$date)

# keep observations between 1/1/2013 (and the week before) and 31/12/2019 
start_date <- as.Date("24-12-2012", format="%d-%m-%Y") #2 extra days to be safe as.Date("24-12-2012", format="%d-%m-%Y")
end_date <- as.Date("31-12-2019", format="%d-%m-%Y")
df_temp_updated <- combined[combined$date >= start_date & combined$date <= end_date, ]
summary(df_temp_updated$date)

# summarise dataset
summary(df_temp_updated$Tmed) # 4.80   20.20   22.80   22.43   25.00   34.30 
summary(df_temp_updated$Hmed) # 19.10   61.70   72.70   70.54   81.50  100.00
quantile(df_temp_updated$Tmed, c(.05, .50, .95)) # 16.0  22.8  27.6  


############################################################
# Calculate population-weighted percentiles #
############################################################

# Format timeseries and exposures datasets 
names(birth_count)[names(birth_count) == 'dob'] <- 'date' # rename column
names(birth_count)[names(birth_count) == 'birth_municip'] <- 'municipality_code' # rename column
birth_count$date <- as.Date(birth_count$date) # ensure correct format
birth_count$municipality_code <- as.character(birth_count$municipality_code) # ensure correct format
df_temp_updated$municipality_code <- substr(df_temp_updated$municipality_code, 1, 
                                            nchar(df_temp_updated$municipality_code) - 1) # remove last character from municipality code in combined exposures dataset

# Calculate total births per municipality over study period
tbirths <- aggregate(birth_count ~ municipality_code, data = birth_count, FUN = sum)

# keep koppen zone
koppen_info <- aggregate(Koppen ~ municipality_code, data = birth_count, FUN = function(x) x[1])
tbirths <- merge(tbirths, koppen_info, by = "municipality_code", all.x = TRUE)

# Merge with daily Tmed data
Tmed2 <- merge(df_temp_updated, tbirths, by = "municipality_code", all.x = TRUE)

# Replace NA values with 0 in tbirths column
Tmed2$birth_count[is.na(Tmed2$birth_count)] <- 0

# Extract temperatures and weights
temperatures <- Tmed2$Tmed
weights <- Tmed2$birth_count

# Sort temperatures and corresponding weights
sorted_indices <- order(temperatures)
sorted_temperatures <- temperatures[sorted_indices]
sorted_weights <- weights[sorted_indices]

# Compute cumulative sum of weights, normalized
cumulative_weights <- cumsum(sorted_weights) / sum(sorted_weights)

# Define percentiles of interest
percentiles <- c(0.05, 0.50, 0.95)

# Find temperature values corresponding to percentiles
weighted_percentiles <- sapply(percentiles, function(p) sorted_temperatures[min(which(cumulative_weights >= p))])

# Print weighted percentiles
print(weighted_percentiles) # 14.8 20.9 26.1


############################################################
# Calculate population-weighted percentiles (Köppen zones) #
############################################################

# create koppen broad category (tropical vs temperate)
Tmed2$koppen_broad <- substr(Tmed2$Koppen, 1, 1)

# boxplot temperature for each koppen zone
boxplot(Tmed ~ Koppen, data = Tmed2, 
        main = "Boxplot of Tmed by Koppen Zone", 
        xlab = "Koppen Zone", 
        ylab = "Tmed", 
        col = "lightblue")
boxplot(Tmed ~ koppen_broad, data = Tmed2, 
        main = "Boxplot of Tmed by Koppen Zone", 
        xlab = "Koppen Zone", 
        ylab = "Tmed", 
        col = "lightblue")

# Function to calculate weighted quantiles
weighted_quantiles <- function(temperatures, weights, percentiles) {
  # Sort temperatures and corresponding weights
  sorted_indices <- order(temperatures)
  sorted_temperatures <- temperatures[sorted_indices]
  sorted_weights <- weights[sorted_indices]
  # Compute cumulative weights
  cumulative_weights <- cumsum(as.numeric(sorted_weights)) / sum(as.numeric(sorted_weights))
  # Find weighted percentiles
  percentile_indices <- sapply(percentiles, function(p) min(which(cumulative_weights >= p)))
  sorted_temperatures[percentile_indices]}

# Define percentiles
percentiles <- c(0.05, 0.50, 0.95)

# Split the data by 'koppen_broad' and apply the weighted_quantiles function to each subset
results <- lapply(split(Tmed2, Tmed2$koppen_broad), function(Tmed2) {
  # Extract temperatures and weights for the current subset
  temperatures <- Tmed2$Tmed
  weights <- Tmed2$birth_count
  # Calculate the weighted quantiles
  quantiles <- weighted_quantiles(temperatures, weights, percentiles)
  # Return the results in a data frame for each koppen_broad
  data.frame(
    koppen_broad = unique(Tmed2$koppen_broad),
    Q05 = quantiles[1],
    Q50 = quantiles[2],
    Q95 = quantiles[3])})

# Print results
print(results)

# A (tropical)
# Q05  Q50   Q95
# 17.4 23.6  28

# C (temperate)
#Q05  Q50  Q95
#14.6 20.6 25.6


