############################################################

# Analysis from:
# 'Acute heat exposure and low APGAR score in newborns: A time-stratified case-crossover analysis 
# in São Paulo state, Brazil'

# Date of code creation:
# 10/2/2025

# Code author:
# Michelle Del Carretto

# Purpose:
# 1) Format health dataset (create controls)
# 2) Format exposure dataset (create lagged values)
# 3) Join and save dataset for for Case-crossover analysis

############################################################
# Load datasets #
############################################################

setwd()
load("Updated_exposures copy.RData")
low_apgar_df <- read.csv("GIT_health_dataset.csv")

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
  low_apgar_dt[, .(unique_id, date, municipality_code,clean_apgar,sex, bin_parity_2, delivery_mod, emer_c_sec, combined_delivery_mod, 
                   bin_prenatal_care, bin_mage, ethnicity, bin_edu, q_IBP, Koppen, # select relevant columns
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

# save
#write.csv(data, file = "GIT_CCdataset.csv", row.names = FALSE)
