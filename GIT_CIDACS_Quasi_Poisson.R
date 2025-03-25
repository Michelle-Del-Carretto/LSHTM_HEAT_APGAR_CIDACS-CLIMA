############################################################

# Analysis from:
# 'Acute heat exposure and low APGAR score in newborns: A time-stratified case-crossover analysis 
# in São Paulo state, Brazil'

# Date of code creation:
# 14/2/2025

# Code author:
# Michelle Del Carretto

# Purpose:
# 1) Format complete timeseries dataset (append exposures)
# 2) Specify & run conditional Quasi-Poisson

############################################################
# Load datasets & packages #
############################################################

setwd()
load("Updated_exposures copy.RData")
complete_ts <- read.csv("GIT_complete_ts_dataset.csv")

library(gnm)
library(dlnm)
library(splines)

############################################################
# Format complete timeseries #
# (series of dates with low APGAR count & birth count in every
# municipality with at least 1 low APGAR between 2013-2019)
############################################################

# format exposure dataset
df_exp <- subset(combined, date >= as.Date("31-12-2012", format="%d-%m-%Y") & date <= as.Date("31-12-2019", format="%d-%m-%Y")) # trim to desired dates
df_exp$municipality_code <- as.character(df_exp$municipality_code) # convert to character
df_exp$municipality_code <- substr(df_exp$municipality_code, 1, nchar(df_exp$municipality_code) - 1) # remove last character in string
df_exp$date = as.Date(df_exp$date, format = "%Y-%m-%d") # set date class

# format health dataset
names(complete_ts)[names(complete_ts) == 'dob'] <- 'date' # rename date column
names(complete_ts)[names(complete_ts) == 'birth_municip'] <- 'municipality_code'# rename municip column
complete_ts$municipality_code <- as.character(complete_ts$municipality_code) # convert to character
complete_ts$date = as.Date(complete_ts$date, format = "%Y-%m-%d") # set date class

# create a rolling average of humidity (lag 0 & lag 1) and append to dataset
df_exp <- df_exp[order(df_exp$municipality_code, df_exp$date), ] # Sort the data by municipality_code and date to ensure proper ordering
df_exp$Hmed_lag1 <- ave(df_exp$Hmed, df_exp$municipality_code, FUN = function(x) c(NA, head(x, -1))) # Create a lagged version of the Hmed column for each municipality
df_exp$Hmed01 <- (df_exp$Hmed + df_exp$Hmed_lag1) / 2 # Calculate the rolling average as the mean of the current and previous value within each municipality

# join exposure and complete timeseries data
complete_ts <- merge(complete_ts, df_exp, by = c("date", "municipality_code"), all.x = TRUE)

# Set the default to na.exclude 
options(na.action="na.exclude")

############################################################
# Specify & run conditional Quasi-Poisson #
############################################################

# specify strata (similar to case-crossover where temperature is compared between sets of days matched on location, month and day-of-week)
complete_ts$stratum <- factor(paste(complete_ts$municipality_code, complete_ts$year, complete_ts$month, complete_ts$dow, sep=":"))

# specify exposure-response relationship using a crossbasis
cb <- crossbasis(complete_ts$Tmed,lag=1,
                 argvar=list(fun = "ns", df = 2),
                 arglag=list(fun="lin",int=T), 
                 group=complete_ts$municipality_code) # group by municipality such that locations can be discerned

# create indicator for whether strata are empty or have at least 1 low APGAR
complete_ts$keep <- as.logical(ave(complete_ts$low_apgar_count, complete_ts$stratum, FUN = function(x) sum(x) > 0))

# run model
mod_qp <- gnm(low_apgar_count ~ cb + ns(Hmed01, df=3), 
              eliminate=stratum, # exploits strata for calculation purposes
              data=complete_ts, 
              offset=log(birth_count + 0.0001),  # allows adjustment for variable birth rate
              family=quasipoisson, # allows adjustment for overdispersion
              subset=keep) # only performs calculation on strata that are non-0 (i.e. strata with 0 observations are eliminated)

# summary
summary(mod_qp)

# PACF
residuals <- residuals(mod_qp, type = "pearson") #anything below 0.1 or 0.2 is usually fine
pacf(residuals, lag.max = 15, plot = TRUE, na.action = na.pass, main = "Partial autocorrelation among residuals of DLNM")

# get predictions
cp <- crosspred(cb, mod_qp, cen = 20.9, by=0.1)
cp$allRRfit["26.1"]
cp$allRRlow["26.1"]
cp$allRRhigh["26.1"]
cp$matRRfit["26.1", "lag0"]
cp$matRRlow["26.1", "lag0"]
cp$matRRhigh["26.1", "lag0"]
cp$matRRfit["26.1", "lag1"]
cp$matRRlow["26.1", "lag1"]
cp$matRRhigh["26.1", "lag1"]

# plot
plot(cp , "overall", ylab="RR", main="Overall week (ns2-ns1 CC)",xlab=expression(paste("Daily mean temperature (",degree,"C)")),cex.axis=1.2, cex.lab=1.2, lwd=2.5, ylim = c(0.7,1.3),ci.arg=list(density=100,col=grey(0.7)))
abline(v = c(14.8, 20.9, 26.1), col = c("blue", "black", "red"), lty = "dashed", lwd = 2)
plot(cp, "slices", lag=0, ylab="RR", main= "Lag 0 (ns2-ns1 CC)",xlab=expression(paste("Daily mean temperature (",degree,"C)")),cex.axis=1.2, cex.lab=1.2, lwd=2.5, ylim = c(0.7,1.3),ci.arg=list(density=100,col=grey(0.7)))
abline(v = c(14.8, 20.9, 26.1), col = c("blue", "black", "red"), lty = "dashed", lwd = 2)
plot(cp, "slices", lag=1, ylab="RR", main= "Lag 0 (ns2-ns1 CC)",xlab=expression(paste("Daily mean temperature (",degree,"C)")),cex.axis=1.2, cex.lab=1.2, lwd=2.5, ylim = c(0.7,1.3),ci.arg=list(density=100,col=grey(0.7)))
abline(v = c(14.8, 20.9, 26.1), col = c("blue", "black", "red"), lty = "dashed", lwd = 2)

