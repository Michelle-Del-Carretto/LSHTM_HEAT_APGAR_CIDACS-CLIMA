############################################################

# Analysis from:
# 'hort-term ambient heat exposure and low APGAR score in newborns: A time-stratified case-crossover analysis 
# in São Paulo state, Brazil'

# Date of code creation:
# 14/2/2025

# Code author:
# Michelle Del Carretto

# Purpose:
# 1) Format timeseries datasets
# 2) Visualise seasonal patterns in rate of low APGAR
# 3) Unadjusted relationships between temperature, relative humidity and low APGAR (lag 0 & 1)

############################################################
# Load datasets & packages #
############################################################

setwd()
load("Updated_exposures copy.RData")
birth_count <- read.csv("GIT_ts_dataset.csv")

library(ggplot2)

############################################################
# Format datasets #
############################################################

# format exposure dataset
df_exp <- subset(combined, date >= as.Date("31-12-2012", format="%d-%m-%Y") & date <= as.Date("31-12-2019", format="%d-%m-%Y")) # trim to desired dates
df_exp$municipality_code <- as.character(df_exp$municipality_code) # convert to character
df_exp$municipality_code <- substr(df_exp$municipality_code, 1, nchar(df_exp$municipality_code) - 1) # remove last character in string
df_exp$date = as.Date(df_exp$date, format = "%Y-%m-%d") # set date class

# format health dataset
names(birth_count)[names(birth_count) == 'dob'] <- 'date' # rename date column
names(birth_count)[names(birth_count) == 'birth_municip'] <- 'municipality_code'# rename municip column
birth_count$municipality_code <- as.character(birth_count$municipality_code) # convert to character
birth_count$date = as.Date(birth_count$date, format = "%Y-%m-%d") # set date class

# join datasets
data <- merge(birth_count, df_exp, by = c("municipality_code", "date"), all.x = TRUE)
data$doy <- as.numeric(format(data$date, "%j")) # create new day-of-year column

### create aggregated dataset 1 (rate of low APGAR for every day in period 2013-2019)
aggregated_data <- transform(
  aggregate(cbind(low_apgar_count, birth_count) ~ date,  # aggregate data (summing) by date
            data = data, 
            FUN = sum, na.rm = TRUE),
  rate = low_apgar_count / birth_count) # calculate rate
aggregated_data$Tmed <- tapply(data$Tmed, data$date, FUN = mean, na.rm = TRUE)
aggregated_data$Hmed <- tapply(data$Hmed, data$date, FUN = mean, na.rm = TRUE)
# Create lagged versions of Tmed and Hmed (previous day values)
aggregated_data$Tmed_lag1 <- c(NA, head(aggregated_data$Tmed, -1))  # Shift Tmed values by 1 row
aggregated_data$Hmed_lag1 <- c(NA, head(aggregated_data$Hmed, -1))  # Shift Hmed values by 1 row

### create aggregated dataset 2 (rate of low APGAR by day of year)
# for APGAR
aggregated_data2 <- transform(
  aggregate(cbind(low_apgar_count, birth_count) ~ doy,# aggregate data (summing) by day-of-year
            data = data, 
            FUN = sum, na.rm = TRUE),
  rate = low_apgar_count / birth_count) # calculate rate

# for temperature
municips_represented <- unique(data$municipality_code) # Get unique municipality codes from health data
filtered_exp_df <- df_exp[df_exp$municipality_code %in% municips_represented, ] # Filter the exposure dataframe to include only represented municipalities
filtered_exp_df$doy <- as.numeric(format(filtered_exp_df$date, "%j")) # add day-of-year column
aggregated_temp <- aggregate(Tmed ~ doy, data = filtered_exp_df, FUN = mean, na.rm = TRUE) # Aggregate temperature by day-of-year

# join temperature and health
aggregated_data2 <- merge(aggregated_data2, aggregated_temp, by = "doy")


############################################################
# Seasonal patterns #
############################################################

# timeseries plot
ggplot()+
  geom_point(data = aggregated_data, aes(x = date, y = rate *100), size = 1.5, alpha = 0.6) +
  ylim(c(0, 3)) +
  geom_smooth(data = aggregated_data, aes(x = date, y = rate * 100), method = "loess", se = FALSE, color = "red", linetype = "solid") +
  labs(title = "", x = "Year", y = "Mean low APGAR (≤7) rate (%)") +
  theme_minimal() +  
  theme(
    panel.grid.major = element_line(color = "gray90"),  
    panel.grid.minor = element_line(color = "gray80", size = 0.5),  
    minor_breaks = seq(min(aggregated_data$date), max(aggregated_data$date), by = "1 year"),  # Set minor gridline breaks, adjust 'by' to control frequency
    plot.margin = margin(10, 10, 10, 10)  # Add margins for better spacing
  ) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year")  # Format x-axis labels for date and set breaks


# preparation for day-of-year plots (Brazilian summer in the middle)
# Define month breakpoints and labels
month_days <- cumsum(c(0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)) - 182 # list of total days in each month and then 182 days subtracted from each value. 
                                                                                 # This shifts the dates so that July (the beginning of the Brazilian winter) aligns around the middle of the plot
month_midpoints <- month_days[-13] + diff(month_days) / 2 #calculating the midpoint of each month to use as x-axis ticks
month_names <- c("Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun")
summershade <- data.frame(xmin=-30,xmax=61) # create shading to indicate start and end of summer season

# ggplot
a<-ggplot(aggregated_data2, aes(x = (doy - 1 + 182) %% 365 - 182, y = rate * 100)) + # shifts the data to center it around the Brazilian summer
  geom_point(size = 2) +
  geom_smooth(method = "loess", se = TRUE, col = "royalblue", size = 1.5) +
  ylim(c(0.6, 1.6)) +
  geom_rect(data = summershade, aes(xmin = xmin, xmax = xmax), 
            ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.1, inherit.aes = FALSE) +
  labs(
    title = "",
    x = "Month",
    y = "Mean low APGAR (≤7) rate (%)"
  ) +
  scale_x_continuous(
    breaks = month_midpoints, # Centered labels for each month
    labels = month_names,     # Corresponding month names
    limits = c(-182, 182)     # Center January
  ) +
  Theme

# 
b<-ggplot(aggregated_data2, aes(x = (doy - 1 + 182) %% 365 - 182, y = Tmed)) +
  geom_point(size = 2) +
  geom_smooth(method = "loess", se = TRUE, col = "royalblue", size = 1.5) +
  ylim(c(16, 26)) +
  geom_rect(data = summershade, aes(xmin = xmin, xmax = xmax), 
            ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.1, inherit.aes = FALSE) +
  labs(
    title = "",
    x = "Month",
    y = "Mean temperature (°C)"
  ) +
  scale_x_continuous(
    breaks = month_midpoints, # Centered labels for each month
    labels = month_names,     # Corresponding month names
    limits = c(-182, 182)     # Center January
  ) +
  Theme

#combined_plot <- a + b + plot_annotation(tag_levels = 'A')
#combined_plot
#ggsave("xxx.png", plot = combined_plot, width = 14, height = 6, dpi = 300)

############################################################
# Unadjusted relationships #
############################################################

# Temperature-humidity
ggplot(aggregated_data, aes(Tmed, Hmed)) +
  geom_point() +
  geom_smooth() +
  labs(title = "",
       x = expression(paste("Daily mean temperature (" ,degree,"C)")), 
       y = "Relative humidity (%)")

# Temperature-Low APGAR (lag 0)
ggplot(aggregated_data, aes(Tmed, rate*100)) +
  geom_point() +
  geom_smooth() +
  labs(title = "",
       x = expression(paste("Daily mean temperature (" ,degree,"C)")), 
       y = "Low Apgar Rate") +
  theme_minimal()

# Temperature-Low APGAR (lag 1)
ggplot(aggregated_data, aes(Tmed_lag1, rate*100)) +
  geom_point() +
  geom_smooth() +
  labs(title = "",
       x = expression(paste("Daily mean temperature (" ,degree,"C)")), 
       y = "Low Apgar Rate") +
  theme_minimal()

# Humidity-Low APGAR (lag 0)
ggplot(aggregated_data, aes(Hmed, rate*100)) +
  geom_point() +
  geom_smooth() +
  labs(title = "",
       x = expression(paste("Daily mean humidity (" ,degree,"C)")), 
       y = "Low Apgar Rate") +
  theme_minimal()

# Humidity-Low APGAR (lag 1)
ggplot(aggregated_data, aes(Hmed_lag1, rate*100)) +
  geom_point() +
  geom_smooth() +
  labs(title = "",
       x = expression(paste("Daily mean temperature (" ,degree,"C)")), 
       y = "Low Apgar Rate") +
  theme_minimal()




