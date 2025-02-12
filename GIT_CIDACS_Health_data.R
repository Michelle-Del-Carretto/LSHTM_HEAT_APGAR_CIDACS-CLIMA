
############################################################

# Analysis from:
# 'Acute heat exposure and low APGAR score in newborns: A time-stratified case-crossover analysis 
# in São Paulo state, Brazil'

# Date of code creation:
# 3/2/2025

# Code author:
# Michelle Del Carretto

# Purpose:
# 1) Clean and setup health dataset
# 2) Attach municipality-level metrics (IBP, Köppen) 
# 3) Descriptive statistics (counts, percentages and Chi-2)
# 4) Save dataset for case-crossover analysis
# 5) Format & save time series dataset

############################################################
# Load datasets #
############################################################

# set working directory
setwd("~/Documents/LSHTM project/R code")

#load dataset 
load("/Users/michelledelcarretto/Documents/LSHTM project/R code/DNSP_2010_2019 copy.RData")
df_ibp <- read.csv("BDI_Municipalities-Level_Short.csv")
df_koppen <- read.csv("Koppen_municipalities_distinct2.csv")


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


###### Birthweight ######
# convert to numeric
combined_dataframe$b.weight <- as.numeric(as.character(combined_dataframe$PESO))

# restrict to weights between 2500 and 4000
combined_dataframe <- combined_dataframe[combined_dataframe$b.weight >= 2500 & combined_dataframe$b.weight <= 4000 & !is.na(combined_dataframe$b.weight), ]


###### Gestational age ######
summary(combined_dataframe$GESTACAO)

# keep only level "5" (37-41 weeks)
combined_dataframe <- combined_dataframe[combined_dataframe$GESTACAO == "5" & !is.na(combined_dataframe$GESTACAO), ]


###### Congenital anomalies ######
summary(combined_dataframe$IDANOMAL)

# keep only level "2" (no congenital anomalies)
combined_dataframe <- combined_dataframe[combined_dataframe$IDANOMAL == "2" & !is.na(combined_dataframe$IDANOMAL), ]


###### Birth presentation ######
summary(combined_dataframe$TPAPRESENT)

# keep only level "1" (cephalic)
combined_dataframe <- combined_dataframe[combined_dataframe$TPAPRESENT == "1" & !is.na(combined_dataframe$TPAPRESENT), ]


###### Infant sex ######
summary(combined_dataframe$SEXO)

# create new column
combined_dataframe$sex <- as.character(combined_dataframe$SEXO)

# 0, 9, I to NA & reassign M and F to 1 and 2
combined_dataframe$sex[combined_dataframe$sex %in% c("9", "I", "0")] <- NA
combined_dataframe$sex[combined_dataframe$sex == "M"] <- "1"
combined_dataframe$sex[combined_dataframe$sex == "F"] <- "2"
combined_dataframe$sex <- factor(combined_dataframe$sex, levels = c("1", "2"), labels = c("Male", "Female"))
# check
summary(combined_dataframe$sex)


###### Parity (previous vaginal and c-section births are combined) ######
summary(combined_dataframe$QTDPARTNOR)
summary(combined_dataframe$QTDPARTCES)

# Vaginal births - remove '99' and 'NA'
combined_dataframe$v_births <- combined_dataframe$QTDPARTNOR
combined_dataframe$v_births[combined_dataframe$v_births == "99"] <- NA

# Caesarian births - remove '99' and 'NA'
combined_dataframe$c_births <- combined_dataframe$QTDPARTCES
combined_dataframe$c_births[combined_dataframe$c_births == "99"] <- NA

# Combine vaginal and c-section births to calculate parity
combined_dataframe$parity <- as.numeric(as.character(combined_dataframe$v_births)) + 
  as.numeric(as.character(combined_dataframe$c_births))

# Remove all observations above 10
combined_dataframe$parity[combined_dataframe$parity > 10] <- NA

# Create 3-category variable
combined_dataframe$bin_parity_2 <- ifelse(combined_dataframe$parity == 0, 0, 
                                          ifelse(combined_dataframe$parity == 1, 1, 2))

# Convert to factor with labels
combined_dataframe$bin_parity_2 <- factor(combined_dataframe$bin_parity_2, 
                                          levels = c(0, 1, 2), 
                                          labels = c("Nulliparous", "Primiparous", "Multiparous"))
# Check
summary(combined_dataframe$bin_parity_2)


###### Delivery mode (vaginal, c-section) ######
summary(combined_dataframe$PARTO)

# Convert 9 to NA
combined_dataframe$delivery_mod <- combined_dataframe$PARTO
combined_dataframe$delivery_mod[combined_dataframe$delivery_mod == 9] <- NA

# label factor
combined_dataframe$delivery_mod <- factor(combined_dataframe$delivery_mod, levels = c(1, 2), labels = c("Vaginal", "C-section"))

# Check 
summary(combined_dataframe$delivery_mod)


###### Delivery mode (C-section initiation) ######
summary(combined_dataframe$STCESPARTO)

# Convert 9 to NA
combined_dataframe$emer_c_sec <- combined_dataframe$STCESPARTO
combined_dataframe$emer_c_sec[combined_dataframe$emer_c_sec == 9] <- NA

# label factor
combined_dataframe$emer_c_sec <- factor(combined_dataframe$emer_c_sec, levels = c("1", "2", "3"), labels = c("Pre-Labour initiation", "Post-Labour initiation", "Not applicable")) 

# check
summary(combined_dataframe$emer_c_sec)


###### Prenatal care initiation (access proxy) ######
summary(combined_dataframe$MESPRENAT)

# add 99 to NAs
combined_dataframe$bin_prenatal_care <- combined_dataframe$MESPRENAT
combined_dataframe$bin_prenatal_care[combined_dataframe$bin_prenatal_care == 99] <- NA

# binary (timely/adequate start of antenatal care - no: first antenatal appointment >1st trimester; yes: first antenatal appointment ≤1st trimester)
combined_dataframe$bin_prenatal_care <- factor(ifelse(as.numeric(as.character(combined_dataframe$bin_prenatal_care)) <= 3, 0, 1), 
                                               levels = c(0, 1), 
                                               labels = c("Access in 1st trimester", "Access after 1st trimester"))
#check
summary(combined_dataframe$bin_prenatal_care)


###### Maternal age ######
summary(combined_dataframe$IDADEMAE)

# add 99 to NAs
combined_dataframe$bin_mage<-combined_dataframe$IDADEMAE
combined_dataframe$bin_mage[combined_dataframe$bin_mage == 99] <- NA

# Create bins and assign labels
combined_dataframe$bin_mage <- factor(
  cut(as.numeric(as.character(combined_dataframe$bin_mage)), 
      breaks = c(-Inf, 19, 34, Inf), 
      labels = c("<20", "20-34", ">=35"), 
      right = TRUE))

#check
summary(combined_dataframe$bin_mage)


###### Maternal ethnicity ######
summary(combined_dataframe$RACACORMAE)

# add 9 to NAs
combined_dataframe$ethnicity <-combined_dataframe$RACACORMAE 
combined_dataframe$ethnicity[combined_dataframe$ethnicity  == 9] <- NA

# relabel
combined_dataframe$ethnicity <- factor(combined_dataframe$ethnicity, levels = c("1", "2", "3", "4", "5"), labels = c("White", "Black", "Asian", "Mixed", "Indigenous"))

# check
summary(combined_dataframe$ethnicity)


###### Maternal education ######
summary(combined_dataframe$ESCMAE)

# bin education
combined_dataframe$bin_edu <- factor(
  ifelse(as.numeric(as.character(combined_dataframe$ESCMAE)) %in% 1:4, 0,
         ifelse(as.numeric(as.character(combined_dataframe$ESCMAE)) == 5, 1, NA)), 
  levels = c(0, 1), 
  labels = c("Highschool or below", "Further education"))

# check
summary(combined_dataframe$bin_edu)


############################################################
# Attach municipality-level metrics (IBP, Köppen) #
############################################################

###### IBP ######
# check municip column character number
municipality_code_ibp_char <- as.character(df_ibp$codmun6)
all(nchar(municipality_code_ibp_char) == 6) # TRUE

# rename municip col names 
names(df_ibp)[names(df_ibp) == 'codmun6'] <- 'birth_municip'
names(df_ibp)[names(df_ibp) == 'q_measure_1f_12'] <- 'q_IBP'

# select only the required column from df_ibp
df_ibp_selected <- df_ibp[, c("birth_municip", "q_IBP")]

# Convert columns to factor
df_ibp_selected$birth_municip <- factor(df_ibp_selected$birth_municip)
df_ibp_selected$q_IBP <- factor(df_ibp_selected$q_IBP)

# Join datasets by municipality of birth
combined_dataframe <- merge(combined_dataframe, df_ibp_selected, by = "birth_municip", all.x = TRUE)

# examine 
summary(combined_dataframe$q_IBP)



###### Koppen ######
# check municip column character number
municipality_code_kop_char <- as.character(df_koppen$IBGE.Code6)
all(nchar(municipality_code_kop_char) == 6) # TRUE

# rename municip col name (6 character code)
names(df_koppen)[names(df_koppen) == 'IBGE.Code6'] <- 'birth_municip'

# select only the required columns 
df_kop_selected <- df_koppen[, c("birth_municip", "Koppen")]

# Convert columns to factor
df_kop_selected$birth_municip <- factor(df_kop_selected$birth_municip)
df_kop_selected$Koppen <- factor(df_kop_selected$Koppen)

# Join datasets by municipality of birth
combined_dataframe <- merge(combined_dataframe, df_kop_selected, by = "birth_municip", all.x = TRUE)

# examine
summary(combined_dataframe$Koppen)


############################################################
# Descriptive statistics #
############################################################

###### Tidy delivery mode ######
# Cross tabulation of delivery mode and emergency c-section
table(combined_dataframe$emer_c_sec, combined_dataframe$delivery_mod,useNA = "ifany")

# Create a new column for combined delivery mode 
combined_dataframe$combined_delivery_mod <- with(combined_dataframe, 
                                            ifelse(is.na(emer_c_sec) & delivery_mod == "Vaginal", "Vaginal", # missing in emer_c_sec and present in delivery_mod as vaginal, then vaginal
                                            ifelse(is.na(emer_c_sec) & delivery_mod == "C-section", "caesarean (Other)", # missing in emer_c_sec and present in delivery_mod as caesarean, then caesarean
                                            ifelse(is.na(delivery_mod), NA, # missing in delivery_mod then NA
                                            ifelse(delivery_mod == "Vaginal", "Vaginal", # if delivery_mod vaginal then always vaginal
                                            ifelse(delivery_mod == "C-section" & emer_c_sec == "Pre-Labour initiation", "caesarean (Pre-Labour initiation)", # if delivery_mod C-section and emer_c_sec is pre-labour, then pre-labour 
                                            ifelse(delivery_mod == "C-section" & emer_c_sec == "Post-Labour initiation", "caesarean (Post-Labour initiation)", # if delivery_mod C-section and emer_c_sec is post-labour, then post-labour
                                            ifelse(delivery_mod == "C-section" & emer_c_sec == "Not applicable", "caesarean (Other)", NA)))))))) # if delivery_mod C-section and emer_c_sec isnot applicable, then caesarean (other). Also if data fits no category then NA

table(as.factor(combined_dataframe$combined_delivery_mod), combined_dataframe$bin_apgar,useNA = "ifany")


###### Tidy Koppen ######
# remove unused levels
print(levels(combined_dataframe$Koppen))
combined_dataframe$Koppen <- droplevels(combined_dataframe$Koppen)


###### Function to loop through all variables ######
# Define a function to perform the analysis 
analyze_variables <- function(var1, var2, data) { #var1 will be the category, var2 will be bin_apgar
  # Count (with NAs)
  count_table <- table(data[[var1]], data[[var2]], useNA = "ifany")
  print("Count with NAs:")
  print(count_table)
  
  # Count 
  count_table_noNA <- table(data[[var1]], data[[var2]])
  print("Count without NAs:")
  print(count_table_noNA)
  
  # Percentage (excluding NAs)
  prop_table <- prop.table(count_table, margin = 2) * 100
  print("Percentage (excluding NAs):")
  print(round(prop_table, 2))  # Round to 2 decimal places
  
  # Chi-squared test (automatically excludes NAs)
  chi_sq_test <- chisq.test(count_table_noNA, correct=F)
  print("Chi-squared test:")
  print(chi_sq_test)}

# List of variables to analyze 
variables <- c("bin_mage", "bin_edu", "ethnicity", "bin_parity_2", "bin_prenatal_care", "combined_delivery_mod", "sex", "Koppen")  

# Loop through the list
for (var in variables) {
  analyze_variables(var, "bin_apgar", combined_dataframe)}


###### Descriptive analysis IBP (Chi2 cannot be done when there are categories of 5 and below observations) ######
# examine cross table as count and percentages
table(combined_dataframe$q_IBP, combined_dataframe$bin_apgar,useNA = "ifany")
round(prop.table(table(combined_dataframe$q_IBP, combined_dataframe$bin_apgar), margin = 2) * 100, 2)

# Create new 'q_IBP' column and replace '5' with NA
combined_dataframe$q_IBP_chi <- combined_dataframe$q_IBP
combined_dataframe$q_IBP_chi[combined_dataframe$q_IBP_chi == "5"] <- NA

# Drop unused levels 
combined_dataframe$q_IBP_chi <- droplevels(combined_dataframe$q_IBP_chi)

# Chi-squared test
chisq.test(table(combined_dataframe$q_IBP_chi, combined_dataframe$bin_apgar), correct=F)

# Remove 'q_IBP_chi' column
combined_dataframe$q_IBP_chi <- NULL


############################################################
# Save dataset (APGAR ≤7) #
############################################################

# filter for low Apgar births
low_apgar_df <- combined_dataframe[combined_dataframe$bin_apgar == "Low Apgar", ]

# Select specific columns
low_apgar_df <- low_apgar_df[, c("birth_municip", "dob", "clean_apgar", "bin_apgar","sex", "bin_parity_2", "delivery_mod",  "emer_c_sec", "combined_delivery_mod", "bin_prenatal_care", "bin_mage", "ethnicity", "bin_edu", "q_IBP", "Koppen"  )]

# reoder by date of birth
low_apgar_df <- low_apgar_df[order(low_apgar_df$dob), ]

# save
# write.csv(low_apgar_df, file = "GIT_health_dataset.csv", row.names = FALSE)


############################################################
# Format and save timeseries dataset #
############################################################

# Aggregate counts of low APGAR by municipality and date (count of low apgar births in each municipality on each day)
ts_df <- aggregate(bin_apgar ~ birth_municip + dob, 
                data = combined_dataframe[combined_dataframe$bin_apgar == "Low Apgar", ], 
                FUN = length)
names(ts_df)[3] <- "low_apgar_count"

# Aggregate total birth counts by municipality and date
birth_counts <- aggregate(bin_apgar ~ birth_municip + dob, 
                          data = combined_dataframe, 
                          FUN = length)
names(birth_counts)[3] <- "birth_count"

# Merge the two data frames (to get df with separate count of low apgar births and total births) - currently sorted by municpality code
result_df <- merge(ts_df, birth_counts, by = c("birth_municip", "dob"), all = TRUE)

#replace NA with 0 (i.e. if there were no low apgars on a given day)
result_df[is.na(result_df)] <- 0

# Extract date components
result_df$month <- format(result_df$dob, "%m")
result_df$month_year <- format(result_df$dob, "%Y-%m")
result_df$year <- format(result_df$dob, "%Y")

# How many unique birth municipalities
result_df$birth_municip <- droplevels(result_df$birth_municip)
nlevels(result_df$birth_municip) #448 

# Range of dates
range(result_df$dob)

# Count days with 0 births or 0 low Apgar births
sum(result_df$birth_count == 0) #no days with 0 births
sum(result_df$low_apgar_count == 0) #387006 days with 0 low APGAR 

# are all days represented
all_dates <- seq.Date(from = as.Date("2013-01-01"), to = as.Date("2019-12-31"), by = "day")
setdiff(all_dates, result_df$dob) #yes

# append Koppen zones
merged_result_df <- merge(result_df, df_kop_selected, by = "birth_municip", all.x = TRUE)

# order by date and municipality code
merged_result_df <- merged_result_df[order( merged_result_df$birth_municip, merged_result_df$dob), ]

# write csv
#write.csv(merged_result_df, "GIT_ts_dataset.csv", row.names = FALSE)

