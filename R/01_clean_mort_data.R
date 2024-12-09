# Loading necessary packages
library(lubridate)
library(here)
library(data.table)
library(tidyverse)

# Do not rely on this to completely clean your environment
# Better to do a full restart of R before running
rm(list=ls())

# Set working directory
i_am('R/01_clean_mort_data.R')

in_dir <- here('data')
out_dir <- here('output')

# We define a function to read the CDC WONDER data. We read the different
# files and combine then by appending rows. We also set missing values for 
# Suppressed, Missing, and Non Applicable values.
read_CDC_data <- function(file) {
  
  data <- fread(
    file,
    na.strings = c("Missing", "Suppressed", "Not Applicable"),
    keepLeadingZeros = TRUE,
    colClasses = c("character")
  )
  
  data
}

# Import death counts and mid-year population estimates by 5-year age groups
# sex, and urbanization codes downloaded from CDC WONDER
# https://wonder.cdc.gov/ucd-icd10.html
# NOTE: Cells with < 10 deaths are censored in these data
mort_data <- read_CDC_data(here(in_dir,'CDCMortByMetro','mortByMetro5YrAgeGroups.txt'))

# Setting intuitive names
mort_data <- mort_data[
  ,
  .(urbanization=`2013 Urbanization`,
    age_group=`Five-Year Age Groups Code`,
    sex=Gender,
    year=Year,
    nDx=Deaths,
    nEx=Population)
]

# Drop data for older age groups. Data on 90+ has to come from another file
# with 10-year age groups because population estimates are not available
# for 5-year groups above age 85.
mort_data <- mort_data[!(age_group %in% c('85-89','90-94','95-99','100+','NS'))]

# Set cols to numeric
numeric_cols <- c('year','nDx','nEx')
mort_data[, (numeric_cols) := lapply(.SD, as.integer), .SDcols = numeric_cols]

# Import death counts and mid-year population estimates by 10-year age groups
# sex, and urbanization codes downloaded from CDC WONDER
# https://wonder.cdc.gov/ucd-icd10.html
# NOTE: Cells with < 10 deaths are censored in these data
mort_data_10 <- read_CDC_data(here(in_dir,'CDCMortByMetro','mortByMetro10YrAgeGroups.txt'))

# Setting intuitive names
mort_data_10 <- mort_data_10[
  ,
  .(urbanization=`2013 Urbanization`,
    age_group=`Ten-Year Age Groups Code`,
    sex=Gender,
    year=Year,
    nDx=Deaths,
    nEx=Population)
]

# Keep only data for 85+ which is what we are missing in the 5-year age groups
# file
mort_data_10 <- mort_data_10[age_group == '85+']
    
# Set cols to numeric
numeric_cols <- c('year','nDx','nEx')
mort_data_10[, (numeric_cols) := lapply(.SD, as.integer), .SDcols = numeric_cols]

# Combine the two datasets so we have a complete set of counts for 5-year age
# groups and 85+ as the open-ended interval.
mort_data <- rbindlist(list(mort_data,mort_data_10)) 

# Create variable for interval starting age. We call it x to conform with 
# standard demographic notation even though it's not exactly a good programming
# practice
mort_data[,x:=as.integer(str_match(age_group,'\\d+'))]
mort_data[,x:=if_else(age_group=='1',0,x)]
setcolorder(mort_data,c('year','urbanization','sex','age_group','x','nDx','nEx'))

# Sort the data to make it easier to see
setorder(mort_data,urbanization,sex,year,x)

# Save data so that we do not have to repeat these operations each time
fwrite(mort_data,here(out_dir,'mort_data_clean.csv'))
