
setwd("~/UNI/NONPARAMETRIC STATISTICS/Project/RefugeeElectricity")
library(readxl)

baseline <- data.frame(read_excel("Data/Refugee_Settlements_Electricity_Access_DB_ver02-0.xlsx", sheet = "Baseline"))
head(baseline)

tier2 <- data.frame(read_excel("Data/Refugee_Settlements_Electricity_Access_DB_ver02-0.xlsx", sheet = "Tier 2"))
head(tier2)

tier3 <- data.frame(read_excel("Data/Refugee_Settlements_Electricity_Access_DB_ver02-0.xlsx", sheet = "Tier 3"))
head(tier3)

