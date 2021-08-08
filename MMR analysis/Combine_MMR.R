#################################################
### 1. Combine High food and Low food MMR spline data 
### 2. Get N per group
#################################################
library(lubridate)
library(tidyverse)
library(ggplot2)

# Load R objects from high food and low food groups
load("Data/MMR_HF")
head(MMR_HF)

load("Data/MMR_LF")
head(MMR_LF)

# Add treatments and bind rows
MMR_LF$Treatment <- "Low food"
MMR_HF$Treatment <- "High food"

All.MMR.data <- rbind(MMR_LF, MMR_HF)
str(All.MMR.data)
save(All.MMR.data, file = "Data/All.MMR.data")

### N in each group
# Summarise N in each family-genotype combination
N_summary <- All.MMR.data %>%
  group_by(Treatment, Family, Call_vgll3, Call_six6) %>%
  summarise(N = length(Ind))

