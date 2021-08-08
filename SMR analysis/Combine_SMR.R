#################################################
### 1. Combine High food and Low food SMR data 
### 2. Compare different SMR metrics
### 3. Get N per group
#################################################
library(lubridate)
library(tidyverse)
library(ggplot2)
library(ggcorrplot)
library(ggpubr)

# Load R objects from high food and low food groups
load("Data/HF.SMR.data")
head(HF.SMR)
str(HF.SMR) 

load("Data/LF.SMR.data")
head(LF.SMR)
str(LF.SMR) # One batch excluded previously

# Add treatments and bind rows
LF.SMR$Treatment <- "Low food"
HF.SMR$Treatment <- "High food"

All.SMR.data <- rbind(LF.SMR, HF.SMR)

str(All.SMR.data)
save(All.SMR.data, file = "Data/All.SMR.data")

### Compare SMR from q0.1, q0.2 and MLND approaches and choose one
# All will correlate with mass, so to evaluate mass-independent data, 
# get residuals from mass-SMR regression for each first

# mass-adjusted values
All.SMR.data<- All.SMR.data %>%
  mutate(
    q0.1.SMR.resid  = resid(lm(log10(SMR.abs.q0.1) ~ log10(Mass))),
    q0.2.SMR.resid  = resid(lm(log10(SMR.abs.q0.2) ~ log10(Mass))),
    mlndSMR.resid  = resid(lm(log10(SMR.abs.mlnd) ~ log10(Mass)))
  )

# Plots of correlations
ggplot(All.SMR.data, aes(y= mlndSMR.resid, x=q0.1.SMR.resid))+
  facet_wrap(~Treatment) +
  geom_point(alpha=0.5) +
  stat_smooth(method = "lm", level = 0.95)+
  theme_bw()

ggplot(All.SMR.data, aes(y= mlndSMR.resid, x=q0.2.SMR.resid))+
  facet_wrap(~Treatment) +
  geom_point(alpha=0.5) +
  stat_smooth(method = "lm", level = 0.95)+
  theme_bw()

ggplot(All.SMR.data, aes(y= q0.1.SMR.resid, x=q0.2.SMR.resid))+
  facet_wrap(~Treatment) +
  geom_point(alpha=0.5) +
  stat_smooth(method = "lm", level = 0.95)+
  theme_bw()

# Use ggcoorplot to get all pairwise correlations and p-values
corrdata <- All.SMR.data[,c("q0.1.SMR.resid", "q0.2.SMR.resid", "mlndSMR.resid")]
corr <- round(cor(corrdata), 2)
# use SMR-MLND for analysis

### N in each group
# Summarise N in each family-genotype combination
N_summary <- All.SMR.data %>%
  group_by(Treatment, Family, Call_vgll3, Call_six6) %>%
  summarise(N = length(PIT))

# N by treatment
All.SMR.data%>%
  group_by(Treatment) %>%
  summarise(N = length(PIT))

