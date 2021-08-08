######################################################################
## Combine SMR and MMR for aerobic scope analysis 
## SMR data: MLND method, MMR data: spline method
## Get mass-corrected values for correlations, and N per group
######################################################################

library(tidyverse)
library(data.table)
library(xlsx)

## Load SMR and MMR data with both High food and Low food treatments, combine, calculate absolute aerobic scope.

### SMR
load("Data/All.SMR.data")
str(All.SMR.data)

## exclude columns and rows that are not used
SMR.data <- All.SMR.data %>%
  dplyr::select(-c(tank, Combined, Chamber.No, SMR.mass.q0.1,
            SMR.abs.q0.1, SMR.mass.q0.2, SMR.abs.q0.2,
            Temp)) %>%
  droplevels()

str(SMR.data)

### MMR
load("Data/All.MMR.data")
str(All.MMR.data)

## exclude columns that are not used
MMR.data <- All.MMR.data %>%
  dplyr::select(c(Ind, MR.abs, MR.mass, Chamber.No)) %>%
  rename("MMR.Chamber" = Chamber.No)

# Combine together when both data are found for the same fish.
all.AS.data <- inner_join(SMR.data, MMR.data, 
                          by = c("PIT" = "Ind")) %>%
  droplevels()

str(all.AS.data)

## Replace genotypes E and L with EE or LL
all.AS.data$Call_vgll3 <- dplyr::recode(all.AS.data$Call_vgll3, E = "EE", L = "LL")
all.AS.data$Call_six6 <- dplyr::recode(all.AS.data$Call_six6, E = "EE", L = "LL")

# Calculate absolute AS, and residual traits adjusted for family and mass (used for calculating correlations). 

all.AS.data <- all.AS.data %>%
  mutate(absAS = MR.abs - SMR.abs.mlnd) %>%
  # filter missing data
  filter(is.na(absAS) == FALSE) %>%
  mutate(
    all.SMR.resid  = resid(lmer(log10(SMR.abs.mlnd) ~ log10(Mass)  + (1|Family))),
    all.MMR.resid  = resid(lmer(log10(MR.abs) ~ log10(Mass) + (1|Family))),
    absAS.resid  = resid(lmer(log10(absAS) ~ log10(Mass) + (1|Family))),
  )

head(all.AS.data)

save(all.AS.data, file ="Data/all.AS.data") # This is input data for models.

### Get N per sex and genotype per family and treatment 
ASfishcount <- all.AS.data %>%
  mutate("vgll3_six6" = paste(Call_vgll3, Call_six6))%>%
  group_by(Treatment, vgll3_six6, Family) %>%
  count(Sex) %>%
  pivot_wider(values_from = n, names_from= Sex) %>%
  as.data.frame()

