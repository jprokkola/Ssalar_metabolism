######################################################################
## Statistical models of MMR in high and low food groups
## Spline-MMR data
######################################################################
library(lme4)
library(lmerTest)
library(lmfor)
library(car)
library(ggeffects)
library(tidyverse)
library(merTools)
library(emmeans)
library(sjPlot)
library(ggplot2)
library(openxlsx)
library(lubridate)
library(partR2)
library(MuMIn)

# Load data that was created with Combine_MMR.R.
load("Data/All.MMR.data")
head(All.MMR.data)

All.MMR.data$MR.abs # 2 NA, remove them
All.MMR.data <- filter(All.MMR.data, is.na(MR.abs) == FALSE)

## Edit genotypes E and L to EE or LL (for plots)
All.MMR.data $Vgll3 <- dplyr::recode(All.MMR.data $Call_vgll3, E = "EE", L = "LL")
All.MMR.data $Six6 <- dplyr::recode(All.MMR.data $Call_six6, E = "EE", L = "LL")

## Edit genotypes E and L to -0.5 or 0.5
All.MMR.data $Call_vgll3 <- dplyr::recode(All.MMR.data $Call_vgll3, E = -0.5, L = 0.5)
All.MMR.data $Call_six6 <- dplyr::recode(All.MMR.data $Call_six6, E = -0.5, L = 0.5)

## Edit sex F and M to -0.5 or 0.5
All.MMR.data $Sex <- dplyr::recode(All.MMR.data $Sex, "F" = -0.5, M = 0.5)

## Edit treatment High food and Low food to -0.5 or 0.5
All.MMR.data $Treatment_centred <- dplyr::recode(All.MMR.data $Treatment, "High food" = -0.5, "Low food" = 0.5)

## Removed 1 outlier earlier based on model
All.MMR.data[147,]
All.MMR.data <- filter(All.MMR.data, !(Ind =="A00000D0900001897007276"))
save(All.MMR.data, file ="../Data/All.MMR.data")

#write.table(All.MMR.data, file= "Data/Analysed.MMR.data.txt", row.names = F,
#            quote = F, dec = ".", sep = "\t")

## Plot MMR vs mass
ggplot(All.MMR.data, aes(y= log10(MR.abs), x=log10(Mass)))+
  geom_point(alpha=0.5) +
  stat_smooth(method = "lm", level = 0.95)+
  ylab(expression(paste("Log"[10], " MMR (mg ", O[2], "/h)"))) +
  xlab(expression(paste("Log"[10], " body mass (g)"))) +
  theme_bw()

# Scaling exponents
All.MMR.data %>%
  summarise(b = coef(lm(log10(MR.abs)~log10(Mass)))[2],
            R2 = summary(lm(log10(MR.abs)~log10(Mass)))$r.squared)


#########################################
### Full model with interactions
#########################################

MMR.full.mod<- lmer(na.action = "na.fail" ,
                    log10(MR.abs) ~
                      Treatment_centred +
                      Sex +
                      Call_vgll3 +
                      Call_six6 +
                      scale(log10(Mass))+
                      scale(Order) + #the order of testing pairs of fish 1-8
                      Treatment_centred:scale(log10(Mass))+
                      Treatment_centred:Call_vgll3 +
                      Treatment_centred:Call_six6 +
                      Sex:Call_vgll3 +
                      Sex:Call_six6 +
                      Call_vgll3:Call_six6 +
                      (1|Chamber.No) +
                      (1|initial) +
                      (1|Family),
                    REML = F,
                    data = All.MMR.data)

summary(MMR.full.mod)
outlierTest(MMR.full.mod)
## plot diagnostics
plot_model(MMR.full.mod, type = "diag")

# Residuals and fitted values with lines to help visualization:
MMR.full.mod.res <- cbind(All.MMR.data, resid(MMR.full.mod))
plot(fitted(MMR.full.mod),resid(MMR.full.mod, type="pearson"))
mywhiskers(fitted(MMR.full.mod),resid(MMR.full.mod, type="pearson"), add=T, se=F)

# heteroscedasticity tests
MMR.main.mod.res <- cbind(All.MMR.data, resid(MMR.full.mod))
leveneTest(resid(MMR.full.mod)~ Treatment, data =  MMR.full.mod.res)
leveneTest(resid(MMR.full.mod)~ as.factor(Sex), data =  MMR.full.mod.res)
leveneTest(resid(MMR.full.mod)~ Vgll3, data =  MMR.full.mod.res)
leveneTest(resid(MMR.full.mod)~ Six6, data =  MMR.full.mod.res)

### Export estimates and Type III results
sum.MMR.full.mod<- summary(MMR.full.mod)
aov.MMR.full.mod<-anova(MMR.full.mod)
confs.full <- confint(MMR.full.mod, oldNames = FALSE)

## Make a summary table with estimates and type 3 effects
table.MMR.full <- data.frame( Coefficient = row.names(sum.MMR.full.mod$coefficients),
                              Estimate = sum.MMR.full.mod$coefficients[,1],
                              SE = sum.MMR.full.mod$coefficients[,2],
                              SSq =  c(NA, aov.MMR.full.mod$`Sum Sq`),
                              `Den DF` = c(sum.MMR.full.mod$coefficients[,3][1], aov.MMR.full.mod$DenDF),
                              F = c(sum.MMR.full.mod$coefficients[,4][1], aov.MMR.full.mod$`F value`),
                              p = c(sum.MMR.full.mod$coefficients[,5][1], aov.MMR.full.mod$`Pr(>F)`))

# too many digits, round the numbers
table.MMR.full <- mutate(table.MMR.full, across(-1, round, 4))

table.MMR.full

# Table of random effects:
re.table.MMR.full <- as.data.frame(sum.MMR.full.mod$varcor) %>%
  dplyr::select(-(c(var1, var2, sdcor))) %>% #drop unnecessary cols
  rename("Random effect" = "grp", "Var" = "vcov") %>%
  mutate("C.I.low" = c(confs.full[1,1]^2, confs.full[2,1]^2,confs.full[3,1]^2, NA),
         "C.I.high" = c(confs.full[1,2]^2, confs.full[2,2]^2, confs.full[3,2]^2,NA),) %>%
  mutate(across(2:4, round, 4))

re.table.MMR.full

# export both to xlsx
write.xlsx(table.MMR.full, file = "table.MMR.full.xlsx")
write.xlsx(re.table.MMR.full, file = "re.table.MMR.full.xlsx")

##################################################
## Get a more parsimonious model with MuMIn
##################################################
all.mumin.MMR <-  dredge(MMR.full.mod,rank = "AICc")

subset.mumin.MMR <- subset(all.mumin.MMR, delta <2)

# save as table
subset.mumin.MMR.table <- as.data.frame(subset.mumin.MMR)
subset.mumin.MMR.table <- subset.mumin.MMR.table %>%
  mutate(across(-"df", round, 3))

write.xlsx(subset.mumin.MMR.table, file = "subset.mumin.MMR.table.xlsx")

## Calculate average model
ave.model <- model.avg(subset.mumin.MMR, fit =T)
names(ave.model)
summary(ave.model)
#save(ave.model, file = "ave_model_MMR")

# Take averages using the full model method from the subset of models
ave.model.table <- as.data.frame(summary(ave.model)$coefmat.full)

# Add importance variable
importance <- data.frame("importance" = ave.model$sw)
ave.model.table <- inner_join ( rownames_to_column(ave.model.table), rownames_to_column(importance),
                                by = "rowname")

# round numbers
ave.model.table <- ave.model.table %>%
  mutate(across(3:7, round, 3),
         across(2, round, 4)) %>%
  #rename and organize
  rename(Coef = rowname, SE = `Std. Error`, z = `z value`, P = `Pr(>|z|)`) %>%
  dplyr::select(-4)  %>%
  relocate(importance, .before = P)

ave.model.table

# export to xlsx
write.xlsx(ave.model.table, file = "ave.model.table.MMR.xlsx")

##################################################
## Check variance partitioning from genotypes
##################################################
## need to run a genotype model with mass-corrected data (i.e. the relevant response)

All.MMR.data$resid_MMR <- resid(lm(log10(MR.abs) ~ log10(Mass), data =All.MMR.data))

MMR.mod.part <- lmer(na.action = "na.fail" ,
                  resid_MMR ~
                   Call_vgll3 +
                   Call_six6 +
                   Call_vgll3:Call_six6+
                   (1|Chamber.No) +
                   (1|initial) +
                   (1|Family),
                 REML = F,
                 data = All.MMR.data)

partR2_MMR <- partR2(MMR.mod.part, partvars = c("Call_vgll3", "Call_six6", "Call_vgll3:Call_six6"), nboot = 1000)

R2s <- partR2_MMR$R2
R2s <- mutate(R2s, across(-1, round, 5))
              
write.xlsx(R2s, file ="MMR_R2s.xlsx")

##################################################
## Pairwise comparison of genotype effects:
##################################################

## Use the top model from dredge, because average model not compatible with emmeans
prwise <- emmeans(MMR.full.mod,  ~ Call_vgll3 | Call_six6)
pairs(prwise) # vgll3 diff

prwise <- emmeans(MMR.full.mod,  ~ Call_six6 | Call_vgll3)
pairs(prwise) # six6 diff

##################################################
## Predicted means
##################################################

MMR.predict <- ggpredict(MMR.full.mod, terms = c("Call_vgll3","Call_six6"), type = "fe", ci.lvl = 0.9,
                          back.transform = FALSE, condition =c("Mass" = 3.7))
MMR.predict
plot(MMR.predict)

# For presenting results, calculate MMR per kg 
# Manually back-transform the values to mg/O2 (per 3.71g)
MMR.predict[,2:5]<- 10^ MMR.predict[,2:5]
# Calculate in mg/O2/kg
MMR.predict[,2:5]<- MMR.predict[,2:5]/0.00371

## Plot made separately after combining SMR, MMR, and AS. Save prediction table.
save(MMR.predict, file = "MMR.predict")

##################################################
## Calculating Ns
##################################################
all_mmr_N <- All.MMR.data %>%
  group_by(Treatment, Family, Call_vgll3, Call_six6) %>%
  summarise(N = length(MR.abs[is.na(MR.abs)=="FALSE"]))

