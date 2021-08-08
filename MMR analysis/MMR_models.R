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
library(viridis)
library(ggpubr)
library(openxlsx)

# Load data that was created with Combine_MMR.R.
load("All.MMR.data")
head(All.MMR.data)

All.MMR.data$MR.abs # 2 NA, remove them
All.MMR.data <- filter(All.MMR.data, is.na(MR.abs) == FALSE)

## Edit genotypes E and L to EE or LL
All.MMR.data $Call_vgll3 <- dplyr::recode(All.MMR.data $Call_vgll3, E = "EE", L = "LL")
All.MMR.data $Call_six6 <- dplyr::recode(All.MMR.data $Call_six6, E = "EE", L = "LL")

## Removed 1 outlier earlier based on model
All.MMR.data[147,]
All.MMR.data <- filter(All.MMR.data, !(Ind =="A00000D0900001897007276"))
save(All.MMR.data, file ="Data/All.MMR.data")

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


## Model with only main effects (for safekeeping)

MMR.main.mod <- lmer(na.action = "na.fail" ,
               log10(MR.abs) ~
                 Treatment +
                 Sex +
                 Call_vgll3 +
                 Call_six6 +
                 log10(Mass)+
                 (1|Family) +
                 (1|Chamber.No),
              REML = F,
              data = All.MMR.data)

outlierTest(MMR.main.mod)

plot_model(MMR.main.mod, type = "diag")
MMR.main.mod.res <- cbind(All.MMR.data, resid(MMR.main.mod))
leveneTest(resid(MMR.main.mod)~ Treatment, data =  MMR.main.mod.res)
leveneTest(resid(MMR.main.mod)~ Sex, data =  MMR.main.mod.res)
leveneTest(resid(MMR.main.mod)~ Call_vgll3, data =  MMR.main.mod.res)
leveneTest(resid(MMR.main.mod)~ Call_six6, data =  MMR.main.mod.res)
summary(MMR.main.mod)
anova(MMR.main.mod)

#### Full model with interactions
MMR.full.mod<- lmer(na.action = "na.fail" ,
     log10(MR.abs) ~
       Treatment +
       Sex +
       Call_vgll3 +
       Call_six6 +
       Call_vgll3:Call_six6+
       log10(Mass)+
       Treatment:log10(Mass)+
       Treatment:Call_vgll3 +
       Treatment:Call_six6 +
       Sex : Call_vgll3+
       Sex: Call_six6 +
       (1|Family) +
       (1|Chamber.No),
     REML = F,
     data = All.MMR.data)

anova(MMR.full.mod)

confs.full <- confint(MMR.full.mod, oldNames=FALSE)

sum.MMR.full.mod<- summary(MMR.full.mod)
aov.MMR.full.mod<-anova(MMR.full.mod)

## Make a summary table with estimates and type 3 effects
table.MMR.full <- data.frame( Coefficient = row.names(sum.MMR.full.mod$coefficients),
                              Estimate = sum.MMR.full.mod$coefficients[,1],
                              SE = sum.MMR.full.mod$coefficients[,2],
                              SSq =  c(NA, aov.MMR.full.mod$`Sum Sq`),
                              `Den DF` = c(sum.MMR.full.mod$coefficients[,3][1], aov.MMR.full.mod$DenDF),
                              F = c(sum.MMR.full.mod$coefficients[,4][1], aov.MMR.full.mod$`F value`),
                              p = c(sum.MMR.full.mod$coefficients[,5][1], aov.MMR.full.mod$`Pr(>F)`))

# too many digits, round the numbers
table.MMR.full <- mutate(table.MMR.full, across(-1, round, 3))

table.MMR.full

# Table of random effects:
re.table.MMR.full <- as.data.frame(sum.MMR.full.mod$varcor) %>%
  dplyr::select(-(c(var1, var2, sdcor))) %>% #drop unnecessary cols
  rename("Random effect" = "grp", "Var" = "vcov") %>%
  mutate("C.I.low" = c(confs.full[1,1], confs.full[2,1], NA),
         "C.I.high" = c(confs.full[1,2], confs.full[2,2], NA),) %>%
  mutate(across(2:4, round, 4))

re.table.MMR.full

# export both into xlsx

write.xlsx(table.MMR.full, file = "table.MMR.full.xlsx")
write.xlsx(re.table.MMR.full, file = "re.table.MMR.full.xlsx")


### Simplify model, remove Treatment:Call_six6 
MMR.mod.2<- lmer(na.action = "na.fail" ,
                    log10(MR.abs) ~
                      Treatment +
                      Sex +
                      Call_vgll3 +
                      Call_six6 +
                      Call_vgll3:Call_six6+
                      log10(Mass)+
                      Treatment:log10(Mass)+
                      Treatment:Call_vgll3 +
                      Sex : Call_vgll3+
                      Sex: Call_six6 +
                      (1|Family) +
                      (1|Chamber.No),
                    REML = F,
                    data = All.MMR.data)

anova(MMR.mod.2)
confint(MMR.mod.2)

### Simplify model, remove Sex:Call_six6 
MMR.mod.2<- lmer(na.action = "na.fail" ,
                 log10(MR.abs) ~
                   Treatment +
                   Sex +
                   Call_vgll3 +
                   Call_six6 +
                   Call_vgll3:Call_six6+
                   log10(Mass)+
                   Treatment:log10(Mass)+
                   Treatment:Call_vgll3 +
                   Sex : Call_vgll3+
                   (1|Family) +
                   (1|Chamber.No),
                 REML = F,
                 data = All.MMR.data)

anova(MMR.mod.2)
confint(MMR.mod.2)


### Simplify model, remove Treatment:Call_vgll3 
MMR.mod.2<- lmer(na.action = "na.fail" ,
                 log10(MR.abs) ~
                   Treatment +
                   Sex +
                   Call_vgll3 +
                   Call_six6 +
                   Call_vgll3:Call_six6+
                   log10(Mass)+
                   Treatment:log10(Mass)+
                   Sex : Call_vgll3+
                   Sex: Call_six6 +
                   (1|Family) +
                   (1|Chamber.No),
                 REML = F,
                 data = All.MMR.data)

anova(MMR.mod.2)
confint(MMR.mod.2)

### Simplify model, remove Sex:Call_six6
MMR.mod.2<- lmer(na.action = "na.fail" ,
                 log10(MR.abs) ~
                   Treatment +
                   Sex +
                   Call_vgll3 +
                   Call_six6 +
                   Call_vgll3:Call_six6+
                   log10(Mass)+
                   Treatment:log10(Mass)+
                   Sex : Call_vgll3+
                   (1|Family) +
                   (1|Chamber.No),
                 REML = F,
                 data = All.MMR.data)

anova(MMR.mod.2)
summary(MMR.mod.2)
confint(MMR.mod.2)

### Simplify model, remove Sex:Call_vgll3 (similar p-value es Treatment-mass but much smaller effect)
MMR.mod.2<- lmer(na.action = "na.fail" ,
                 log10(MR.abs) ~
                   Treatment +
                   Sex +
                   Call_vgll3 +
                   Call_six6 +
                   Call_vgll3:Call_six6+
                   log10(Mass)+
                   Treatment:log10(Mass)+
                   (1|Family) +
                   (1|Chamber.No),
                 REML = F,
                 data = All.MMR.data)

anova(MMR.mod.2)
summary(MMR.mod.2)
confint(MMR.mod.2)

## Remove also Treatment:log10(Mass)

MMR.mod.2<- lmer(na.action = "na.fail" ,
                 log10(MR.abs) ~
                   Treatment +
                   Sex +
                   Call_vgll3 +
                   Call_six6 +
                   Call_vgll3:Call_six6+
                   log10(Mass)+
                   (1|Family) +
                   (1|Chamber.No),
                 REML = F,
                 data = All.MMR.data)

anova(MMR.mod.2)
anova(MMR.mod.2, type = 2) # are type 2 ssq consistent?
summary(MMR.mod.2)

confs.mod2 <- confint(MMR.mod.2, oldNames=FALSE)

# Final model, significant gene interaction
sum.MMR.mod.2 <- summary(MMR.mod.2)
aov.MMR.mod.2 <- anova(MMR.mod.2)

## Make a summary table with estimates and type 3 effects (report this model)
table.MMR.mod2 <- data.frame( Coefficient = row.names(sum.MMR.mod.2$coefficients),
                              Estimate = sum.MMR.mod.2$coefficients[,1],
                              SE = sum.MMR.mod.2$coefficients[,2],
                              SSq =  c(NA, aov.MMR.mod.2$`Sum Sq`),
                              `Den DF` = c(sum.MMR.mod.2$coefficients[,3][1], aov.MMR.mod.2$DenDF),
                              F = c(sum.MMR.mod.2$coefficients[,4][1], aov.MMR.mod.2$`F value`),
                              p = c(sum.MMR.mod.2$coefficients[,5][1], aov.MMR.mod.2$`Pr(>F)`))

table.MMR.mod2 <- mutate(table.MMR.mod2, across(-1, round, 3))

table.MMR.mod2

# Table of random effects:
re.table.MMR.mod2 <- as.data.frame(sum.MMR.full.mod$varcor) %>%
  dplyr::select(-(c(var1, var2, sdcor))) %>% #drop unnecessary cols
  rename("Random effect" = "grp", "Var" = "vcov") %>%
  mutate("C.I.low" = c(confs.mod2[1,1], confs.mod2[2,1], NA),
         "C.I.high" = c(confs.mod2[1,2], confs.mod2[2,2], NA)) %>%
  mutate(across(2:4, round, 4))

re.table.MMR.mod2

# export both into xlsx

write.xlsx(table.MMR.mod2, file = "table.MMR.final.xlsx")
write.xlsx(re.table.MMR.mod2, file = "re.table.MMR.final.xlsx")

##################################################
## Pairwise comparison of genotype effects:
##################################################
prwise <- emmeans(MMR.mod.2,  ~ Call_vgll3 | Call_six6)
pairs(prwise) # vgll3 diff

prwise <- emmeans(MMR.mod.2,  ~ Call_six6 | Call_vgll3)
pairs(prwise) # six6 diff

##################################################
## Predicted means
##################################################

MMR.predict <- ggpredict(MMR.mod.2, terms = c("Call_vgll3","Call_six6"), type = "fe", ci.lvl = 0.9,
                            back.transform = FALSE)
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
