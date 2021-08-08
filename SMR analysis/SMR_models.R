##############################################
### SMR statistical models, High food and Low food data
### Input data from Combine_SMR.R
#############################################
library(lme4)
library(lmerTest)
library(car)
library(ggeffects)
library(tidyverse)
library(merTools)
library(emmeans)
library(sjPlot)
library(ggplot2)
library(ggpubr)
library(viridis)
library(openxlsx)

#############################
## Load and plot data
############################
### SMR
load("Data/All.SMR.data")
str(All.SMR.data)

## Replace genotypes E and L with EE or LL
All.SMR.data$Call_vgll3 <- dplyr::recode(All.SMR.data$Call_vgll3, E = "EE", L = "LL")
All.SMR.data$Call_six6 <- dplyr::recode(All.SMR.data$Call_six6, E = "EE", L = "LL")

#Box plot
ggplot(All.SMR.data, aes(y= SMR.abs.mlnd, x=Family, color = Treatment))+
  geom_boxplot() +
  geom_jitter(alpha =0.5)+
  theme_bw()+
  theme(panel.spacing.x=unit(0, "lines"),
        axis.title.x = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())

#histogram
hist(All.SMR.data$SMR.abs.mlnd)

#########################################
### Models
#########################################
# Main effect model
SMR.main.mod <-lmer(log10(SMR.abs.mlnd) ~   
                      Treatment +
                      Call_vgll3 + 
                      Call_six6 + 
                      Sex +
                      log10(Mass)+
                      (1|Family) +
                      (1|Batch),
                    REML = F, 
                    data= All.SMR.data,
                    na.action = na.fail
)

summary(SMR.main.mod)
anova(SMR.main.mod)
outlierTest(SMR.main.mod)
# Strong outlier:
All.SMR.data[67,]
# Filter this fish from the data
All.SMR.data<- filter(All.SMR.data,
                   !(PIT == "A00000D0900001897006087"))

# rerun model
save(All.SMR.data, file= "Data/All.SMR.data")

# Full model with interactions
SMR.full.mod <-lmer(log10(SMR.abs.mlnd) ~  
                  Treatment +
                  Call_vgll3 + 
                  Call_six6 + 
                  Sex +
                  log10(Mass)+
                  log10(Mass):Treatment +
                  Call_vgll3:Treatment+
                  Call_six6:Treatment+
                  Call_vgll3:Sex+
                  Call_six6:Sex+
                  Call_vgll3:Call_six6+
                  (1|Family) +
                  (1|Batch),
                REML = F, 
                data= All.SMR.data,
                na.action = na.fail
)

summary(SMR.full.mod)
anova(SMR.full.mod)
confs.full <- confint(SMR.full.mod, oldNames = FALSE)

sum.SMR.full.mod <- summary(SMR.full.mod)
aov.SMR.full.mod <- anova(SMR.full.mod)

## Make a summary table with estimates and type 3 effects
table.SMR.full.mod <- data.frame(Coefficient = row.names(sum.SMR.full.mod$coefficients),
                                Estimate = sum.SMR.full.mod$coefficients[,1],
                                SE = sum.SMR.full.mod$coefficients[,2],
                                SSq =  c(NA, aov.SMR.full.mod$`Sum Sq`),
                                `Den DF` = c(sum.SMR.full.mod $coefficients[,3][1], aov.SMR.full.mod$DenDF),
                                F = c(sum.SMR.full.mod $coefficients[,4][1], aov.SMR.full.mod$`F value`),
                                p = c(sum.SMR.full.mod $coefficients[,5][1], aov.SMR.full.mod$`Pr(>F)`))

table.SMR.full.mod <- mutate(table.SMR.full.mod, across(-1, round, 3))
table.SMR.full.mod

# Table of random effects:
re.table.SMR.full <- as.data.frame(sum.SMR.full.mod$varcor) %>%
  dplyr::select(-(c(var1, var2, sdcor))) %>% #drop unnecessary cols
  rename("Random effect" = "grp", "Var" = "vcov") %>%
  mutate("C.I.low" = c(confs.full[1,1], confs.full[2,1], NA),
         "C.I.high" = c(confs.full[1,2], confs.full[2,2], NA),) %>%
  mutate(across(2:4, round, 4))

re.table.SMR.full

# export both into xlsx

write.xlsx(table.SMR.full.mod, file = "table.SMR.full.mod.xlsx")
write.xlsx(re.table.SMR.full, file = "re.table.SMR.full.xlsx")


# Simplify model, remove Call_vgll3:Sex
SMR.mod.2<-lmer(log10(SMR.abs.mlnd) ~  
                       Treatment +
                       Call_vgll3 + 
                       Call_six6 + 
                       Sex +
                       log10(Mass)+
                       log10(Mass):Treatment +
                       Call_vgll3:Treatment+
                       Call_six6:Treatment+
                       Call_six6:Sex+
                       Call_vgll3:Call_six6+
                       (1|Family) +
                       (1|Batch),
                     REML = F, 
                     data= All.SMR.data,
                     na.action = na.fail
)

summary(SMR.mod.2)
anova(SMR.mod.2)
confint(SMR.mod.2)

# Simplify model, remove Call_six6:sex
SMR.mod.2<-lmer( log10(SMR.abs.mlnd) ~  
                       Treatment +
                       Call_vgll3 + 
                       Call_six6 + 
                       Sex +
                       log10(Mass)+
                       log10(Mass):Treatment +
                       Call_vgll3:Treatment+
                       Call_six6:Treatment+
                       Call_vgll3:Call_six6+
                       (1|Family) +
                       (1|Batch),
                     REML = F, 
                     data= All.SMR.data,
                     na.action = na.fail
)

summary(SMR.mod.2)
anova(SMR.mod.2)
confint(SMR.mod.2)

# remove Treatment:Call_vgll3 
SMR.mod.2<-lmer( log10(SMR.abs.mlnd) ~  
                       Treatment +
                       Call_vgll3 + 
                       Call_six6 + 
                       Sex +
                       log10(Mass)+
                       log10(Mass):Treatment +
                       Call_six6:Treatment+
                       Call_vgll3:Call_six6+
                       (1|Family) +
                       (1|Batch),
                     REML = F, 
                     data= All.SMR.data,
                     na.action = na.fail
)

summary(SMR.mod.2)
anova(SMR.mod.2)
confint(SMR.mod.2)

# remove Call_vgll3:Call_six6
SMR.mod.2<-lmer( log10(SMR.abs.mlnd) ~  
                       Treatment +
                       Call_vgll3 + 
                       Call_six6 + 
                       Sex +
                       log10(Mass)+
                       log10(Mass):Treatment +
                       Call_six6:Treatment+
                       (1|Family) +
                       (1|Batch),
                     REML = F, 
                     data= All.SMR.data,
                     na.action = na.fail
)

summary(SMR.mod.2)
anova(SMR.mod.2)

confs <- confint(SMR.mod.2, oldNames=FALSE)
# Pairwise differences between six6 at food levels
emm.smr <- emmeans(SMR.mod.2, specs = c("Treatment", "Call_six6"))
pairs(emm.smr)
# non-significant

# Remove  Call_six6:Treatment
SMR.mod.2 <-lmer( log10(SMR.abs.mlnd) ~  
                       Treatment +
                       Sex +
                       Call_vgll3 + 
                       Call_six6 + 
                       log10(Mass)+
                       log10(Mass):Treatment +
                       (1|Family) +
                       (1|Batch),
                     REML = F, 
                     data= All.SMR.data,
                     na.action = na.fail
)

sum.SMR.mod.2 <- summary(SMR.mod.2)
aov.SMR.mod.2 <- anova(SMR.mod.2)
plot_model(SMR.mod.2, type = "diag")

SMR.mod.res <- cbind(All.SMR.data, resid(SMR.mod.2))
#Residuals and fitted values with lines to help visualization:
plot(fitted(SMR.mod.2),resid(SMR.mod.2, type="pearson"))
mywhiskers(fitted(SMR.mod.2),resid(SMR.mod.2, type="pearson"), add=T, se=F)

leveneTest(resid(SMR.mod.2)~ Treatment, data =  SMR.mod.res)
boxplot(SMR.mod.res$`resid(SMR.mod.2)`~SMR.mod.res$Treatment) # a little heteroscedastic, but sample
# sizes not equal, and effect not strong

leveneTest(resid(SMR.mod.2)~ Sex, data =  SMR.mod.res)
leveneTest(resid(SMR.mod.2)~ Call_vgll3, data =  SMR.mod.res)
leveneTest(resid(SMR.mod.2)~ Call_six6, data =  SMR.mod.res)

## Make a summary table with estimates and type 3 effects
table.SMR.mod.2 <- data.frame(Coefficient = row.names(sum.SMR.mod.2$coefficients),
                                Estimate = sum.SMR.mod.2$coefficients[,1],
                                SE = sum.SMR.mod.2$coefficients[,2],
                                SSq =  c(NA, aov.SMR.mod.2$`Sum Sq`),
                                `Den DF` = c(sum.SMR.mod.2$coefficients[,3][1], aov.SMR.mod.2$DenDF),
                                F = c(sum.SMR.mod.2$coefficients[,4][1], aov.SMR.mod.2$`F value`),
                                p = c(sum.SMR.mod.2$coefficients[,5][1], aov.SMR.mod.2$`Pr(>F)`))

table.SMR.mod.2 <- mutate(table.SMR.mod.2, across(-1, round, 3))

# Table of random effects with 95% confint:
re.table.SMR.mod.2 <- as.data.frame(sum.SMR.mod.2$varcor) %>%
  dplyr::select(-(c(var1, var2, sdcor))) %>% #drop unnecessary cols
  rename("Random effect" = "grp", "Var" = "vcov") %>%
  mutate("C.I.low" = c(confs[1,1], confs[2,1], NA),
         "C.I.high" = c(confs[1,2], confs[2,2], NA),) %>%
  mutate(across(2:4, round, 4))

re.table.SMR.mod.2

# export both into xlsx

write.xlsx(table.SMR.mod.2, file = "table.SMR.mod.2.xlsx")
write.xlsx(re.table.SMR.mod.2, file = "re.table.SMR.mod.2.xlsx")

#########################
## Scaling exponents
#########################
All.SMR.data %>%
  dplyr::group_by(Treatment) %>%
  summarise(b = coef(lm(log10(SMR.abs.mlnd)~log10(Mass)))[2],
            R2 = summary(lm(log10(SMR.abs.mlnd)~log10(Mass)))$r.squared)

##################################################
## SMR-mass scaling plot using partial residuals
################################################## 
interactions::interact_plot(SMR.mod.2, pred = Mass, modx = Treatment, interval = T, 
                            int.type = "confidence", plot.points = T, partial.residuals = T, 
                            point.alpha = 0.4, colors = c("#F57B17FF", "#FCAE12FF")) +
  coord_trans(x = "log") +
  scale_color_viridis(discrete=TRUE, option = "inferno", begin = 0.1, end = 0.7) +
  #scale_colour_manual(values = c("#F57B17FF", "grey"), breaks = waiver()) +
  ylab(expression(paste("Log"[10], " SMR (mg ", O[2], "/h)"))) +
  xlab(expression(paste("Body mass (g)"))) +
  theme(legend.position = c(0.8,0.2),
        legend.title = element_blank())

##################################################
## Predicted means
##################################################
SMR.predict <- ggpredict(SMR.mod.2, terms =  c("Call_vgll3", "Call_six6"),
                          ci.lvl = 0.9,
                          back.transform = FALSE)
SMR.predict
plot(SMR.predict)

# Manually back-transform the values to mg/O2 (per 3.68g)
SMR.predict[,2:5]<- 10^SMR.predict[,2:5]
# Calculate in mg/O2/kg
SMR.predict[,2:5]<- SMR.predict[,2:5]/0.00368

# save object, combine later
save(SMR.predict, file = "SMR.predict")


