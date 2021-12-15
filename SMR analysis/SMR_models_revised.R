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
library(openxlsx)
library(interactions)
library(lmfor)
library(MuMIn)

#############################
## Load and plot data
############################
### SMR
load("Data/All.SMR.data")
str(All.SMR.data)

## Edit genotypes E and L to EE or LL (for plots)
All.SMR.data $Vgll3 <- dplyr::recode(All.SMR.data $Call_vgll3, E = "EE", L = "LL")
All.SMR.data $Six6 <- dplyr::recode(All.SMR.data $Call_six6, E = "EE", L = "LL")

## Edit genotypes E and L to -0.5 or 0.5
All.SMR.data $Call_vgll3 <- dplyr::recode(All.SMR.data $Call_vgll3, E = -0.5, L = 0.5)
All.SMR.data $Call_six6 <- dplyr::recode(All.SMR.data $Call_six6, E = -0.5, L = 0.5)

## Edit sex F and M to -0.5 or 0.5
All.SMR.data $Sex <- dplyr::recode(All.SMR.data $Sex, "F" = -0.5, M = 0.5)

## Edit treatment High food and Low food to -0.5 or 0.5
All.SMR.data $Treatment_centred <- dplyr::recode(All.SMR.data $Treatment, "High food" = -0.5, "Low food" = 0.5)

## Remove unused columns
All.SMR.data <- All.SMR.data %>%
  select(-c(q0.1.SMR.resid, q0.2.SMR.resid, Temp, tank))

#write.table(All.SMR.data, file= "Data/Analysed.SMR.data.txt", row.names = F,
#            quote = F, dec = ".", sep = "\t")

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
### Full model with interactions
#########################################
# Chamber removed because of singularity (cannot be estimated)
SMR.full.mod <-lmer(log10(SMR.abs.mlnd) ~  
                  Treatment_centred +
                  Sex +
                  Call_vgll3 + 
                  Call_six6 + 
                  scale(log10(Mass))+
                  scale(log10(Mass)):Treatment_centred +
                  Treatment_centred:Call_vgll3 +
                  Treatment_centred:Call_six6 +
                  Sex:Call_vgll3 +
                  Sex:Call_six6 +
                  Call_vgll3:Call_six6+
                  (1|Family) +
                  (1|Batch) ,
                REML = F, 
                data= All.SMR.data,
                na.action = na.fail
)

outlierTest(SMR.full.mod)
# Strong outlier:
All.SMR.data[67,]
# Filter this fish from the data
All.SMR.data<- filter(All.SMR.data,
                      !(PIT == "A00000D0900001897006087"))

# rerun model
#save(All.SMR.data, file= "Data/All.SMR.data")

summary(SMR.full.mod)
anova(SMR.full.mod)

## Check assumptions
SMR.mod.res <- cbind(All.SMR.data, resid(SMR.full.mod))

# Residuals and fitted values with lines to help visualization:
plot(fitted(SMR.full.mod),resid(SMR.full.mod, type="pearson"))
mywhiskers(fitted(SMR.full.mod),resid(SMR.full.mod, type="pearson"), add=T, se=F)

# heteroscedasticity tests
leveneTest(resid(SMR.full.mod)~ Treatment, data =  SMR.mod.res)
boxplot(SMR.mod.res$`resid(SMR.full.mod)`~SMR.mod.res$Treatment) # a little heteroscedastic

leveneTest(resid(SMR.full.mod)~ as.factor(Sex), data =  SMR.mod.res)
leveneTest(resid(SMR.full.mod)~ Vgll3, data =  SMR.mod.res)
leveneTest(resid(SMR.full.mod)~ Six6, data =  SMR.mod.res)

## plot diagnostics
plot_model(SMR.full.mod, type = "diag")

### Export estimates and Type III results
sum.SMR.full.mod <- summary(SMR.full.mod)
aov.SMR.full.mod <- anova(SMR.full.mod)
confs.full <- confint(SMR.full.mod, oldNames = FALSE)

## Make a summary table with estimates and type 3 effects
table.SMR.full.mod <- data.frame(Coefficient = row.names(sum.SMR.full.mod$coefficients),
                                Estimate = sum.SMR.full.mod$coefficients[,1],
                                SE = sum.SMR.full.mod$coefficients[,2],
                                SSq =  c(NA, aov.SMR.full.mod$`Sum Sq`),
                                `Den DF` = c(sum.SMR.full.mod $coefficients[,3][1], aov.SMR.full.mod$DenDF),
                                F = c(sum.SMR.full.mod $coefficients[,4][1], aov.SMR.full.mod$`F value`),
                                p = c(sum.SMR.full.mod $coefficients[,5][1], aov.SMR.full.mod$`Pr(>F)`))

table.SMR.full.mod <- mutate(table.SMR.full.mod, across(-1, round, 5))
table.SMR.full.mod

# Table of random effects:
re.table.SMR.full <- as.data.frame(sum.SMR.full.mod$varcor) %>%
  dplyr::select(-(c(var1, var2, sdcor))) %>% #drop unnecessary cols
  rename("Random effect" = "grp", "Var" = "vcov") %>%
  mutate("C.I.low" = c(confs.full[1,1]^2, confs.full[2,1]^2, NA),
         "C.I.high" = c(confs.full[1,2]^2, confs.full[2,2]^2, NA),) %>%
  mutate(across(2:4, round, 4))

re.table.SMR.full

# export both into xlsx
write.xlsx(table.SMR.full.mod, file = "table.SMR.full.mod.xlsx")
write.xlsx(re.table.SMR.full, file = "re.table.SMR.full.xlsx")

# Pairwise differences between six6 at food levels
emm.smr <- emmeans(SMR.full.mod, specs = c("Treatment_centred", "Call_six6"))
pairs(emm.smr)

##################################################
## Get a more parsimonious model with MuMIn
##################################################
all.mumin.SMR <-  dredge(SMR.full.mod,rank = "AICc")

subset.mumin.SMR <- subset(all.mumin.SMR, delta <2)

# save as table
subset.mumin.SMR.table <- as.data.frame(subset.mumin.SMR)
subset.mumin.SMR.table <- subset.mumin.SMR.table %>%
  mutate(across(-"df", round, 3))

subset.mumin.SMR.table

write.xlsx(subset.mumin.SMR.table, file = "subset.mumin.SMR.table.xlsx")

## Calculate average model
ave.model <- model.avg(subset.mumin.SMR, fit =T)
names(ave.model)
summary(ave.model)
#save(ave.model, file = "ave_model_AS")
#load("ave_model_AS")

# Take averages using the full model method
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
write.xlsx(ave.model.table, file = "ave.model.table.SMR.xlsx")

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
plot_scaling.smr <- interactions::interact_plot(SMR.full.mod, pred = Mass, modx = Treatment_centred, modx.values = c(-0.5, 0.5),
                            modx.labels = c("High food", "Low food"), interval = T, 
                            int.type = "confidence", plot.points = T, partial.residuals = T, 
                            point.alpha = 0.4, colors = c("#2A1A3C", "#F3771A")) +
  coord_trans(x = "log") +
  ylab(expression(paste("Log"[10], " SMR (mg ", O[2], "/h)"))) +
  xlab(expression(paste("Body mass (g)"))) +
  theme(legend.position = c(0.78,0.2),
        legend.title = element_blank())

plot_scaling.smr

ggsave(plot_scaling.smr, file= "plot_scaling_smr.tiff",  width = 8.5, height = 7, units = "cm")

##################################################
## Predicted means
##################################################
mean(All.SMR.data$Mass)

SMR.predict <- ggpredict(SMR.full.mod, terms =  c("Call_vgll3", "Call_six6" ),
                          ci.lvl = 0.9, condition =c("Mass" = 3.68),
                          back.transform = FALSE)
SMR.predict
plot(SMR.predict)

# change predicted SMR to linear scale (mg/O2/h)
SMR.predict[,2:5]<- 10^SMR.predict[,2:5]

# divide predictions by mass in kg
SMR.predict[,2:5]<- SMR.predict[,2:5]/0.00368

# Export predictions, plot together with MMR and AS
save(SMR.predict, file = "SMR.predict") 


