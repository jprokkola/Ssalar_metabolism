######################################################################
## Linear mixed models for absolute aerobic scope in salmon
## Jenni Prokkola 2021
######################################################################

library(tidyverse)
library(data.table)
library(ggplot2)
library(lme4)
library(lmerTest)
library(car)
library(ggeffects)
library(merTools)
library(sjPlot)
library(ggcorrplot)
library(ggpubr)
library(viridis)
library(openxlsx)

# Data from AS_dataprep: all.AS.data
load("Data/all.AS.data")
str(all.AS.data)

## Plot AS vs mass
ggplot(all.AS.data, aes(y= log10(absAS), x=log10(Mass)))+
  facet_wrap(~Treatment) +
  geom_point(alpha=0.5) +
  stat_smooth(method = "lm", level = 0.95)+
  ylab(expression(paste("Log"[10], " abs AS (mg ", O[2], "/h)"))) +
  xlab(expression(paste("Log"[10], " body mass (g)"))) +
  theme_bw()

# Scaling exponents
all.AS.data %>%
  group_by(Treatment) %>%
  summarise(b = coef(lm(log10(absAS)~log10(Mass)))[2],
            R2 = summary(lm(log10(absAS)~log10(Mass)))$r.squared)
  
## Main effects of genotypes, sex, treatment and mass on AS.
absAS.main.mod <- lmer(na.action = "na.fail" ,
                  log10(absAS) ~
                    Treatment_centred +
                    Sex +
                    Call_vgll3 +
                    Call_six6 +
                    scale(Order)+
                    scale(log10(Mass))+
                    (1|Family) +
                    (1|MMR.Chamber) +
                    (1|initial),
                  REML = F,
                  data = all.AS.data)

plot_model(absAS.main.mod, type = "diag")
outlierTest(absAS.main.mod)

summary(absAS.main.mod)
anova(absAS.main.mod)

### test a full model including interactions, and exclude interaction terms based on p-values

absAS.mod.full <- lmer(na.action = "na.fail" ,
                  log10(absAS) ~
                    Treatment_centred +
                    Sex +
                    Call_vgll3 +
                    Call_six6 +
                    Call_vgll3:Call_six6+
                    scale(log10(Mass))+
                    scale(Order)+
                    Treatment_centred:scale(log10(Mass))+
                    Treatment_centred:Call_vgll3 +
                    Treatment_centred:Call_six6 +
                    Sex : Call_vgll3+
                    Sex: Call_six6 +
                    (1|Family) +
                    (1|MMR.Chamber) +
                    (1|initial),
                  REML = F,
                  data = all.AS.data)

anova(absAS.mod.full)

# Save table
sum.absAS.full <- summary(absAS.mod.full)
aov.absAS.full <- anova(absAS.mod.full)
confs.full <- confint(absAS.mod.full, oldNames=FALSE)

## Make a summary table with estimates and type 3 effects
table.absAS.full <- data.frame(Coefficient = row.names(sum.absAS.full$coefficients),
                           Estimate = sum.absAS.full$coefficients[,1],
                           SE = sum.absAS.full$coefficients[,2],
                           SSq =  c(NA, aov.absAS.full$`Sum Sq`),
                           `Den DF` = c(sum.absAS.full$coefficients[,3][1], aov.absAS.full$DenDF),
                           F = c(sum.absAS.full$coefficients[,4][1], aov.absAS.full$`F value`),
                           p = c(sum.absAS.full$coefficients[,5][1], aov.absAS.full$`Pr(>F)`))

table.absAS.full <- mutate(table.absAS.full, across(-1, round, 5))

table.absAS.full

# Table of random effects 
re.table.absAS.full <- as.data.frame(sum.absAS.full$varcor) %>%
  dplyr::select(-(c(var1, var2, sdcor))) %>% #drop unnecessary cols
  rename("Random effect" = "grp", "Var" = "vcov") %>%
  mutate("C.I.low" = c(confs.full[1,1]^2, confs.full[2,1]^2, confs.full[3,1]^2, NA),
       "C.I.high" = c(confs.full[1,2]^2, confs.full[2,2]^2,confs.full[3,2]^2, NA),) %>%
  mutate(across(2:4, round, 4))

re.table.absAS.full

# export both into xlsx

write.xlsx(table.absAS.full, file = "table.absAS.full.xlsx")
write.xlsx(re.table.absAS.full, file = "re.table.absAS.full.xlsx")

## Simplify model. Remove Treatment:Call_six6
absAS.mod <- lmer(na.action = "na.fail" ,
                  log10(absAS) ~
                    Treatment_centred +
                    Sex +
                    Call_vgll3 +
                    Call_six6 +
                    Call_vgll3:Call_six6+
                    scale(log10(Mass))+
                    scale(Order)+
                    Treatment_centred:scale(log10(Mass))+
                    Treatment:Call_vgll3 +
                    Sex : Call_vgll3+
                    Sex: Call_six6 +
                    (1|Family) +
                    (1|MMR.Chamber)+
                    (1|initial),
                  REML = F,
                  data = all.AS.data)

anova(absAS.mod)
confint(absAS.mod)

# Remove Sex:Call_six6

absAS.mod <- lmer(na.action = "na.fail" ,
                  log10(absAS) ~
                    Treatment_centred +
                    Sex +
                    Call_vgll3 +
                    Call_six6 +
                    Call_vgll3:Call_six6+
                    scale(log10(Mass))+
                    scale(Order)+
                    Treatment_centred:scale(log10(Mass))+
                    Treatment:Call_vgll3 +
                    Sex : Call_vgll3+
                    (1|Family) +
                    (1|MMR.Chamber)+
                    (1|initial),
                  REML = F,
                  data = all.AS.data)

anova(absAS.mod)
confint(absAS.mod)

# Remove Treatment:Call_vgll3
absAS.mod <- lmer(na.action = "na.fail" ,
                  log10(absAS) ~
                    Treatment_centred +
                    Sex +
                    Call_vgll3 +
                    Call_six6 +
                    Call_vgll3:Call_six6+
                    scale(log10(Mass))+
                    scale(Order)+
                    Treatment_centred:scale(log10(Mass))+
                    Sex : Call_vgll3+
                    (1|Family) +
                    (1|MMR.Chamber)+
                    (1|initial),
                  REML = F,
                  data = all.AS.data)

anova(absAS.mod)
summary(absAS.mod)

# Remove genotype interaction
absAS.mod <- lmer(na.action = "na.fail" ,
                  log10(absAS) ~
                    Treatment_centred +
                    Sex +
                    Call_vgll3 +
                    Call_six6 +
                    scale(log10(Mass))+
                    scale(Order)+
                    Treatment_centred:log10(Mass)+
                    Sex : Call_vgll3+
                    (1|Family) +
                    (1|MMR.Chamber)+
                    (1|initial),
                  REML = F,
                  data = all.AS.data)

anova(absAS.mod)
confint(absAS.mod)

# remove Sex:Call_vgll3

absAS.mod <- lmer(na.action = "na.fail" ,
                  log10(absAS) ~
                    Treatment_centred +
                    Sex +
                    Call_vgll3 +
                    Call_six6 +
                    scale(Order)+
                    scale(log10(Mass))+
                    Treatment_centred:scale(log10(Mass))+
                    (1|Family) +
                    (1|MMR.Chamber) +
                    (1|initial),
                  REML = F,
                  data = all.AS.data)

anova(absAS.mod)

# remove order
absAS.mod.final <- lmer(na.action = "na.fail" ,
                  log10(absAS) ~
                    Treatment_centred +
                    Sex +
                    Call_vgll3 +
                    Call_six6 +
                    scale(log10(Mass))+
                    Treatment_centred:scale(log10(Mass))+
                    (1|Family) +
                    (1|MMR.Chamber) +
                    (1|initial),
                  REML = F,
                  data = all.AS.data)

sum.absAS.final <- summary(absAS.mod.final)
aov.absAS.final <- anova(absAS.mod.final)
confs.final <- confint(absAS.mod.final, oldNames=FALSE)

## Check variance partitioning, need to run again with mass-corrected data (which is the relevant response)
all.AS.data$resid_AS <- resid(lm(log10(absAS) ~ log10(Mass), data = all.AS.data))

AS.mod.part <- lmer(na.action = "na.fail" ,
                     resid_AS ~
                      Treatment_centred +
                      Sex +
                      Call_vgll3 +
                      Call_six6 +
                      (1|Family) +
                      (1|MMR.Chamber) +
                      (1|initial),
                     REML = F,
                     data = all.AS.data)

partR2_AS <- partR2(AS.mod.part, partvars = c("Call_vgll3"), nboot = 1000)
R2s_AS <- partR2_AS$R2
write.xlsx(R2s_AS, file ="AS_R2s.xlsx")

## Make a summary table with estimates and type 3 effects
table.absAS.final <- data.frame(Coefficient = row.names(sum.absAS.final$coefficients),
                               Estimate = sum.absAS.final$coefficients[,1],
                               SE = sum.absAS.final$coefficients[,2],
                               SSq =  c(NA, aov.absAS.final$`Sum Sq`),
                               `Den DF` = c(sum.absAS.final$coefficients[,3][1], aov.absAS.final$DenDF),
                               F = c(sum.absAS.final$coefficients[,4][1], aov.absAS.final$`F value`), #First value is for intercept
                               p = c(sum.absAS.final$coefficients[,5][1], aov.absAS.final$`Pr(>F)`))

table.absAS.final <- mutate(table.absAS.final, across(-1, round, 5))

table.absAS.final

# Table of random effects:
re.table.absAS.final <- as.data.frame(sum.absAS.final$varcor) %>%
  dplyr::select(-(c(var1, var2, sdcor))) %>% #drop unnecessary cols
  rename("Random effect" = "grp", "Var" = "vcov") %>%
  mutate("C.I.low" = c(confs.final[1,1]^2, confs.final[2,1]^2, confs.final[3,1]^2, NA),
         "C.I.high" = c(confs.final[1,2]^2, confs.final[2,2]^2,confs.final[3,2]^2, NA)) 

re.table.absAS.final

# export both into xlsx

write.xlsx(table.absAS.final, file = "table.absAS.final.xlsx")
write.xlsx(re.table.absAS.final, file = "re.table.absAS.final.xlsx")


### Does adding rSMR to the model affect genotype result?
all.AS.data$rSMR.mass <- resid(lm(log10(SMR.abs.mlnd) ~ log10(Mass), data = all.AS.data))

absAS.mod.smr <- lmer(na.action = "na.fail" ,
                        log10(absAS) ~
                        scale(all.SMR.resid) + 
                        Treatment_centred +
                        Sex +
                        Call_vgll3 +
                        Call_six6 +
                        scale(log10(Mass))+
                        Treatment_centred:scale(log10(Mass))+
                        (1|Family) +
                        (1|MMR.Chamber) +
                        (1|initial),
                        REML = F,
                        data = all.AS.data)

anova(absAS.mod.smr)
summary(absAS.mod.smr)

sum.absAS.smr <- summary(absAS.mod.smr)
aov.absAS.smr <- anova(absAS.mod.smr)
confs.as.smr <- confint(absAS.mod.smr, oldNames=FALSE)


## Make a summary table with estimates and type 3 effects
table.absAS.smr <- data.frame(Coefficient = row.names(sum.absAS.smr$coefficients),
                                Estimate = sum.absAS.smr$coefficients[,1],
                                SE = sum.absAS.smr$coefficients[,2],
                                SSq =  c(NA, aov.absAS.smr$`Sum Sq`),
                                `Den DF` = c(sum.absAS.smr$coefficients[,3][1], aov.absAS.smr$DenDF),
                                F = c(sum.absAS.smr$coefficients[,4][1], aov.absAS.smr$`F value`), #First value is for intercept
                                p = c(sum.absAS.smr$coefficients[,5][1], aov.absAS.smr$`Pr(>F)`))

table.absAS.smr <- mutate(table.absAS.smr, across(-1, round, 5))

table.absAS.smr

# Table of random effects:
re.table.absAS.smr <- as.data.frame(sum.absAS.smr$varcor) %>%
  dplyr::select(-(c(var1, var2, sdcor))) %>% #drop unnecessary cols
  rename("Random effect" = "grp", "Var" = "vcov") %>%
  mutate("C.I.low" = c(confs.as.smr[1,1]^2, confs.as.smr[2,1]^2,confs.as.smr[3,1]^2, NA),
         "C.I.high" = c(confs.as.smr[1,2]^2, confs.as.smr[2,2]^2,confs.as.smr[3,2]^2, NA))

re.table.absAS.smr

# export both into xlsx

write.xlsx(table.absAS.smr, file = "table.absAS.smr.xlsx")
write.xlsx(re.table.absAS.smr, file = "re.table.absAS.smr.xlsx")

### Levene's tests
absAS.mod.res <- cbind(all.AS.data, resid(absAS.mod.final))

leveneTest(resid(absAS.mod.final)~ Treatment, data =  absAS.mod.res)
leveneTest(resid(absAS.mod.final)~ as.factor(Sex), data =  absAS.mod.res)
leveneTest(resid(absAS.mod.final)~ Vgll3, data =  absAS.mod.res)
leveneTest(resid(absAS.mod.final)~ Six6, data =  absAS.mod.res)

##################################################
##SMR-mass scaling plot using partial residuals
################################################## 

plot_scaling <- interactions::interact_plot(absAS.mod.final, pred = Mass, modx = Treatment_centred, 
                                            modx.values = c(-0.5, 0.5),
                                            modx.labels = c("High food", "Low food"), interval = T, 
                                            int.type = "confidence", plot.points = T, partial.residuals = T, 
                                            point.alpha = 0.4, colors = c("#2A1A3C", "#F3771A"))  +
  coord_trans(x = "log") +
  ylab(expression(paste("Log"[10], " AS (mg ", O[2], "/h)"))) +
  xlab(expression(paste("Body mass (g)"))) +
  theme(legend.position = c(0.78,0.2),
        legend.title = element_blank())

ggsave(plot_scaling, file= "plot_scaling_AS.tiff",  width = 8.5, height = 7, units = "cm")

##################################################
## Predicted means
##################################################
### The means of vgll3 EE and LL
mean(all.AS.data$Mass)
predict.vgll3 <- ggpredict(absAS.mod.final, terms = c("Call_vgll3"), type = "fe", ci.lvl = 0.9,
                           condition =c("Mass" = 3.7), back.transform = FALSE)

predict.vgll3

predict.vgll3[,2:5]<- 10^predict.vgll3[,2:5]
# Calculate in mg/O2/kg
predict.vgll3[,2:5]<- predict.vgll3[,2:5]/0.0037
predict.vgll3

# Predictions from final model for six6 and vgll3
absAS.predict <- ggpredict(absAS.mod.final, terms = c("Call_vgll3", "Call_six6"), type = "fe", 
                           ci.lvl = 0.9,condition =c("Mass" = 3.7),
                           back.transform = FALSE)

plot(absAS.predict)

# Manually back-transform the values to mg/O2 (per 3.7g)
absAS.predict[,2:5]<- 10^absAS.predict[,2:5]
# Calculate in mg/O2/kg
absAS.predict[,2:5]<- absAS.predict[,2:5]/0.0037

# save predictions for plotting separately
save(absAS.predict, file = "absAS.predict") 

##################################################
## Correlation matrix of metabolic traits
##################################################
# Separated by treatment
# metabolic trait correlations (from family- and mass-adjusted residuals)
corrdata1 <- all.AS.data %>%
  filter(Treatment == "High food") %>%
  dplyr::select(c("all.SMR.resid", "all.MMR.resid", "absAS.resid"))


names(corrdata1) <- c("rSMR", "rMMR", "rAbsAS")
corr1 <- round(cor(corrdata1), 2)
# p-values for correlations
p.mat <- cor_pmat(corrdata1)

corrdata2 <- all.AS.data %>%
  filter(Treatment == "Low food") %>%
  dplyr::select(c("all.SMR.resid", "all.MMR.resid", "absAS.resid"))

names(corrdata2) <- c("rSMR", "rMMR", "rAbsAS")
corr2 <- round(cor(corrdata2), 2)
# p-values for correlations
p.mat2 <- cor_pmat(corrdata2)

# Scatter plot of SMR and MMR
smr_mmr <- ggplot(all.AS.data, aes(y= all.SMR.resid, x=all.MMR.resid))+
  facet_wrap(~Treatment)+
  geom_point(alpha=0.4) +
  ylab("rSMR") +
  xlab("rMMR") +
  scale_x_continuous(breaks=seq(-0.10,0.10,0.1))+
  theme_bw()+
  theme(legend.position = "none",
        axis.title = element_text(size=9))
ggsave(smr_mmr, file= "Fig3.png", width = 9, height = 7, units = "cm")

##################################################
# Get N in each group
##################################################


AAS_N <- all.AS.data %>%
  group_by(Call_vgll3, Call_six6) %>%
  summarise(N = length(absAS[is.na(absAS)=="FALSE"])) %>%
  rename(x = Call_vgll3, group = Call_six6)

AAS_allfamN<-all.AS.data %>%
  group_by(Family,Treatment,Call_vgll3, Call_six6) %>%
  summarise(N = length(absAS[is.na(absAS)=="FALSE"])) %>%
  rename(x = Call_vgll3, group = Call_six6)

all.AS.data %>%
  group_by(Treatment) %>%
  summarise(N = length(absAS[is.na(absAS)=="FALSE"]))

## write table
AS.data.analysed <- all.AS.data %>%
  select(-c(all.SMR.resid, all.MMR.resid, absAS.resid, q0.1.SMR.resid,
            q0.2.SMR.resid, mlndSMR.resid , SMR.mass.mlnd, SMR.abs.mlnd,
            MR.abs, MR.mass )) %>%
  rename("absAS_mgO2_h" = "absAS")

write.table(AS.data.analysed, file= "Data/Analysed.AbsAS.data.txt", row.names = F,quote = F, dec = ".", sep = "\t")

##################################################
## Plot mass-adjusted data
##################################################

all.AS.data$Vgll3 <- dplyr::recode(all.AS.data$Vgll3, "EE" = "Early",  "LL" = "Late")
all.AS.data$Six6 <- dplyr::recode(all.AS.data$Six6, "EE" = "Early",  "LL" = "Late")

AS_plot <- ggplot(all.AS.data, aes(y= absAS.resid, x= Vgll3, shape = Six6))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(.3), 
             size =2, alpha = 0.5, color = "#420A68FF") +
  # y-axis title is the plot title
  ggtitle(label = "Residual aerobic scope")+
  xlab("vgll3")+
 # labs(colour = "six6")+
 # scale_color_viridis(begin = 0.19, end = 0.7, discrete = T, option = "inferno") +
  theme_bw()+
  theme(panel.spacing.x=unit(0, "lines"),
        plot.title = element_text(size =10),
        plot.subtitle =  element_text(size =9),
        axis.text.y = element_text(size=10),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=10, vjust=0.6),
        panel.grid.minor.x =element_blank(),
        panel.grid.major.x=element_blank(),
        legend.position = "none")

ggsave(AS_plot, filename = "AS_plot.png", width = 7, height = 12, units = "cm", dpi = 300)

