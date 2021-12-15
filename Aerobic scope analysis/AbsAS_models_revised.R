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
library(openxlsx)
library(MuMIn)
library(emmeans)

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
  
#########################################
### Full model with interactions
#########################################

absAS.mod.full <- lmer(na.action = "na.fail" ,
                  log10(absAS) ~
                    Treatment_centred +
                    Sex +
                    Call_vgll3 +
                    Call_six6 +
                    scale(log10(Mass))+
                    scale(Order)+
                    Treatment_centred:scale(log10(Mass))+
                    Treatment_centred:Call_vgll3 +
                    Treatment_centred:Call_six6 +
                    Sex:Call_vgll3+
                    Sex:Call_six6 +
                    Call_vgll3:Call_six6 +
                    (1|Family) +
                    (1|MMR.Chamber) +
                    (1|initial),
                  REML = F,
                  data = all.AS.data)

## plot diagnostics
plot_model(absAS.mod.full, type = "diag")

## outlier test
outlierTest(absAS.mod.full)

# Residuals and fitted values with lines to help visualization:
absAS.mod.res <- cbind(all.AS.data, resid(absAS.mod.full))
plot(fitted(absAS.mod.full),resid(absAS.mod.full, type="pearson"))
mywhiskers(fitted(absAS.mod.full),resid(absAS.mod.full, type="pearson"), add=T, se=F)

# heteroscedasticity tests
leveneTest(resid(absAS.mod.full)~ Treatment, data =  absAS.mod.res)
leveneTest(resid(absAS.mod.full)~ as.factor(Sex), data =  absAS.mod.res)
leveneTest(resid(absAS.mod.full)~ Vgll3, data =  absAS.mod.res)
leveneTest(resid(absAS.mod.full)~ Six6, data =  absAS.mod.res)

### Export estimates and Type III results
sum.AS.full.mod<- summary(absAS.mod.full)
aov.AS.full.mod<-anova(absAS.mod.full)
confs.full <- confint(absAS.mod.full, oldNames = FALSE)

## Make a summary table with estimates and type 3 effects
table.AS.full <- data.frame( Coefficient = row.names(sum.AS.full.mod$coefficients),
                              Estimate = sum.AS.full.mod$coefficients[,1],
                              SE = sum.AS.full.mod$coefficients[,2],
                              SSq =  c(NA, aov.AS.full.mod$`Sum Sq`),
                              `Den DF` = c(sum.AS.full.mod$coefficients[,3][1], aov.AS.full.mod$DenDF),
                              F = c(sum.AS.full.mod$coefficients[,4][1], aov.AS.full.mod$`F value`),
                              p = c(sum.AS.full.mod$coefficients[,5][1], aov.AS.full.mod$`Pr(>F)`))

# too many digits, round the numbers
table.AS.full <- mutate(table.AS.full, across(-1, round, 5))

table.AS.full

# Table of random effects:
re.table.AS.full <- as.data.frame(sum.AS.full.mod$varcor) %>%
  dplyr::select(-(c(var1, var2, sdcor))) %>% #drop unnecessary cols
  rename("Random effect" = "grp", "Var" = "vcov") %>%
  mutate("C.I.low" = c(confs.full[1,1]^2, confs.full[2,1]^2,confs.full[3,1]^2, NA),
         "C.I.high" = c(confs.full[1,2]^2, confs.full[2,2]^2, confs.full[3,2]^2,NA),) %>%
  mutate(across(2:4, round, 4))

re.table.AS.full

# export both to xlsx

write.xlsx(table.AS.full, file = "table.AS.full.xlsx")
write.xlsx(re.table.AS.full, file = "re.table.AS.full.xlsx")

##################################################
## Get a more parsimonious model with MuMIn
##################################################
all.mumin.AS <-  dredge(absAS.mod.full,rank = "AICc")
  
subset.mumin.AS <- subset(all.mumin.AS, delta <2)

# save as table
sub.mumin.AS.table <- as.data.frame(subset.mumin.AS)
sub.mumin.AS.table <- sub.mumin.AS.table %>%
  mutate(across(-"df", round, 3))
        
sub.mumin.AS.table

write.xlsx(sub.mumin.AS.table, file = "sub.mumin.AS.table.xlsx")

## Calculate average model
ave.model <- model.avg(subset.mumin.AS, fit =T)
names(ave.model)
summary(ave.model)
#save(ave.model, file = "ave_model_AS")
#load("ave_model_AS")

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
write.xlsx(ave.model.table, file = "ave.model.table.AS.xlsx")

##################################################
## Does adding rSMR to the model affect genotype result?
##################################################
all.AS.data$rSMR.mass <- resid(lm(log10(SMR.abs.mlnd) ~ log10(Mass), data = all.AS.data))

absAS.mod.smr <- lmer(na.action = "na.fail" ,
                                        log10(absAS) ~
                                          Treatment_centred +
                                          Sex +
                                          Call_vgll3 +
                                          Call_six6 +
                                          scale(log10(Mass))+
                                          scale(Order)+
                                          scale(rSMR.mass) +
                                          Treatment_centred:scale(log10(Mass))+
                                          Treatment_centred:Call_vgll3 +
                                          Treatment_centred:Call_six6 +
                                          Sex:Call_vgll3+
                                          Sex:Call_six6 +
                                          Call_vgll3:Call_six6 +
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

# export to xlsx

write.xlsx(table.absAS.smr, file = "table.absAS.smr.xlsx")
write.xlsx(re.table.absAS.smr, file = "re.table.absAS.smr.xlsx")

##################################################
## Check variance partitioning
##################################################
## need to run model with mass-corrected data (the relevant response)

all.AS.data$resid_AS <- resid(lm(log10(absAS) ~ log10(Mass), data = all.AS.data))

AS.mod.part <- lmer(na.action = "na.fail" ,
                    resid_AS ~
                      Call_vgll3 +
                      (1|Family) +
                      (1|MMR.Chamber) +
                      (1|initial),
                    REML = F,
                    data = all.AS.data)

partR2_AS <- partR2(AS.mod.part, partvars = c("Call_vgll3"), nboot = 1000)
R2s_AS <- partR2_AS$R2
write.xlsx(R2s_AS, file ="AS_R2s.xlsx")

##################################################
##SMR-mass scaling plot using partial residuals
################################################## 

plot_scaling <- interactions::interact_plot(absAS.mod.full, pred = Mass, modx = Treatment_centred, 
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
## Predicted means for genotypes
##################################################
mean(all.AS.data$Mass) # Use fixed mass
vgll3.predict <- ggpredict(absAS.mod.full, terms = c("Call_vgll3"), type = "fe", 
                           ci.lvl = 0.9,condition =c("Mass" = 3.7),
                           back.transform = FALSE)

# Manually back-transform the values to mg/O2 (per 3.7g)
vgll3.predict[,2:5]<- 10^vgll3.predict[,2:5]
# Calculate in mg/O2/kg
vgll3.predict[,2:5]<- vgll3.predict[,2:5]/0.0037

# Predictions for six6 and vgll3
absAS.predict <- ggpredict(absAS.mod.full, terms = c("Call_vgll3", "Call_six6"), type = "fe", 
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
ggsave(smr_mmr, file= "rSMR_rMMR.png", width = 9, height = 7, units = "cm")

##################################################
# Get N in each group
##################################################

AS_N <- all.AS.data %>%
  group_by(Call_vgll3, Call_six6) %>%
  summarise(N = length(absAS[is.na(absAS)=="FALSE"])) %>%
  rename(x = Call_vgll3, group = Call_six6)

AS_allfamN<-all.AS.data %>%
  group_by(Family,Treatment,Call_vgll3, Call_six6) %>%
  summarise(N = length(absAS[is.na(absAS)=="FALSE"])) %>%
  rename(x = Call_vgll3, group = Call_six6)

all.AS.data %>%
  group_by(Treatment) %>%
  summarise(N = length(absAS[is.na(absAS)=="FALSE"]))

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




