#######################################################################################
## Compilation and annotation of Figure with SMR, MMR and AS predicted means and CI's.
#######################################################################################

library(ggplot2)
library(tidyverse)
library(viridis)
library(ggpubr)

##### Make one figure with all 3 together using prediction tables
load("SMR.predict")
SMR.predict
load("MMR.predict")
MMR.predict
load("absAS.predict")
absAS.predict

## Combine and recode genotypes as Early and Late
all.pred <- bind_rows(SMR.predict, MMR.predict, absAS.predict, .id = "Var") 
all.pred$Var <- dplyr::recode(all.pred$Var, "1" = "SMR",  "2" = "MMR", "3" =  "AS")
all.pred$x <- dplyr::recode(all.pred$x, "-0.5" = "Early",  "0.5" = "Late")
all.pred$group <- dplyr::recode(all.pred$group, "-0.5" = "Early",  "0.5" = "Late")

## Make data frame for annotating plot with significance values
# Vgll3 effect in six6 LL
annot_text_vg_MMR <- data.frame(
  group = "Late",
  Var = "MMR",
  group1= "Early",
  group2 ="Late",
  p = c("p = 0.0001") 
)

# six6 effect in vgll3 LL
annot_text_s6_MMR <- data.frame(
  x = "Late",
  group = "Late",
  Var = "MMR",
  group1= "Early",
  group2 ="Late",
  p = c("p = 0.020")
)

#Other p-values added in powerpoint because adding them to vgll3 main effect not feasible :/

All.pred.plot<-
  ggplot(all.pred, aes(y= predicted, x=x, shape = group, colour = Var))+
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(.3),
                width=0.2) +
  geom_point(position = position_dodge(.3), 
             stat = "identity",size =2) +
  # y-axis title is the plot title
  ggtitle(label = expression(paste("Predicted mean, mg ", O[2], "/kg/h (90% CI)"))) +
  xlab("vgll3")+
  labs(shape = "six6")+
  scale_y_continuous(limits = c(100,700), breaks = scales::extended_breaks(10))+
  scale_color_viridis(begin = 0.2, end = 0.7, discrete = TRUE, option = "inferno") +
  
  stat_pvalue_manual(annot_text_vg_MMR, 
    y.position = 690, step.increase = 0.1,
    label = "p", tip.length = 0.02, size = 3,
    position = position_nudge(0.09), xmin = "group1") +
  
  stat_pvalue_manual(annot_text_s6_MMR, 
    y.position = 656, x = "group2",
    label = "p", size = 3, tip.length = 0.5,
    position = position_dodge(0.3)) +
  
  annotate("text", label = "SMR", x = 2.4, y = 125, colour = "#F3771AFF")+
  annotate("text", label = "MMR", x = 2.4, y = 590, colour = "#A82E5FFF")+
  annotate("text", label = "AS", x = 2.4, y = 430, colour = "#420A68FF")+
  guides(color = FALSE, shape = guide_legend(override.aes = list(linetype = "blank")))+
  theme_bw()+
  theme(panel.spacing.x=unit(0, "lines"),
        plot.title = element_text(size =10),
        plot.subtitle =  element_text(size =9),
        axis.text.y = element_text(size=10),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=10, vjust=0.6),
        panel.grid.minor.x =element_blank(),
        panel.grid.major.x=element_blank(),
        legend.position = "right")

ggsave(All.pred.plot, filename = "All.pred.plot.pdf", width = 9, height = 12, units = "cm")
ggsave(All.pred.plot, filename = "All.pred.plot.png", width = 9, height = 12, units = "cm", dpi = 300)




