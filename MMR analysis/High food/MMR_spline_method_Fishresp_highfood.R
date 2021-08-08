######################################################################
## MMR calculation from slopes using FishResp
## Slopes obtained using tangent from the beginning of 
## a spline slope (spline method)
######################################################################
library(FishResp)
library(tidyverse)
library(ggplot2)
library(lubridate)

load("Data/MMR_HF_data5") # A list of all O2 data for the MMR slope for each fish
head(HF_data5[[1]])

## Plot all slopes with a spline slope. We will use the tangent from the beginning of the slope as MMR. 
## This is better than slopes based on e.g. 10 points as less sensitive to random variation around the slope.

# Function for plot:
rawdat.multiplot <- function(data) {
  lapply(data, function(x){
    plotdat = x
    DD = 1:length ( plotdat$O2.correct )
    MM = plotdat$O2.correct 
    spline <- smooth.spline(plotdat$O2.correct  ~ 1:length ( plotdat$O2.correct)  , df = 10  )
    plotdat$spl <- spline$y
    p <- ggplot(data = plotdat, 
                aes(x = Time, y = O2.correct ))+
      geom_point(cex = 0.4)+
      geom_line(aes(x = Time, y = spl ), color = "green")+ 
      ylab("O2 mg")+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 90, size = 8),
            axis.title.x = element_blank())
    
    print(p)
    rm(plotdat)
    rm(slope)
    rm(i)
  })
}

# Plot and check plots for slope quality (if too much noise, exclude. Some can have pump malfunctions or bubbles)
pdf(file = "MMR_HF_slopes_splines.pdf", width = 5, height = 5)
rawdat.multiplot(MMR_HF_data5)
dev.off()

## Bad ones: 23, 126
head(HF_data5[[23]]) #known low pump speed, not comparable
head(HF_data5[[126]]) #nothing known, but not good

### Spline method to get tangents at 1 and 20 seconds (T1 and T2). This also plots all slopes.
MMR_HF_data6 = t(sapply ( 1: length(HF_data5) , function(i) {
  DD = 1:length ( HF_data5[[i]]$O2.correct )
  MM = ( HF_data5[[i]]$O2.correct )
  spl <- smooth.spline(MM ~ DD , df = 10)
  plot ( MM~DD, ylab=i, cex = 0.4, xaxp = c(0, 900, 10))
  lines(spl, col="red" , lwd=3)
  T1 = predict(spl, x=1, deriv=1)$y
  T2 = predict(spl, x=20, deriv=1)$y
  c(T1,T2)
}))

plot(MMR_HF_data6[,1]) ## slope b. all look reasonable. 
plot(MMR_HF_data6[,2]) 
unlist ( lapply(HF_data5, nrow)  )  ## number of data points per slope

# Compare slopes from T1 and T2.
plot(MMR_HF_data6[,1],MMR_HF_data6[,2])
abline(0,1) 
#Almost identical.

## Exclude poor data identified earlier.
MMR_HF_data6 [c(23,126),] = NA

# Make a new data frame with slopes and original cols (this only takes the first row for each fish,
# i.e. the BR data will not be relevant). Slopes from T1 and T2 in 1 and 2.
MMR_HF_data7 = cbind ( do.call ( rbind , lapply(HF_data5 ,function(x) {x[1,]} ) ) , MMR_HF_data6 )
str(MMR_HF_data7)
head(MMR_HF_data7)
names(MMR_HF_data7)
colnames(MMR_HF_data7)[c(18:19)] = c("T1slp","T20slp")
head(MMR_HF_data7)                  

# Add index to each fish
MMR_HF_data7$index = 1:nrow(MMR_HF_data7)

## Continue with FishResp for MMR calc
#save(MMR_HF_data7, file = "MMR_HF_data7" )
 
# Change format into Fishresp slope format:

# Cols that need to be in the data (empty or not) "Chamber.No"    "Ind"           "Mass"          "Volume"        "Date.Time"     "Phase"  "Temp"          "Slope.with.BR" "Slope"         "SE"            "R2"            "DO.unit"     

names(MMR_HF_data7)

MMR_HF_T1 <- MMR_HF_data7 %>%
  dplyr::select("Chamber.No",   "Ind","Mass","Volume",
         "Date.Time","Phase", "Temp","T1slp", 
        "index",  "DO.unit") %>%
  mutate(SE = NA,
         R2 = NA,
         BR = NA,
         Slope.with.BR =NA) %>%
  dplyr::rename("Slope" = T1slp) %>%
  mutate(Ind = as.factor(Ind))

head(MMR_HF_T1)


### Add batch for each fish to identify problematic chambers from notes. 
batch <- read.table(file = "Data/All_MMR_batch_info.txt", stringsAsFactors = T , h=T, dec = ".", sep = "\t")
head(batch)

batchonly <- dplyr::select(batch, c(PIT_MMR,Batch))

MMR_HF_T1 <- dplyr::inner_join(MMR_HF_T1, batchonly, by = c("Ind" = "PIT_MMR"))
head(MMR_HF_T1)
str(MMR_HF_T1)

#small bubble in chambers 11 and 16
a <- dplyr::filter(batch, Batch  == "C3" & Chamber == 11)$PIT_MMR
b <- dplyr::filter(batch, Batch  == "C4" & Chamber == 16)$PIT_MMR

MMR_HF_T1$Volume <- ifelse(
  MMR_HF_T1$Ind == a | MMR_HF_T1$Ind == b, MMR_HF_T1$Volume -0.03, MMR_HF_T1$Volume )

rm(a)
rm(b)
subset(MMR_HF_T1, Batch  == "C3")$Volume

## Remove two individuals with very large bubbles:
a <- dplyr::filter(batch, Batch  == "C3" & Chamber == 9)$PIT
b <- dplyr::filter(batch, Batch  == "C4" & Chamber == 2)$PIT

MMR_HF_T1 <- dplyr::filter(MMR_HF_T1, !(Ind ==a | Ind ==b)) %>%
  droplevels()

rm(a)
rm(b)

MMR_HF_T1$Volume

# Calculate MMR from slope
MRdat <- calculate.MR(MMR_HF_T1, density=1000, 
                      plot.BR = F,
                      plot.MR.abs = F,
                      plot.MR.mass = F)
head(MRdat)

#Remove NA cols
MRdat <- MRdat %>%
  dplyr::select(-c(MR.abs.with.BR, Slope.with.BR, BR, SE)) %>%
  droplevels()

## Add genotype information
Fish_genot <- read.table(file = "Data/All_family_info.txt", header=T, dec = "\t")
str(Fish_genot)

MRdat <- dplyr::left_join(MRdat, dplyr::select(Fish_genot, c(PIT, "Call_six6", "Call_vgll3", "Sex","Family")),
                                 by = c("Ind" = "PIT"))

head(MRdat)

hist(MRdat$MR.mass)

# Drop heterozygotes
MRdat<- MRdat %>% 
  dplyr::filter(Call_vgll3 %in% c("E", "L") & Call_six6 %in% c("E", "L")) %>%
  droplevels()

str(MRdat)
hist(MRdat$MR.mass)

### Combine with low food and SMR data in another script

MMR_HF <- MRdat
save(MMR_HF, file="Data/MMR_HF")



