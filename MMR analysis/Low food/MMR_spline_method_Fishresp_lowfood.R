######################################################################
## MMR calculation from slopes using FishResp
## Slopes from Tutku's analysis using tangent from the beginning of 
## a quadratic slope
######################################################################
library(FishResp)
library(tidyverse)
library(ggplot2)
library(lubridate)

load("Data/MMR_LF_data5") # A list of all O2 data for the MMR slope for each fish
head(MMR_LF_data5[[1]])

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
pdf(file = "MMR_LF_slopes_splines.pdf", width = 5, height = 5)
rawdat.multiplot(MMR_LF_data5)
dev.off()


### Spline method to get tangents at 1 and 20 seconds (T1 and T2). This also plots all slopes.
## Df means the number of turns possible on the slope
MMR_LF_data6 = t(sapply ( 1: length(MMR_LF_data5) , function(i) {
  DD = 1:length ( MMR_LF_data5[[i]]$O2.correct )
  MM = ( MMR_LF_data5[[i]]$O2.correct )
  spl <- smooth.spline(MM ~ DD , df = 10 )
  plot ( MM~DD, ylab=i, cex = 0.4, xaxp = c(0, 900, 10))
  lines(spl, col="red" , lwd=3)
  T1 = predict(spl, x=1, deriv=1)$y
  T2 = predict(spl, x=20, deriv=1)$y
  c(T1,T2)
}))

plot.new()
plot(MMR_LF_data6[,1])
plot(MMR_LF_data6[,2])
plot(MMR_LF_data6[,1],MMR_LF_data6[,2])
abline(0,1) # 1 outlier (T1 and T2 don't agree)

MMR_LF_data6[which(MMR_LF_data6[,2] < -0.008),] = NA
head(MMR_LF_data6)

unlist ( lapply(MMR_LF_data5, nrow)  )  ## number of data points per slope

# Make a new data frame with slopes and original cols (this only takes the first row for each fish,
# i.e. the BR data will not be relevant). Slopes from T1 and T2 in 1 and 2.
MMR_LF_data7 = cbind ( do.call ( rbind , lapply(MMR_LF_data5 ,function(x) {x[1,]} ) ) , MMR_LF_data6 )
str(MMR_LF_data7)
head(MMR_LF_data7)
names(MMR_LF_data7)
colnames(MMR_LF_data7)[c(18:19)] = c("T1slp","T20slp")
head(MMR_LF_data7)                  

# Add index to each fish
MMR_LF_data7$index = 1:nrow(MMR_LF_data7)

## Continue with FishResp for MMR calc
#save(MMR_LF_data7, file = "MMR_LF_data7" )

# Change format into Fishresp slope format:

# Cols that need to be in the data (empty or not) "Chamber.No"    "Ind"           "Mass"          "Volume"        "Date.Time"     "Phase"
# "Temp"          "Slope.with.BR" "Slope"         "SE"            "R2"            "DO.unit"     

names(MMR_LF_data7)

MMR_LF_T1 <- MMR_LF_data7 %>%
  dplyr::select("Chamber.No",   "Ind","Mass","Volume",
                "Date.Time","Phase", "Temp","T1slp", 
                "index",  "DO.unit") %>%
  mutate(SE = NA,
         R2 = NA,
         BR = NA,
         Slope.with.BR =NA) %>%
  dplyr::rename("Slope" = T1slp)


head(MMR_LF_T1)

### Add batch for each fish to identify problematic chambers from notes. Note that chambers in batch info refer to SMR, not MMR, but batches are the same for both.
batch <- read.table(file = "Data/All_MMR_batch_info.txt", stringsAsFactors = T , h=T, dec = ".", sep = "\t")
head(batch)

batchonly <- dplyr::select(batch, c(PIT_MMR,Batch))

MMR_LF_T1 <- dplyr::inner_join(MMR_LF_T1, batchonly, by = c("Ind" = "PIT_MMR"))
head(MMR_LF_T1)

MRdat <- calculate.MR(MMR_LF_T1, density=1000, 
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

## Exclude B13 data (not enough time in acclimation)
MRdat <- MRdat %>%
  filter(!(Batch == "B13"))
head(MRdat)

# Exclude chambers with other problems (big air bubble etc)
subset(MRdat, Batch == "B15" & Chamber.No == "CH7")
subset(MRdat, Batch == "B15" & Chamber.No == "CH8")
subset(MRdat, Batch == "B24" & Chamber.No == "CH13")
subset(MRdat, Batch == "B14" & Chamber.No == "CH15")
subset(MRdat, Batch == "B17" & Chamber.No == "CH14")


MRdat <- MRdat %>%
  filter(!(Ind %in% c("A00000D0900001897005743", "A00000D0900001897005830", "A00000D0900001897007624", "A00000D0900001897005837", "A00000D0900001897005987")))

# Drop heterozygotes
MRdat<- MRdat %>% 
  filter(Call_vgll3 %in% c("E", "L") & Call_six6 %in% c("E", "L")) %>%
  droplevels()
#163 fish left

str(MRdat)

hist(MRdat$MR.abs)
MMR_LF <- MRdat

### Combine with high food and SMR data in another script
save(MMR_LF, file="Data/MMR_LF")
