#######################################
## SMR analysis for salmon Meas2 SMR data 
## Jenni Prokkola Dec 2020
#######################################
library(data.table)
library(tidyverse)
library(FishResp)
library(chron)
library(ggplot2)
library(lubridate)
library(mclust)

# Data will be processed through filters to select slopes that are 
# decreasing linearly for SMR calculation. See separate scripts for filter functions.
# Outputs are in lists, one quadrant from one batch per element (in total 4 quadrants per batch).

# Import background-corrected raw data 
# Each measurement in separate row
# The correct O2 data to use is in O2.correct column. PIT tag number is in Ind.
LF_data <-read.table ("Data/Lowfood_bgCorrected_alldata_SMR.txt" , stringsAsFactors = T , h=T, dec = ".", sep = "\t")
str(LF_data)

# Add batch and quadrant information to data 
info.batch <- read.table(file="Data/All_batch_info.txt", h=T, dec = ".", sep = "\t")
str(info.batch)

# Select columns and levels
info.batch <- info.batch %>%
  select(PIT, Batch, Quadrant) %>%
  filter(Batch %in% c(paste("B", c(13:24), sep = ""))) %>%
  droplevels()

# Combine
LF_data <- left_join(LF_data, info.batch, by = c("Ind" = "PIT"))
head(LF_data)

# Make new column combining batch and quadrant
LF_data$bq <- paste(LF_data$Batch, LF_data$Quadrant, sep = ".")
head(LF_data)
tail(LF_data)
str(LF_data)
## Make a vector of batches and quadrants, these will be list names
bqs <- unique(LF_data$bq)
str(bqs)

##################################################
##### Function 1. r-sq. value and outlier-filtering
##################################################
# Outputs a list of data in f1.output and a list of raw data in datalist, 
#  a list of r2 value for each slope in r2values, and all excluded slopes in excl.slop.f1
# Run function on data, bqs vector makes element names
rsq.outliers.f(LF_data, bqs)

str(datalist[[1]])
names(datalist)
str(f1.output[[1]])
names(f1.output[[1]])
head(f1.output[[1]][,10:20])
names(f1.output)

#remove large data file
rm(LF_data)

#  Summarise how many slopes have poor r2
lapply(r2values, function(x) x%>% 
         filter(rsq <0.95) %>%
         group_by(Chamber.No) %>%
         summarize(count = length(Phase)))

# Counts of DO values excluded as outliers are saved in NAcounts. Summarise the number of outliers.
lapply(NAcounts, function(x) filter(x, count.o >0))
# Check that the values of specific slope are 0
filter(f1.output[["B22.Q1"]], Chamber.No =="CH4", Phase =="M2")$O2.correct

# Summary of deleted slopes from both r2 and outlier filters is in excl.slop.f1
excl.slop.f1

##################################################
## Function 2. slope linearity filtering
##################################################
# Automatically returns the number of deleted slopes based on linearity
# Produces a list of data in f2.output
slope.lin.f(f1.output)

#sanity check
summary(f2.output[[1]]$O2.correct)

##################################################
## Function 3. Identify and exclude individuals with 
## too many missing slopes after function 1 and 2
##################################################
# The function returns all individuals with a 1 for excluded ones
# Produces a list of data in filt.final.dat
NAslope.f(f2.output)

# Counts of removed and total slopes for each fish 
#check that removed ind matches function output 
# only excluded fish:
rbindlist(lapply(NAslope.counts, function(x) x %>%
                   filter(ToomanyNA == "yes")), idcol = "Batch") 

##################################################
# Make plots of slopes before and after filtering
# linear regression line shown in orange
# note that individuals may have different y-axis scales
##################################################
# raw data, separate pages for individuals (this is pretty slow, go get coffee)
pdf(file="LFSlopePlotRawdata.pdf", width = 14, height = 9)
rawdat.multiplot(datalist)
dev.off()

#filtered data
pdf(file="LFSlopePlotFiltered.pdf", width = 14, height = 9)
rawdat.multiplot(filt.final.dat)
dev.off()

##################################################
#### Calculate slopes and oxygen consumption using Fishresp 
##################################################
## Droplevels is essential for this to run!!
slope.list <- lapply(filt.final.dat, function(x) extract.slope(droplevels(x), method="all", length = 720))
head(slope.list[[1]])

# Change the volume of chambers with bubbles
# B24 Ch15 reduced by 30ul = 0.03 ml
fix.vol <-slope.list[["B24.Q4"]]
fix.vol$Volume<- ifelse(fix.vol$Chamber.No == "CH15",
    (fix.vol$Volume-0.03), fix.vol$Volume)

slope.list[["B24.Q4"]] <- fix.vol
rm(fix.vol)

# B15 Ch2 reduced by 60ul = 0.06 ml
fix.vol <-slope.list[["B15.Q1"]]
fix.vol$Volume<- ifelse(fix.vol$Chamber.No == "CH2",
                        (fix.vol$Volume-0.06), fix.vol$Volume)

slope.list[["B15.Q1"]] <- fix.vol
rm(fix.vol)

slope.list[["B15.Q1"]]$Volume

# B16 Ch9 reduced by 30ul = 0.03 ml
fix.vol <-slope.list[["B16.Q3"]]
fix.vol$Volume<- ifelse(fix.vol$Chamber.No == "CH9",
                        (fix.vol$Volume-0.03), fix.vol$Volume)

slope.list[["B16.Q3"]] <- fix.vol
rm(fix.vol)

slope.list[["B16.Q3"]]$Volume

# B18 Ch8 reduced by 30ul = 0.03 ml
fix.vol <-slope.list[["B18.Q2"]]
fix.vol$Volume<- ifelse(fix.vol$Chamber.No == "CH8",
                        (fix.vol$Volume-0.03), fix.vol$Volume)

slope.list[["B18.Q2"]] <- fix.vol
rm(fix.vol)

slope.list[["B18.Q2"]]$Volume


# Calculate oxygen consumption in mg/h for each slope (if plots are needed, run on each element
# of slope.list separately)
#calculate.MR(slope.list[[1]], density=1000, 
#             plot.BR = T,
#             plot.MR.abs = F,
#             plot.MR.mass = T)


MR.list<-lapply(slope.list, function(x) calculate.MR(x, density=1000, 
                                                     plot.BR = F,
                                                     plot.MR.abs = F,
                                                     plot.MR.mass = F))
str(MR.list[[1]])

# Plot MO2 over day for each fish
pdf(file="All_batches_MR_plot_LF.pdf", width = 8, height = 9)
MR.dat.plot(MR.list)
dev.off()

# Determine if some fish or phases need to be excluded based on the plot, e.g. missing long period of slopes or otherwise (errors etc)
rem <- as.character(unique(subset(filt.final.dat[["B13.Q4"]], Chamber.No == "CH13")$Ind))
rem <- c(rem, as.character(unique(subset(filt.final.dat[["B19.Q3"]], Chamber.No == "CH11")$Ind)))
rem <- c(rem, as.character(unique(subset(filt.final.dat[["B13.Q4"]], Chamber.No == "CH13")$Ind)))
rem <- c(rem, as.character(unique(subset(filt.final.dat[["B21.Q3"]], Chamber.No == "CH12")$Ind)))
rem <- c(rem, as.character(unique(subset(filt.final.dat[["B23.Q4"]], Chamber.No == "CH13")$Ind)))

MR.list.SMR <- lapply(MR.list, function(x) filter(x, !(Ind %in% rem)))
unique(MR.list.SMR[["B13.Q4"]]$Chamber.No)
unique(MR.list.SMR[["B21.Q3"]]$Chamber.No)

################################################################
## Use calcSMR function by Chabot et al. 2016 to extract SMR variables
################################################################
# Check the supplements of Chabot et al.The determination of standard metabolic rate in fishes. Journal of Fish Biology, 88, 81-121, for the calcSMR R function.

SMR.result <- lapply(MR.list.SMR, function(x) x %>%
                       group_by(Chamber.No) %>%
                       summarise(
                         SMR.mass.mlnd    = calcSMR(MR.mass)$mlnd,
                         SMR.abs.mlnd = calcSMR(MR.abs)$mlnd,
                         SMR.mass.q0.1 = calcSMR(MR.mass)$quant[1],
                         SMR.abs.q0.1 = calcSMR(MR.abs)$quant[1],
                         SMR.mass.q0.2 = calcSMR(MR.mass)$quant[3],
                         SMR.abs.q0.2 = calcSMR(MR.abs)$quant[3],
                         PIT    = unique(Ind),
                         Mass   = unique(Mass),
                         Temp = mean(Temp),
                         # Date time needs to be parsed to show only the start date of measurement
                         Date = unique(floor_date(ymd_hms(Date.Time), unit = "day"))[1]
                       ))


# Plot mlnd method for finding SMR. 
# Input MR.list with all phases, 
# SMR.result for showing SMR in the plot and 
# number of individuals to plot from data

pdf(file="Mcclust_plot_LF.pdf", width = 9, height = 10)
par(mfcol = c(2,2))
plot.mcclust.f(MR.list.SMR, SMR.result, 190)
dev.off()

# Plot mass-SMR relationship in different metrics (total fish-specific SMR, i.e. SMR.abs)
vars <- c("SMR.abs.mlnd", "SMR.abs.q0.2", "SMR.abs.q0.1")
SMR.weight.plot(SMR.result, vars)

##################################################
## Bind SMR data to fish genotype, sex and family data
##################################################
# Bind results into a dataframe and add Batch column
# idcol makes an additional column in the data from list names, keep only Batch column
result.SMR <- rbindlist(SMR.result, idcol = "Batch.Q") %>%
  separate(Batch.Q, c("Batch", NA))

# Read genotypes, sex, and family
Fish.genot <- read.table(file="Data/All_family_info.txt", sep = "\t", header = T)
str(Fish.genot)

LF.SMR <- inner_join(Fish.genot, result.SMR, by = "PIT")
head(LF.SMR)
str(LF.SMR)

# remove large lists
rm(f1.output)
rm(f2.output)

## Exclude B13 data from LF, because had a shorter acclimation time (more food left) 
LF.SMR <- LF.SMR %>%
  filter(!(Batch == "B13")) %>%
  droplevels()

# Exclude H-genotypes
LF.SMR<- LF.SMR %>% 
  filter(Call_vgll3 %in% c("E", "L") & Call_six6 %in% c("E", "L")) %>%
  droplevels()
str(LF.SMR)

# Export object to use in statistical and aerobic scope analysis:
save(LF.SMR, file ="Data/LF.SMR.data")

##################################################
### Exploratory data plotting 
##################################################
# Different SMR variables
plot(LF.SMR$SMR.mass.mlnd~LF.SMR$SMR.mass.q0.2)
abline(0,1)
cor(LF.SMR$SMR.mass.mlnd,LF.SMR$SMR.mass.q0.2)

plot(LF.SMR$SMR.mass.mlnd~LF.SMR$SMR.mass.q0.1)
abline(0,1)
cor(LF.SMR$SMR.mass.mlnd,LF.SMR$SMR.mass.q0.1)

## Mass-specific by batch (fish get bigger towards the end)
ggplot(LF.SMR, aes(y= SMR.mass.mlnd, x=Batch))+
  geom_boxplot() +
  geom_jitter(alpha =0.5)+
  ylab("SMR mg O2/kg/h")+
  theme_bw()

## Mass-specific by chamber 
ggplot(LF.SMR, aes(y= SMR.mass.mlnd, x=Chamber.No))+
  geom_boxplot() +
  geom_jitter(alpha =0.5)+
  ylab("SMR mg O2/kg/h")+
  theme_bw()

## Mass-specific by family
ggplot(LF.SMR, aes(y= SMR.mass.mlnd, x=Family))+
  geom_boxplot() +
  geom_jitter(alpha =0.5)+
  ylab("SMR mg O2/kg/h")+
  theme_bw()

# Clear workspace and Continue in Combine_SMR and then SMR_models scripts
