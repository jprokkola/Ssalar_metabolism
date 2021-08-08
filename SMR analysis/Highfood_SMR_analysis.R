###################################################
## SMR analysis for juvenile salmon high food group 
##  Prokkola et al. 2021
###################################################
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
HF_data <-read.table ("Data/Highfood_bgCorrected_alldata_SMR.txt" , stringsAsFactors = T , h=T, dec = ".", sep = "\t")
str(HF_data)

# Add batch and quadrant information to data 
info.batch <- read.table(file="Data/All_batch_info.txt", h=T, dec = ".", sep = "\t")
range(info.batch$weight)
str(info.batch)

# Select columns and levels
info.batch <- info.batch %>%
  select(PIT, Batch, Quadrant) %>%
  filter(Batch %in% c(paste("C", c(1:8), sep = ""))) %>%
  droplevels()

# Combine
HF_data <- left_join(HF_data, info.batch, by = c("Ind" = "PIT"))
head(HF_data)

# Make new column combining batch and quadrant
HF_data$bq <- paste(HF_data$Batch, HF_data$Quadrant, sep = ".")
head(HF_data)
tail(HF_data)
str(HF_data)

## Make a vector of batches and quadrants, these will be list names
bqs <- unique(HF_data$bq)
str(bqs)

##################################################
##### Function 1. r-sq. value and outlier-filtering
##################################################
# Outputs a list of data in f1.output and a list of raw data in datalist, 
#  a list of r2 value for each slope in r2values, and all excluded slopes in excl.slop.f1
# Run function on data, bqs vector makes element names
rsq.outliers.f(HF_data, bqs)

#sanity checks
str(datalist[[1]])
names(datalist)
str(f1.output[[1]])
names(f1.output[[1]])
head(f1.output[[1]][,10:20]) # When whole slope excluded, O2.correct = 0
names(f1.output)

#remove large data file
rm(HF_data)

#  Summarise how many slopes have poor r2
lapply(r2values, function(x) x%>% 
         filter(rsq <0.95) %>%
         group_by(Chamber.No) %>%
         summarize(count = length(Phase)))

# Counts of DO values excluded as outliers are saved in NAcounts. Summarise the number of outliers.
lapply(NAcounts, function(x) filter(x, count.o >0))

# Summary of deleted slopes from both r2 and outlier filters is in excl.slop.f1
excl.slop.f1

##################################################
## Function 2. slope linearity filtering
##################################################
# outputs lists linearity.data, and f2.output, and 
# returns the number of excluded slopes for each fish
slope.lin.f(f1.output)

#sanity check
summary(f2.output[[1]]$O2.correct)
head(linearity.data)
lapply(linearity.data, function (x) filter(x, SLOPCV >15))

##################################################
## Function 3. Identify and exclude individuals with 
## too many missing slopes after function 1 and 2
##################################################
# The function returns all individuals with a 1 for excluded ones
# Produces a list of data in filt.final.dat
NAslope.f(f2.output)

# Counts of removed and total slopes for each fish.
# Check that removed ind matches function output 
# only excluded fish:
rbindlist(lapply(NAslope.counts, function(x) x %>%
  filter(ToomanyNA == "yes")), idcol = "Batch") 

##################################################
# Make plots of slopes before and after filtering
# linear regression line shown in orange
# note that individuals may have different y-axis scales
##################################################
# raw data, separate pages for individuals (this is pretty slow, go get coffee)
pdf(file="HFSlopePlotRawdata.pdf", width = 14, height = 9)
rawdat.multiplot(datalist)
dev.off()

#filtered data
pdf(file="HFSlopePlotFiltered.pdf", width = 14, height = 9)
rawdat.multiplot(filt.final.dat)
dev.off()

#######################################################################
#### Calculate slopes and oxygen consumption using Fishresp and Mcclust
#######################################################################
## Droplevels is essential for this to run!!
# Note that this version of Fishresp (1.0.4) has been fixed not to complain about NA data points.
# This takes a while:
slope.list <- lapply(filt.final.dat, function(x) 
  extract.slope(droplevels(x), method="all", length = 720))

head(slope.list[[1]])

# Reduce volume in the chambers with air bubble 
# BC4 Ch9 reduced by 30ul = 0.03 ml
fix.vol <-slope.list[["C4.Q3"]]
fix.vol$Volume<- ifelse(fix.vol$Chamber.No == "CH9",
                        (fix.vol$Volume-0.03), fix.vol$Volume)

slope.list[["C4.Q3"]] <- fix.vol
rm(fix.vol)

slope.list[["C4.Q3"]]$Volume

# BC8 Ch2 and 3 reduced by 60ul = 0.06 ml
fix.vol <-slope.list[["C8.Q1"]]
fix.vol$Volume<- ifelse(fix.vol$Chamber.No == "CH2" | fix.vol$Chamber.No == "CH3",
                        (fix.vol$Volume-0.06), fix.vol$Volume)

slope.list[["C8.Q1"]] <- fix.vol
rm(fix.vol)

slope.list[["C8.Q1"]]$Volume

# Calculate oxygen consumption in mg/h for each slope (if plots are needed, run on each element
# of slope.list separately)
MR.list<-lapply(slope.list, function(x) calculate.MR(x, density=1000, 
                                                     plot.BR = F,
                                                     plot.MR.abs = F,
                                                     plot.MR.mass = F))
str(MR.list[[1]])

# Plot MO2 over day for each fish
pdf(file="All_batches_MR_plot_HF.pdf", width = 8, height = 9)
MR.dat.plot(MR.list)
dev.off()

# Determine if some fish or phases need to be excluded based on the plot or otherwise
rem <- as.character(unique(subset(filt.final.dat[["C8.Q4"]], Chamber.No == "CH13")$Ind)) #due to big bubble
rem <- c(rem, as.character(unique(subset(filt.final.dat[["C3.Q1"]], Chamber.No == "CH2")$Ind))) #strangely high baseline, not a large or small fish

# Exclude these:
MR.list.SMR <- lapply(MR.list, function(x) filter(x, !(Ind %in% rem)))
unique(MR.list.SMR[["C8.Q4"]]$Chamber.No)

# Make a table from the list if needed
#MR.dat.HF.allslopes<- rbindlist(MR.list.SMR, idcol = "bq")
#write.table(MR.dat.HF.allslopes, file = "MR.dat.HF.allslopes.txt", dec = ".", 
#            quote = F, row.names = F, sep = "\t")

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

pdf(file="Mcclust_plot_HF.pdf", width = 9, height = 10)
par(mfcol = c(2,2))
plot.mcclust.f(MR.list, SMR.result, 190)
dev.off()


# Plot mass-SMR relationship in different metrics (total fish-specific SMR, i.e. SMR.abs)
vars <- c("SMR.abs.mlnd", "SMR.abs.q0.2", "SMR.abs.q0.1")
SMR.weight.plot(SMR.result, vars)

##################################################
## Bind SMR data to fish genotype, sex and family data
##################################################
# Bind results into a data frame and add Batch column
# idcol makes an additional column in the data from list names, keep only Batch column
result.SMR <- rbindlist(SMR.result, idcol = "Batch.Q") %>%
  separate(Batch.Q, c("Batch", NA))

# Read genotypes, sex, and family
Fish.genot <- read.table(file="Data/All_family_info.txt", sep = "\t", header = T)
str(Fish.genot)

HF.SMR <- inner_join(Fish.genot, result.SMR, by = "PIT")
head(HF.SMR)
str(HF.SMR)

# remove large lists
rm(Fish.genot)
rm(f1.output)
rm(f2.output)
rm(MR.list)
rm(filt.final.dat)
rm(datalist)

# Exclude H-genotypes
HF.SMR<- HF.SMR %>% 
  filter(Call_vgll3 %in% c("E", "L") & Call_six6 %in% c("E", "L")) %>%
  droplevels()
str(HF.SMR)

# Export object to use in statistical and aerobic scope analysis:
save(HF.SMR, file ="Data/HF.SMR.data")

##################################################
### Explorative data plotting 
##################################################
# Different SMR variables
plot(result.SMR$SMR.mass.mlnd~result.SMR$SMR.mass.q0.2)
abline(0,1)
cor(result.SMR$SMR.mass.mlnd,result.SMR$SMR.mass.q0.2)

plot(result.SMR$SMR.mass.mlnd~result.SMR$SMR.mass.q0.1)
abline(0,1)
cor(result.SMR$SMR.mass.mlnd,result.SMR$SMR.mass.q0.1)

## Mass-specific by batch (fish get bigger towards the end)
ggplot(HF.SMR, aes(y= SMR.mass.mlnd, x=Batch))+
  geom_boxplot() +
  geom_jitter(alpha =0.5)+
  ylab("SMR mg O2/kg/h")+
  theme_bw()

## Mass-specific by chamber 
ggplot(HF.SMR, aes(y= SMR.mass.mlnd, x=Chamber.No))+
  geom_boxplot() +
  geom_jitter(alpha =0.5)+
  ylab("SMR mg O2/kg/h")+
  theme_bw()

## Mass-specific by family
ggplot(HF.SMR, aes(y= SMR.mass.mlnd, x=Family))+
  geom_boxplot() +
  geom_jitter(alpha =0.5)+
  ylab("SMR mg O2/kg/h")+
  theme_bw()


# Clear workspace, analyse Low food data, and Continue in Combine_SMR and then SMR_models scripts

##################################################
## Plotting example slopes
##################################################

# Comparison of linear and non-linear slopes (CV <15 or >15)
# Non-linear (CV 15.4) C6 Ch10 M22
# Extract rows from raw data and print new figures
ch10.lin <-filter(linearity.data[["C6.Q3"]], Chamber.No =="CH10" & Phase == "M20" |Chamber.No =="CH10" & Phase == "M22")

ch10.dat <- filter(HF, bq == "C6.Q3" & Chamber.No =="CH10" & Phase == "M20" |bq == "C6.Q3" & Chamber.No =="CH10" & Phase == "M22")

ggplot(data = ch10.dat, 
       aes(x = Time, y = O2.correct ))+
  facet_wrap(~Phase)+
  geom_point(cex = 0.4)+
  geom_line(stat="smooth",method = lm, colour = "darkorange", na.rm = T)+
  ylab(expression(paste("mg ", O[2])))+
  #ggtitle(paste(names(data)[[x]], c))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 8))
