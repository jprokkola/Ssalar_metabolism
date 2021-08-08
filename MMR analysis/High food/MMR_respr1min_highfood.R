############################################################################
## MMR analysis with 1 min sliding windows to find the highest MO2 from a linear slope
## High food treatment
############################################################################
library(respR)
library(plyr)
library(tidyverse)
library(data.table)
library(FishResp)
library(ggplot2)

# Read background corrected data, only relevant points for each slope (one slope per fish)
# Data prepared with TAmethod script.
# Each fish in separate element in a list
load("Data/MMR_HF_data5")

str(HF_data5)
head(HF_data5[[1]]) 

# name the elements
namevec <- c(ldply(HF_data5, function (f)  paste(f$Ind[1], f$Phase[1]) ))
HF_data5<- setNames(HF_data5, namevec$V1)
names(HF_data5)

# Make a function for using auto-rate in each fish, 1 minute window (30 separate measurements)

autorate_all_f <- function(data) {
  mmrlist <<- llply(data, function(x) { # result for each fish and phase in a separate element
  d <- x %>%
        select(Time, O2.correct, Phase,Chamber.No, Ind) 
      # autorate takes the first two columns as default (Time and O2)
      auto_rate(d, method = "max", width = 60) #1-min window (30 data points)
    })
}


## Run respR auto_rate function (prints plots to device)
autorate_all_f(HF_data5)

head(mmrlist)

# To save plots into pdf, run again.
pdf(file= "autorate.plots.1min.TAdat.pdf", width = 11, height = 11)
autorate_all_f(HF_data5)
dev.off()

summary(mmrlist[[1]])

## Data are in a list, each object is in format auto_resp
mmrlist[[1]]
# The max slope is the first slope in the output in column Rate
# Sometimes the max slope has a low R2. First filter all data by R2, then pick the first remaining slope.
# Note that the plots are only made for the best slope, 
# so the selected slopes need to be plotted again.
# mmrlist[[1]]$summary$rsq #example
# mmrlist[[1]]$rate #example

# Take only slopes and rsq from the auto-rate objects
## This only runs with plyr, not dplyr, summarise. Llply preserves element names
mmrlist.slop.rsq <- llply(mmrlist, function(x) x %>%
                            plyr::summarise(max  = rate,
                                            rsq = x$summary$rsq,
                                            startrow = x$summary$row,
                                            endrow = x$summary$endrow))
head(mmrlist.slop.rsq[[1]])

# filter based on rsq 
mmrlist.slop.rsq <- llply(mmrlist.slop.rsq, function(x) x %>%
                            filter(rsq >0.95)) 

# select the top slope and rsq                         
max.slop.list <- llply(mmrlist.slop.rsq, function(x) x %>%
                         summarise(max.slop  = max[1],
                                   rsq = rsq[1],
                                   startrow = startrow[1],
                                   endrow = endrow[1]))
names(max.slop.list)
max.slop.list[[1]]

## Make a table
max.slop.df <- rbindlist(max.slop.list, idcol = "ID")
max.slop.df

# To plot the slopes with the region highlighted in the plot, need to make a df with data and start row and end row
## Separate ID into two cols
max.slop.df <- max.slop.df %>%
  separate(ID, into = c("PIT", "Phase"), sep =" ")

## Combine to data using IDs
HF_data5_df <- rbindlist(HF_data5)
max.slop.meas <- right_join(HF_data5_df, max.slop.df, by = c("Ind"= "PIT", "Phase" = "Phase"))
head(max.slop.meas)

# Use start row and end row to highlight the selected section
pdf(file = "MMR_slopes_HF_1min.pdf", width = 8, height = 8)

for (i in unique(max.slop.meas$Ind)) {
  plotdat = filter(max.slop.meas, Ind == i)
  start = plotdat$startrow[1]
  end = plotdat$endrow[1]
  slope = slice(plotdat, start:end)
  p <- ggplot(data = plotdat, 
              aes(x = Time, y = O2.correct ))+
    geom_point(cex = 0.4)+
    layer(geom = "point", stat = "identity", position = "identity",
          data = slope, show.legend = F,
          params = list(color = "green", cex = 0.4)) + 
    ylab("O2 mg")+
    ggtitle(i)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, size = 8),
          axis.title.x = element_blank())
  print(p)
  rm(plotdat)
  rm(start)
  rm(end)
  rm(slope)
  rm(i)
}

dev.off()
# Check the slopes. Remove bad ind later: A00000D0900001897005852 and A00000D0900001897005195,
# and A00000D0900001897005042, A00000D0900001897007208, A00000D0900001897007035

## Add relevant information from the input file to use in FishResp. 

mmr_info <- HF_data5_df %>% 
  select(c(Date, Chamber.No, Ind, Mass, Volume, Temp)) %>%
  group_by(Ind, Date, Chamber.No, Mass, Volume) %>%
  summarise(Temp = mean(Temp))

head(mmr_info)

# Use right join to include only all rows from y
max.slop.df.full <- right_join(mmr_info, max.slop.df, by = c("Ind" = "PIT"))
head(max.slop.df.full)

rm(mmr_info)
## Input this into Fishresp analysis
# Rename columns to match Fishresp slope file, create empty columns for other variables

names(max.slop.df.full)[8:9] <- c("Slope", "R2")
max.slop.df.full$Slope.with.BR <- 0
max.slop.df.full$SE <- 0

# Adjust volumes for chambers with medium or small air bubbles:
# Read batch and family file to identify individuals
batch <- read.table(file = "Data/All_batch_info.txt", stringsAsFactors = T , h=T, dec = ".", sep = "\t")
head(batch)

#small bubble in chambers 11 and 16
a <- filter(batch, Batch  == "C3" & Chamber == 11)$PIT
b <- filter(batch, Batch  == "C4" & Chamber == 16)$PIT

max.slop.df.full$Volume <- ifelse(
  max.slop.df.full$Ind == a | max.slop.df.full$Ind == b, max.slop.df.full$Volume -0.03, max.slop.df.full$Volume )

rm(a)
rm(b)

## Remove two individuals with very large bubbles:
a <- filter(batch, Batch  == "C3" & Chamber == 9)$PIT
b <- filter(batch, Batch  == "C4" & Chamber == 2)$PIT

max.slop.df.full <- filter(max.slop.df.full, !(Ind ==a | Ind ==b))

rm(a)
rm(b)

## Calculate MR in mg/O2/h for each phase
max.mr.df <- calculate.MR(max.slop.df.full, density=1000, 
                          plot.BR = F,
                          plot.MR.abs = F,
                          plot.MR.mass = F)

# Exclude chambers with other problems (poor data quality in the beginning or whole slope)

max.mr.df <- max.mr.df %>%
  filter(!(Ind %in% c("A00000D0900001897005852", "A00000D0900001897005195",
            "A00000D0900001897005042", "A00000D0900001897007035")))


plot(max.mr.df$MR.mass) #sanity check

### Add batch and family information 
batch <- select(batch, -c(weight, Quadrant, Channel, Chamber, Date))
max.mr.df.genot <- left_join(max.mr.df, batch, by = c("Ind" = "PIT"))
head(max.mr.df.genot)

# Read genotypes + sex
Fish.genot <- read.table(file="Data/All_family_info.txt", sep = "\t", header = T)
str(Fish.genot)

Fish.genot <- select(Fish.genot, c(PIT, Call_vgll3, Call_six6, Sex))

max.mr.df.genot <- left_join(max.mr.df.genot, Fish.genot, by = c("Ind" = "PIT"))
head(max.mr.df.genot)
str(max.mr.df.genot)
max.mr.df.genot <- droplevels(max.mr.df.genot) #Drop extra batches

## Plot by batch
boxplot(MR.mass ~ Batch, data = max.mr.df.genot)

# Drop heterozygotes
max.mr.df.genot <- max.mr.df.genot %>% 
  filter(Call_vgll3 %in% c("E", "L") & Call_six6 %in% c("E", "L")) %>%
  select(-c(SE, MR.abs.with.BR))%>%
  droplevels()

str(max.mr.df.genot)

MMR_respr_HF <- max.mr.df.genot

save(MMR_respr_HF, file = "Data/MMR_respr_HF")
