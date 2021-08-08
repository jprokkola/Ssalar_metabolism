######################################################################
## MMR analysis with 1 min windows to find the highest MO2 with a linear slope
## Low food treatment
######################################################################
library(respR)
library(plyr)
library(tidyverse)
library(data.table)
library(FishResp)
library(ggplot2)

# Read background corrected data, only relevant points for each slope (one slope per fish)
# each fish in separate element in a list
load("Data/MMR_LF_data5") # A list of all O2 data for the MMR slope for each fish
head(MMR_LF_data5[[1]])

# name the elements
namevec <- c(ldply(MMR_LF_data5, function (f)  paste(f$Ind[1], f$Phase[1]) ))
MMR_LF_data5<- setNames(MMR_LF_data5, namevec$V1)
names(MMR_LF_data5)

# Make a function for using auto-rate in each fish, 1 minute window (30 separate measurements)

autorate_all_f <- function(data) {
  mmrlist <<- llply(data, function(x) { # result for each fish and phase in a separate element
  d <- x %>%
        select(Time, O2.correct, Phase,Chamber.No, Ind) 
      # autorate takes the first two columns as default (Time and O2)
      auto_rate(d, method = "max", width = 60) #1-min window (30 data points)
    })
}


## Run auto_rate function (prints plots to device)
autorate_all_f(MMR_LF_data5)

head(mmrlist)

# To save plots into pdf, run again.
pdf(file= "LF.autorate.plots.1min.TAdat.pdf", width = 11, height = 11)
autorate_all_f(MMR_LF_data5)
dev.off()

summary(mmrlist[[1]])

## Data are in a list, each object is in format auto_resp
mmrlist[[1]]
# The max slope is the first slope in the output in column Rate
# Sometimes the max slope has a low R2. First filter all data by R2, then pick the first slope left.
# Note that the plots are only made for the best slope, 
# so the selected slopes need to be plotted again.
# Take only slopes and rsq from the auto-rate objects
# mmrlist[[1]]$summary$rsq #example
# mmrlist[[1]]$rate #example

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

## Separate ID into two cols
max.slop.df <- max.slop.df %>%
  separate(ID, into = c("PIT", "Phase"), sep =" ")

## Combine to data 
MMR_LF_data5_df <- rbindlist(MMR_LF_data5)
max.slop.meas <- right_join(MMR_LF_data5_df, max.slop.df, by = c("Ind"= "PIT", "Phase" = "Phase"))

# Use start row and end row to highlight the selected section
pdf(file = "MMR_slopes_LF_1min_TAdat.pdf", width = 8, height = 8)

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
    ylab("O2 mg / L")+
    xlab("Time (s)")+
    ggtitle(i)+
    theme_bw()+
    theme(axis.text.x = element_text(size = 8),
          axis.title.x = element_blank())
  print(p)
  rm(plotdat)
  rm(start)
  rm(end)
  rm(slope)
  rm(i)
}

dev.off()

## Add relevant information from the input file.

mmr_info <- MMR_LF_data5_df %>% 
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

# Exclude chambers with problems (big air bubble etc, same ind as TA analysis)
max.slop.df.full <- max.slop.df.full %>%
  filter(!(Ind %in% c("A00000D0900001897005743", 
                      "A00000D0900001897005830", "A00000D0900001897007624", 
                      "A00000D0900001897005837", "A00000D0900001897005987",
                      "A00000D0900001897007276")))

## Calculate MR in mg/O2/h for each phase
max.mr.df <- calculate.MR(max.slop.df.full, density=1000, 
                          plot.BR = F,
                          plot.MR.abs = F,
                          plot.MR.mass = F)


# Exclude chambers with other problems (poor data quality in the beginning or whole slope)

max.mr.df <- max.mr.df %>%
  filter(!(Ind %in% c("A00000D0900001897005401", "A00000D0900001897005837",
                      "A00000D0900001897005987", "A00000D0900001897007772")))


plot(max.mr.df$MR.mass) #sanity check

### Add batch and family information 
batch <- select(batch, -c(weight, Quadrant, Channel, Chamber, Date))
max.mr.df.genot <- left_join(max.mr.df, batch, by = c("Ind" = "PIT"))
head(max.mr.df.genot)

## Read genotypes + sex
Fish.genot <- read.table(file="Data/All_family_info.txt", sep = "\t", header = T)
str(Fish.genot)

Fish.genot <- select(Fish.genot, c(PIT, Call_vgll3, Call_six6, Sex))

max.mr.df.genot <- left_join(max.mr.df.genot, Fish.genot, by = c("Ind" = "PIT"))
head(max.mr.df.genot)
str(max.mr.df.genot)
max.mr.df.genot <- droplevels(max.mr.df.genot) #Drop extra batches

## Plot by batch
boxplot(MR.mass ~ Batch, data = max.mr.df.genot)

# Drop hets and B13
max.mr.df.genot <- max.mr.df.genot %>% 
  filter(!(Batch == "B13")) %>%
  filter(Call_vgll3 %in% c("E", "L") & Call_six6 %in% c("E", "L")) %>%
  select(-c(SE, MR.abs.with.BR, BR))%>%
  droplevels()

MMR_respr_LF <- max.mr.df.genot

save(MMR_respr_LF, file = "MMR_respr_LF")


