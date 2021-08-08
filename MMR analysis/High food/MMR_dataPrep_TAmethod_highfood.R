########################################################################
## Identify the beginning of MMR slope, High food treatment
## Tutku Aykanat, Jenni Prokkola (TA method)
########################################################################
library(plyr)
library(tidyverse)
library(ggplot2)

# import high food data. Background respiration -corrected data is in O2.correct.
HF_data = read.table("Data/Highfood_bgCorrected_alldata_MMR.txt" ,
                         h=T, sep = "\t" , stringsAsFactors = F)

head(HF_data)
nrow(HF_data)/900

## Find the the start and end of the first MMR phase. (OK, I know there are easier ways (i.e using "Phase" column but I just saw it late)  )
first_record = t(sapply ( unique(HF_data$Ind) , function(i) { 
  HF_data2 = HF_data [ which ( HF_data$Ind == i ) , ]
  # use the absolute change in O2.correct to identify measurements that have a fish in the chamber
  diff_mmr = HF_data2$O2.correct - c(HF_data2$O2.correct[-1],max(HF_data2$O2.correct) ) 
  END = which ( diff_mmr < -1 )[1]+1
  STR = END-900
  c(STR,END)
  
}))
head(first_record)

# The number of phases for each fish varies 1-5 (entered chambers at different times)
num_phase2 = sapply ( unique(HF_data$Ind) , function(i) { 
  HF_data2 = HF_data [ which ( HF_data$Ind == i ) , ]
  diff_mmr = HF_data2$O2.correct - c(HF_data2$O2.correct[-1],max(HF_data$O2.correct) ) 
  length(which( ( diff_mmr < -1 )))          
})

head(num_phase2) 
range(num_phase2)

table(num_phase2) ## there are two phases that does not have distinct slopes. exclude them.
mmrIDs = names(which ( num_phase2 != 0)) 
length(mmrIDs)# 126 fish have slopes

# Make a list of the first slope for each fish
HF_data4 = (sapply  ( mmrIDs , simplify=F , function(i) {
  pit_ind1 = which ( HF_data$Ind == i )
  pit_ind2 = first_record [ which ( rownames(first_record)==i ) , ]
  HF_data2 = HF_data [  pit_ind1 [ pit_ind2[1]:pit_ind2[2]]   , ]  
  # Exclude last row, belongs to next phase
  HF_data3 = HF_data2[ -nrow(HF_data2), ]
  HF_data3
}))          

str(HF_data4)

## This code excludes the initial seconds of the phase before we put the fish in and close the chamber.
## The beginning is expected to have no oxygen consumption, so the slope is close to zero. 
## When the fish placed in and the chamber closed, O2 consumption starts.
## We do it by a spline that has df=25, and quantify the drop in O2. 
## Essentially, we try to find the minimum point of first and second derivative. (technically that should be the point when the chamber is closed but I had to do a few tricks to avoid noise to get captured sometimes). Later the initial seconds of the slope are also excluded, as the O2 level was affected by closing the chamber.

trunkIND_df25 = sapply (1: length(HF_data4) , function(i) { print(i)
  dat = HF_data4[[i]]$O2.correct
  AA = spline(dat)
  DD =  smooth.spline (dat , df=25 )
  plot(predict(DD))
  points(dat, col="red")
  plot(predict(DD, deriv=1)$y[1:500])
  plot(predict(DD, deriv=2)$y[1:500])
  abline(h=mean(predict(DD, deriv=1)$y))
  d1 = which.min(predict(DD, deriv=1)$y[1:100])
  d2 = which.min(predict(DD, deriv=2)$y[1:100])
  c(d1,d2)
}) 

# we exclude the first 10 seconds that can be really noisy and misleading 
trunkIND = sapply ( 1:ncol(trunkIND_df25)  , function(i) { min ( trunkIND_df25 [ trunkIND_df25 [  , i ] > 10 , i ] ) })

# Plotting the start of the slope for each fish. Check all plots.
pdf (  "MMR_HF_phases_splines.pdf" ) 
par(mfrow=c(3,2))
for (i in 1 : length(trunkIND) ) { print(i)
  dat = HF_data4[[i]]$O2.correct
  DD =  smooth.spline (dat , df=25 )
  plot ( predict(DD)$y [1:900] , ylab=i, cex = 0.4, xaxp = c(0, 900, 10))
  abline ( v = trunkIND[i]+20 ) # set slope to begin 20s after the beginning (as there's some noise)
}
dev.off()

## The figures of 125, 126, 7 ,8 the slope start is not correct. In all, fish placed very late (less time in measurement).
## Adjust them by hand:

#7
dat = HF_data4[[7]]$O2.correct
DD =  smooth.spline (dat , df=25 )
plot ( predict(DD)$y [1:900] , ylab=i)
plot ( predict(DD)$y [500:550] , ylab=i)
abline(v=40)
#8
dat = HF_data4[[8]]$O2.correct
DD =  smooth.spline (dat , df=25 )
plot ( predict(DD)$y [1:900] , ylab=i)
plot ( predict(DD)$y [500:550] , ylab=i)
abline(v=50)
#125
dat = HF_data4[[125]]$O2.correct
DD =  smooth.spline (dat , df=25 )
plot ( predict(DD)$y [1:900] , ylab=i)
plot ( predict(DD)$y [400:450] , ylab=i)
abline(v=50)
#126
dat = HF_data4[[126]]$O2.correct
DD =  smooth.spline (dat , df=25 )
plot ( predict(DD)$y [1:900] , ylab=i)
plot ( predict(DD)$y [420:470] , ylab=i)
abline(v=40)

## now filling these 
trunkIND [ c(7,8,125,126) ] = c(540,550,450,460)

## each phase in a list, only the relevant portion of the MMR measurement
## Omit the 20 rows after the beginning, as in the plots above

MMR_HF_data5 = list()
for(i in 1:length(HF_data4)) {
  MMR_HF_data5[[i]] = HF_data4[[i]] [( trunkIND[i]+20):nrow(HF_data4[[i]]) , ]
}

# Save this to input into MMR_spline_method and MMR_respr1min scripts.
save(MMR_HF_data5, file = "Data/MMR_HF_data5" )

