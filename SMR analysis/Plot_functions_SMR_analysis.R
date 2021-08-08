#### Plot functions for multi-batch respirometry analysis.
### Jenni Prokkola, Eirik Åsheim

# Plot 1 - rawdat.multiplot. Makes a multipanel plot of slopes (O2 concentration vs. Time) that have data. Adds a linear regression line in green.
# Input: list that can have data for multiple (4) individuals in each element, in long format

rawdat.multiplot <- function(data) {
	# using seq_along which allows using list element names in plot title
  lapply(seq_along(data), function(x){
    # exclude fish with missing data:
    dat <- data[[x]] %>%
      group_by(Chamber.No) %>%
      filter(!(all(O2.correct %in% c(NA, 0))))
    dat <- droplevels(dat)
    chamber = levels(dat$Chamber.No)
    for (c in chamber){
      plotdat = subset(dat, Chamber.No == c)
      #Exclude missing phases from plots
      plotdat = plotdat %>%
        group_by(Phase) %>%
        filter(!(all(O2.correct %in% c(NA, 0))))
        
      p <- ggplot(data = plotdat, 
                  aes(x = Time, y = O2.correct ))+
        facet_wrap(~Phase)+ # Plots different phases of one fish on the same page
        geom_point(cex = 0.4)+
        geom_line(stat="smooth",method = lm, colour = "green", na.rm = T)+
        ylab("O2 mg")+
        ggtitle(paste(names(data)[[x]], c))+ # combine element name and chamber number in title
        theme_bw()+
        theme(axis.text.x = element_text(angle = 90, size = 8),
              axis.title.x = element_blank())
      
      print(p)
    }
  })
}

## Plot 2 - MR.dat.plot. Final mass-specific data for each chamber over time
# input is list of MO2 in all phases
# For multiple batches:
MR.dat.plot <- function(data) {
  lapply(seq_along(data), function(x){
    dat <- data[[x]] 
    dat <- droplevels(dat)
    # Helpful to get date and time in the right format (package lubridate)
    dat$Date.Time <- ymd_hms(dat$Date.Time)
    p <-ggplot(data = dat, 
               aes(x = Date.Time, y = MR.mass )) + #plotting mass-specific values 
      facet_wrap(~Chamber.No)+ # four fish on the same page (4 fish in one list element)
      geom_point(cex = 0.5)+
      ylab("O2 mg/kg/h")+
      ggtitle(names(data)[[x]])+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 90, size = 8),
            axis.title.x = element_blank())
    
  })
}

## Plot 3 - plot.mcclust.f. Visualizing mcclust function for selecting SMR
# Input is lists of MO2 from all phases for each individual (data), the value of SMR result (result), and number of individuals to plot (N)
plot.mcclust.f <- function(data, result, N){
  dat<-rbindlist(data) # make the list into a data frame
  res <- rbindlist(result)
  # number of ind plotted:
  for(i in 1:N) {
    mr <- dat %>% filter(Ind == unique(dat$Ind)[i]) %>% pull(MR.mass) # mass-specific for visualizing
    mlnd <- res %>% filter(PIT == unique(dat$Ind)[i]) %>% pull(SMR.mass.mlnd)
    plot(densityMclust(mr,G=1:4), what = "density", type = "hdr", data=mr,breaks=150)
    abline(v = mlnd, col = "blue")
    title(main = res %>% filter(PIT == unique(dat$Ind)[i]) %>% pull(PIT))
    Sys.sleep(time=0.05)
  }
}

## Plot 4 - SMR.weight.plot. SMR vs. fish weight with different metrics, linear and log-scale.
# input is a list of SMR results (one value for each fish, with a column for mass in the same table). varnames is a vector of variables to plot (if using several SMR metrics, for example)
# For multiple batches:
SMR.weight.plot <- function(data, varnames) {
  
  # Make a data frame of result list
  dat<- rbindlist(data)
  plotdat <- dat %>%
    pivot_longer(cols = -c(Chamber.No, Mass, Temp, PIT, Date), names_to = "SMR_var",
                 values_to = "SMR")
  
  plotdat <- filter(plotdat, SMR_var %in% varnames)
  # plot on log scale
  p <- ggplot(data = plotdat, aes(x = Mass, y = SMR ))+
    facet_wrap(~SMR_var)+
    geom_point(cex = 0.5)+
    ylab("log10(SMR abs)")+
    xlab("log10(Mass)")+
    coord_trans(x = 'log10',y = 'log10') + #plot on log scale
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, size = 8))
   # plot on linear scale
  n <- ggplot(data = plotdat, aes(x = Mass, y = SMR ))+
    facet_wrap(~SMR_var)+
    geom_point(cex = 0.5)+
    ylab("SMR abs")+
    xlab("Mass")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, size = 8))
  print(n)
  print(p)
  
}



