### Functions for filtering oxygen concentration data from metabolic rate measurements
### Jenni Prokkola, Eirik Åsheim, Tutku Aykanat

# Previous method used manual curating based on 1. R2 of the total slope, and 2. visual assessment of linearity for every slope. With high-resolution data,
# R2 filtering usually only captures very bad slopes so additional filters are needed.
# Here, functions developed to identify linear slopes using a combination of
# quadratic slopes, linear regression residuals, and R2 filtering. 

# Basic functions
# Function for lm
lmf<- function(x) {
  lm(x ~ c(1:length(x)), na.action = na.exclude)
}

# Function for R2 values
rsqf <-function(x) {
  summary(lmf(x))$r.squared
}


###--- Function 1 = rsq.outliers.f ----###
# Rsq. and outlier-filter function that 1. calculates linear slope r-square value, 
# and 2. filters measurements based on r2 threshold, 3. identifies and replaces outlier-
# data points with NA (defined as residuals of linear regression), and 
# 4. excludes the whole slope (replace with 0's) if too many NAs (>10%)
# 5. saves outputs in lists instead of data frames (four oxygen probes per list element)
# Outputs several lists: r2values, datalist, NAcounts, f1.output and excl.slop.f1

# Input data is table of raw data after correction for background respiration, only M phases, one value/second (original had 0.5/s)
# bq is a vector or batches and probe quadrants (B1.Q1 etc). Batches indicate measurements from different days.

rsq.outliers.f <- function(data, bq){ 
  # Make a new list with r2 values for each quadrant in each batch
  r2values <<- list()
  lapply(bq, function(x){
    r2values[[x]] <<- data %>%
      filter(bq == x) %>%
      group_by(Phase, Chamber.No) %>%
      summarize(
      rsq = rsqf(O2.correct)) %>%
      # remove extra rows
      distinct 
  })
  
  # Make a new list of data from the dataframe
  datalist <<- list()
  lapply(bq, function(x){
  datalist[[x]] <<- data %>%
    filter(bq == x) 
  })
  
  #  Apply rsq filter to datalist
  filt.step.1 <- lapply(datalist, function(x) x %>%
    # Apply to O2.correct (bg-corrected) values of one phase in one chamber:                      
    group_by(Phase, Chamber.No) %>%
    # Replace the O2.correct data with filtered data
    mutate(O2.correct = 
             # If rsq. of the slope <0.95, replace with 0s, else keep data.
             if(rsqf(O2.correct) < 0.95) rep(0, length(O2.correct)) else (O2.correct)))
             
  # Apply an outlier filter (based on absolute values of residuals of a linear slope). Replace outliers with NA.
  filt.step.2 <- lapply(filt.step.1, function(x) x %>%
        group_by(Phase, Chamber.No) %>%
        mutate(O2.correct = 
                 # if the observation's residual from a linear slope is more than 0.05 (mg/O2) 
                 # from the fitted line, replace with NA, else keep data point.
                 ifelse(abs(resid(lmf(O2.correct))) > 0.05,
                              NA_real_, O2.correct)))
  
  # Make a list with counts of NA values from each measurement (can check later)
  NAcounts <<- lapply(filt.step.2, function(x) x %>%
    group_by(Phase, Chamber.No) %>%
      # count the observations that are NA for each slope
    summarize(count.o = sum(is.na(O2.correct))))
  
  # Make a new output with filtered data from rsq filter and add outlier-filter
  f1.output <<- lapply(filt.step.2, function(x) x %>%
    group_by(Phase, Chamber.No) %>%
    mutate(O2.correct = 
             # if the total number of NA in the slope is >10% of the length, replace whole slope with 0, else keep data
             if(sum(is.na(O2.correct)) > (0.1*length(O2.correct))) 
                          rep(0, length(O2.correct)) else
                          (O2.correct)))
 
  # Count how many slopes were removed from each fish
  excl.slop.f1 <<- lapply(f1.output, function(x) x %>%
                  group_by(Phase, Chamber.No) %>%
                    # Makes a column count-excl with 1 for each slope that was excluded
                  summarize(count.excl = sum(all(O2.correct == 0))) %>%
                  group_by(Chamber.No) %>%
                    # count the excluded slopes for each Chamber
                  summarize(count.excl = sum(count.excl))
                  )
}

###--- Function 2. = slope.lin.f ----###
# This function takes the output of function 1 and identifies linear slopes more accurately than R2 value alone. 
# It calculates linear and quadratic slopes, and  c. quadratic slope coefficient of variation across each measurement phase, and excludes slopes that do not pass the linearity threshold.

# input data = output from function 1 (f1.output)
# outputs lists linearity.data, and f2.output, and returns the number of excluded slopes for each fish

slope.lin.f <- function(data){ 
  # Define functions for slopes
  # quadratic:
  slopequad <- function(x,y,z){
    # Some slopes have been excluded, return NA for those:
    if (all(y == 0)) return(NA)
    xsq= x^2
    a = lm(y ~ x + xsq, na.action = na.exclude)$coefficients[2]
    b = lm(y ~ x + xsq, na.action = na.exclude)$coefficients[3]
    return(2*b*z+a)
  }
  
  # linear slope:
  linear <- function(x,y){
    if (all(y == 0)) return(NA)
    return(lm(y ~ x, na.action = na.exclude)$coefficients[2])
  }
  
  # CV% for quadratic slope in specified points
  slopcv <- function(x,y){
    if (all(y == 0)) return(NA)
    # To get a point every 2 min during measurement:
    points <- c(1, (c(1:6)*120)) 
    m = abs(mean(sapply(points,function(i){
      slopequad(x,y,i)
    })))
    sd = sd(sapply(points, function(i) { 
      slopequad(x,y,i)
    }))
    # return the coefficient of variation in %
    return(sd/m*100)
  }
  
  # Calculate slopes into a new list
  linearity.data <<- lapply(data, function(x) x %>%
    group_by(Phase, Chamber.No) %>%
    summarize(
      # CV of quadratic slope
      # in function(x,y) x is the number of rows (i.e. time in seconds, and y is the O2 data)
      SLOPCV = slopcv(1:length(O2.correct),O2.correct),
      # Linear slope (not used, but included for now)
      LINEAR = linear(1:length(O2.correct),O2.correct)
    ))
  
  # Filter slopes: if CV% > 15% or NA, replace slope with 0s (was manually checked to match a visual assessment of linearity)
  f2.output <<- lapply(data, function(x) x %>%
    group_by(Phase, Chamber.No) %>%
    mutate(
      # replace O2.correct
      O2.correct = 
        if(slopcv(1:length(O2.correct),O2.correct) > 15 | 
                        is.na(slopcv(1:length(O2.correct),O2.correct))) 
                    rep(0, length(O2.correct)) 
                    else (O2.correct)))
  
  # return the number of excluded slopes (for keeping track)
  return(lapply(linearity.data, function(x) x %>%
                  group_by(Chamber.No) %>%
                  summarize(count.excl = length(which(SLOPCV > 15))))
  )
}


###--- Function 3. = NAslope.f ---###
# Check if individuals have at least 20 slopes accepted by all filters,
# remove individual (replace values with 0) if not (this was deduced to improve data quality, but some outliers may remain).
# input data = f2.output
# outputs lists NAslope.counts and filt.final.dat, and returns a summary with excl = 1 for fish that were excluded

NAslope.f <- function(data) {
  # Make a list with counts of accepted and total slopes, for checking later.
  NAslope.counts <<- lapply(data, function(x) x %>% 
    group_by(Chamber.No, Phase) %>%
      # note that summarize keep only unique combinations of phase and chamber
      # not.NA is 1 when the data is not 0's
    summarize(not.NA = sum(!(all(O2.correct == 0)))) %>%
      # then regroup to separate only chambers, not phases anymore
    group_by(Chamber.No) %>%
      # Use the previously made not.NA column to count when there are too many missing slopes
    summarize(
               # How many slopes per chamber are kept
              notNAslopes = sum(not.NA),
              # Fish is excluded if it has <20 slopes, return a "yes" for those fish
              ToomanyNA = if (notNAslopes > 20 ) "no" else ("yes"),
              # Total phases for each fish before filtering
              Total = length(unique(Phase)))
  )
  
  # Filter data from individuals with too many missing slopes (replace with 0s)
  filt.final.dat <<- lapply(data, function(x) x %>% 
                            group_by(Phase, Chamber.No) %>%
                              # Making a column with yes for missing slopes
                              # the use of mutate keeps all the original data
                            mutate(is.zero = if(all(O2.correct == 0)) "yes" else ("no")) %>%
                              # regroup by chamber 
                            group_by(Chamber.No) %>%
                              # The sum of no's in is.zero is divided by the length of the slope (i.e. max(Time), same time in all slopes). If there are 20 or more missing slopes, keep data, else replace with 0s.
                            mutate(O2.correct = 
                              if(sum(is.zero == "no")/max(Time) >= 20)
                              O2.correct   
                              else rep(0, length(O2.correct))))
  
  # return a summary with excl = 1 for fish that were excluded (check that it matches NAslope.counts)
  x <- lapply(filt.final.dat, function(x) x %>% 
    group_by(Chamber.No) %>%
    summarize(excl = sum(all(O2.correct == 0))))
  return(x)
  
}


