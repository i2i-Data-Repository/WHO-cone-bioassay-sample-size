# Description: Example of Simulation Framework for the WHO cone bioassay (considering different levels of variability between nets, net pieces and replicates)
# Authors: Gemma Francesca Harvey (gemma.harvey@lstmed.ac.uk), Frank Mechan (frank.mechan@lstmed.ac.uk)


## Outstanding features: 

# smaller (2.5%) increments, 5000 simulations
# greater range of variability

############################################################################################################

#upload packages 

library(devtools)
#devtools::install_github("pcdjohnson/GLMMmisc")
library(GLMMmisc)
library(lme4)
library(foreach)
library(doParallel)
library(dplyr) 

#Input your directory

setwd("")


# How many cores does your CPU have - determining computers processing power
n_cores <- detectCores()
n_cores

#subtracting number of cores (i.e.,  8) -  so computer can operate while simulation is running
cluster <- makeCluster(n_cores - 8)
registerDoParallel(cluster)



#############################################################################################################################

# Core functions
# Just run these once

############################################################################################################################

# This function takes a set of parameters about the size of a cone test dataset and makes a dataframe with the appropriate number of rows.
# (Every time the software estimates power for a given set of conditions, this function will have to be run many times the many different datasets)  )
Cone.InputSheet.Builder <- function(nets, pieces, reps) {
  
  conedata <- expand.grid(NetType  =c("reference","candidate"),
                          NetID    = seq(1,nets,1),
                          PieceID    = seq(1,pieces,1),
                          RepID    = seq(1,reps,1))
  
  conedata$NetType  <- as.factor(conedata$NetType)
  
  conedata$Obs  <- 1:nrow(conedata)
  conedata$n    <- 5
  return(conedata)
  
}

# run this line to see an example of the dataframe builder
#TestDatasheet <- Cone.InputSheet.Builder(nets=4,pieces=2,reps=10)
#TestDatasheet

########################################################################################################################################

# This function does the actual assessment of power.

# There are two layers of this
# The inner layer takes an input sheet (made by the builder function), predicts mortality for each row,  then performs the analysis
# The p-value is stored as the output
# The outer layer performs this inner step many, many times
# The result is a big set of p-values
# It then reports what proportion of p-values were significant. (This is the power!)

Cone.PowerSimulator <- function(InputSheet, RefMort, CanMort, NetVar, PieceVar, RepVar, Nsims) {
  
  LogOdd.TreatMort <- (CanMort / (1-CanMort) ) / (RefMort / (1-RefMort) )
  
  sim.mosdata <- function(...) {
    
    Mosdata.df <- sim.glmm(design.data  = InputSheet , 
                           fixed.eff    = list(intercept = qlogis(RefMort), NetType = log(c(reference = 1, candidate = LogOdd.TreatMort))),
                           rand.V       = c(NetID = NetVar, PieceID = PieceVar, RepID=RepVar), # RepID broke it because there was only one group/value
                           distribution = "binomial")
    
    fit <- glm(cbind(response, n - response) ~ NetType, family = binomial, data = Mosdata.df)
    return(coef(summary(fit))["NetTypecandidate", "Pr(>|z|)"])
    
  }
  
  sim.pvals  <- sapply(1:Nsims, sim.mosdata)
  power.mean <- mean(sim.pvals < 0.05)
  power.CI   <- binom.test(table(factor(sim.pvals < 0.05, c(T, F))))$conf.int
  
  return(list(power.mean,power.CI))
}


# run this line to see an example of the power simulator (number of simulations has been set to a small number to make this example run quickly)
# Note that this example takes the 'TestDatasheet' made in the builder example as an Input
#output <- Cone.PowerSimulator(InputSheet=TestDatasheet, RefMort=0.39, CanMort=0.44, NetVar=0.20, RepVar=0.10, Nsims=100)
#output[1][[1]]
# The output is the power for these conditons (with a 95% confidence interval)



#######################################################################################################

# This function just combines the two previous functions together, and assesses power for all of the different combinations of parameters presented to it
# the code for making this dataframe of combinations is at the bottom of the sheet




####################################################################################################################################

### The full process

# First we have to make a dataframe of all the combinations of parameters we want (sample size, mortality differences, variability)

maxnets  <- 10 # This is the maximum number of nets (in each group) we want to consider power for

#was set with 4 before
nets.seq        <-  seq(2,maxnets,1)   # a sequence of all the numbers of nets to be assessed
pieces.seq      <- seq(2,5,1)
rep.seq         <- seq(2,5,1)    # a sequence of all the numbers of reps to be assessed

RefMort.Values  <- seq(0.50, 0.95, 0.025)    # a sequence of all the mortality values for net 1 aka reference net
CanMort.Values  <- seq(0.55, 0.99, 0.025)    # a sequence of all the mortality values for net 1 aka candidate net


NetVar.values   <- seq(0.10, 1.60, 0.10)   # a sequence of all the variance values for between-net
PieceVar.values  <-seq(0.10, 0.90, 0.10)   # a sequence of all the variance values for between-pieces
RepVar.values   <-seq(0.10, 0.90, 0.10) # a sequence of all the variance values for between-reps


MyCombos.ready <- expand.grid(Nets          =nets.seq,
                              Pieces        =pieces.seq, 
                              Replicates    =rep.seq,
                              ReferenceMort =RefMort.Values,
                              CandidateMort =CanMort.Values,
                              NetVar        =NetVar.values,
                              PieceVar      =PieceVar.values,
                              RepVar        =RepVar.values)

MyCombos.ready <- subset(MyCombos.ready, CandidateMort > ReferenceMort)

View(MyCombos.ready)



#############################################################################################################################################
start_time <- Sys.time()

my.loop <- foreach::foreach(i = 1:nrow(MyCombos.ready),.packages=c("foreach","GLMMmisc","lme4")) %dopar% {
  
  power.output <- Cone.PowerSimulator(InputSheet=Cone.InputSheet.Builder(nets     = MyCombos.ready[i,1],
                                                                         pieces   = MyCombos.ready [i,2],
                                                                         reps     = MyCombos.ready[i,3]), 
                                      RefMort       = MyCombos.ready[i,4], 
                                      CanMort       = MyCombos.ready[i,5],  
                                      NetVar        = MyCombos.ready[i,6],
                                      PieceVar      = MyCombos.ready[i,7],  
                                      RepVar        = MyCombos.ready[i,8],
                                      Nsims         = 5000) 
  return(power.output[1][[1]])
}

end_time <- Sys.time()

my.loop
MyCombos.ready$Power <- unlist(my.loop)

#determine how long simulation ran for
total_time <- end_time - start_time

###################################################################################################################################
# Saving output- into your directory
setwd("")


#write.csv(MyCombos.ready,'File_name.csv')


stopCluster(cl = cluster)
