##########################################
###--- LOAD LIBRARIES AND FUNCTION 
##########################################
rm(list=ls())
library("nimble")
library("nimbleSCR")

setwd("C:/Personal_Cloud/OneDrive/Work/Rovquant/Papers/SpatialSurvivalRates/Wolf/gitHubFiles/MapOfDeath/Functions")

source("SaveStateModel.R")
source("dCatZwhBhPsiHRSxy1Dens.R")
source("dbinomLocal_normalWolf.R")
source("calculateDensity.R")

## ----- LOAD THE DATA ----
load("C:/Personal_Cloud/OneDrive/Work/Rovquant/Papers/SpatialSurvivalRates/Wolf/gitHubFiles/MapOfDeath/Data/Data.RData")

## ----- CREATE NIMBLE OBJETS -----
#nimbleOptions(clearNimbleFunctionsAfterCompiling = TRUE)
model <- nimbleModel( code = modelCode
                      , constants = nimConstants
                      , data = nimData
                      , inits = nimInits
                      , check = FALSE       
                      , calculate = FALSE)

cmodel <- compileNimble(model)
calc <- cmodel$calculate()  
print(calc)

conf <- configureMCMC(model,control = list(reflective = TRUE, adaptScaleOnly = TRUE),
                      useConjugacy = FALSE)

conf$addMonitors(nimParams,"sxy","z","sex")   ##

##ADD THE CATEGORICAL SAMPLER TO z (but not to z[,1])
###
samplerConfList <- unlist(lapply(conf$getSamplers(),function(x) x$target))
zNodes <- samplerConfList[grep("z",samplerConfList)]
#find z nodes from the first year
zNodes <- zNodes[-grep( ", 1]",zNodes)]
conf
#remove samplers
conf$removeSamplers(zNodes)

for(i in 1:length(zNodes)){
  conf$addSampler(target = zNodes[i], type = 'sampler_categorical_general',
                  control=list("numCategories"= 5))
}

## CONFIGURE THE REVERSIBLE JUMP
for(n in 1:5){
  configureRJ(conf,
              targetNodes = paste('betaWCountry[',n+2,']',sep=""),
              indicatorNodes = paste('zW[',n,']',sep=""),
              control = list(mean = 0, scale = .2))
}


for(n in 1:length(targetNodes)){
  configureRJ(conf,
              targetNodes = targetNodes[n],
              indicatorNodes = indicatorNodes[n],
              control = list(mean = 0, scale = .2))
}
#CHECK THE SAMPLERS
conf$printSamplers(c("betaHRoads", "zH[1]"))
conf$printSamplers(c("betaHMoose", "zH[3]"))



Rmcmc <- buildMCMC(conf)
cMCMC <- compileNimble(Rmcmc, project = model, resetFunctions = TRUE)
# RUN THE MCMC
MCMCRuntime <- system.time(myNimbleOutput <- runMCMC( mcmc = cMCMC,
                                                      nburnin = 1,
                                                      niter = 20,#RUN FOR LONGER 
                                                      nchains = 1,
                                                      samplesAsCodaMCMC = TRUE))










