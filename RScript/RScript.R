## ---- LOAD LIBRARIES AND FUNCTION ----
rm(list=ls())
library("nimble")
library("nimbleSCR")
library("coda")
library("tidyverse")
library("sf")
library("viridis")
library("raster")
## ---- LOAD  FUNCTIONs ----
setwd("C:/Personal_Cloud/OneDrive/Work/Rovquant/Papers/SpatialSurvivalRates/Wolf/gitHubFiles/MapOfDeath/Functions")
source("SaveStateModel.R")
source("dCatZwhBhPsiHRSxy1Dens.R")
source("dbinomLocal_normalWolf.R")
source("calculateDensity.R")
source("ProcessCodaOutput.R")

# sourceCpp("GetDensity.cpp")

## ---- LOAD THE DATA ----
load("C:/Personal_Cloud/OneDrive/Work/Rovquant/Papers/SpatialSurvivalRates/Wolf/gitHubFiles/MapOfDeath/Data/Data.RData")
load("C:/Personal_Cloud/OneDrive/Work/Rovquant/Papers/SpatialSurvivalRates/Wolf/gitHubFiles/MapOfDeath/Data/NecessaryFiles.RData")
nYears <- dim(nimData$z)[2]
years <- 2015:2021

### ==== 1. RUN NIMBLE MODEL ====
### ====  1.1 CREATE NIMBLE OBJETS ====
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
                                                      niter = 50,#RUN FOR LONGER 
                                                      nchains = 2,
                                                      samplesAsCodaMCMC = TRUE))


### ==== 2. PROCESS THE POSTERIORS AND PLOT ====
#Use a few posteriors to check 
#setwd("C:/Personal_Cloud/OneDrive/Work/Rovquant/Papers/SpatialSurvivalRates/Wolf/gitHubFiles")
## save(myNimbleOutput,file="MCMCSamples.RData")
#load("MCMCSamples.RData")
#Process the results, get them together 
myResults <- ProcessCodaOutput(myNimbleOutput,params.omit = c("sxy","z"))
### ====  2.1 TABLE 1 ====
#INDICATOR VARIABLES
colnames(myResults$sims.list$zH) <- c("RoadsH","HistPresH", "MooseH", "humanH", "snowH","DensH")
colMeans(myResults$sims.list$zH)
colnames(myResults$sims.list$zW) <- c("RoadsW","HistPresW", "MooseW", "humanW", "snowW","DensW")
colMeans(myResults$sims.list$zW)

myResults$sims.list$zWH <- cbind(myResults$sims.list$zW,myResults$sims.list$zH)
myResults$sims.list$betaWH <- cbind(myResults$sims.list$betaWCountry[,3:7],
                                    myResults$sims.list$betaWDens,
                                    myResults$sims.list$betaHRoads ,
                                    myResults$sims.list$betaHHistPres,
                                    myResults$sims.list$betaHMoose ,
                                    myResults$sims.list$betaHHuman,
                                    myResults$sims.list$betaHSnow ,
                                    myResults$sims.list$betaHDens)

colnames(myResults$sims.list$betaWH ) <- c("RoadsW","HistPresW", "MooseW", "humanW", 
                                           "snowW","DensW","RoadsH","HistPresH", "MooseH", "humanH",
                                           "snowH","DensH")

##GET AVERAGE ESTIMATES WHEN RJ INDICATOR IS 1
betaWH1 <- myResults$sims.list$betaWH
# SET BETA WHEN rj=0 TO NA 
betaWH1[myResults$sims.list$zWH%in%0 ] <- NA

# GET POSTERIOR INCLUSION PROBABILITY
PIP <- colMeans(myResults$sims.list$zWH)
PIP <- sort(PIP, decreasing = T)
PIPTable <- matrix(NA, nrow=length(PIP), ncol=1)
PIPTable[,1] <- round(PIP,digits = 2)
row.names(PIPTable) <- names(PIP)
colnames(PIPTable) <- c("table1")


## Separate COEFF and CI in the table ##
table1 <- matrix(NA,nrow=ncol(betaWH1), ncol=4)
colnames(table1) <- c("Covariate","Coeff","95% CrI", "PIP")
row.names(table1) <- table1[,1] <- colnames(betaWH1)

table1[,"Coeff"] <- round(colMeans(betaWH1, na.rm = T),digits = 2) 
table1[,"95% CrI"] <- paste(round(apply(betaWH1,2,function(x) quantile(x,probs=c(0.025,0.975),na.rm=T)),digits = 2)[1,],
               ";",
               round(apply(betaWH1,2,function(x) quantile(x,probs=c(0.025,0.975),na.rm=T)),digits = 2)[2,],
               sep=""
)
table1[row.names(PIPTable),"PIP"] <- round(PIPTable,digits = 2)
table1 <- table1[row.names(table1),] 
row.names(table1) <- rep(c("W","H"),each=6)
table1 <- table1[order(table1[,"PIP"], row.names(table1),decreasing = T),] 
table1 <- table1[order(row.names(table1)),] 
print(table1)

### ====  2.2 TABLE.S3  ====
# FIND THE BEST MODELS
listBest <- apply(myResults$sims.list$zWH,1, function(x) which(x>0))
listBest1Cause <- list()
for(i in 1:nrow(myResults$sims.list$zH)){
  if(length(names(listBest[[i]]))==0){
  }
  cause <- substr(names(listBest[[i]]),nchar(names(listBest[[i]])),
                  nchar(names(listBest[[i]])) )
  which(cause %in% "W")
  which(cause %in% "H")
  listBest1Cause[[i]] <- c( paste(names(listBest[[i]])[ which(cause %in% "W")],collapse=" + "),
                            paste(names(listBest[[i]])[ which(cause %in% "H")],collapse=" + "))
  
}
# bind cause specific variables
unlistBestMod <- do.call(rbind, listBest1Cause)
UniqueMod <- apply(unlistBestMod,1,function(x) paste(x[1],x[2]))
tabBestModCause <- table(UniqueMod)
bestModelsCause <- sort(tabBestModCause,decreasing = T)

###make table s3
sum((bestModelsCause/sum(tabBestModCause))[1:20])
tab <- (bestModelsCause/sum(tabBestModCause))[1:20]
m <- matrix(NA,nrow=20,ncol=1)
row.names(m) <- names(tab)
colnames(m) <- c("Model weight")
m[,1] <- round(tab,digits=3)

tableS3 <- matrix(NA,nrow = nrow(m),ncol=2)
for(i in 1:nrow(m)){
  tableS3[i,] <-  unlistBestMod[which(UniqueMod %in%  row.names(m)[i])[1],]
}
#bind model weight and variables names
tableS3 <- cbind(tableS3, m)
row.names(tableS3) <- NULL
tableS3



### ====  2.3 TABLE.S4  ====
tableS4 <- matrix(NA, nrow=4,ncol=nYears-1)
colnames(tableS4) <- unlist(years[1:(nYears-1)])
rownames(tableS4) <- c("BetaHWolfZone","BetaHReindeer",
                       "BetaWWolfZone","BetaWReindeer") 

tableS4[1,] <-   paste(round(colMeans(myResults$sims.list$betaHCountry[,1,],na.rm = T),digits = 2), 
                       " (", 
                       round(apply(myResults$sims.list$betaHCountry[,1,],2,function(x) quantile(x,probs=c(0.025,0.975),na.rm=T)),digits = 2)[1,],
                       ";",
                       round(apply(myResults$sims.list$betaHCountry[,1,],2,function(x) quantile(x,probs=c(0.025,0.975),na.rm=T)),digits = 2)[2,],
                       ")", sep=""
)

tableS4[2,3] <-   paste(round(mean(myResults$sims.list$betaHReindeer,na.rm = T),digits = 2), 
                        " (", 
                        round(quantile(myResults$sims.list$betaHReindeer,probs=c(0.025,0.975),na.rm=T),digits = 2)[1],
                        ";",
                        round(quantile(myResults$sims.list$betaHReindeer,probs=c(0.025,0.975),na.rm=T),digits = 2)[2],
                        ")", sep=""
)


tableS4[3,3] <-    paste(round(mean(myResults$sims.list$betaWCountry[,1],na.rm = T),digits = 2), 
                         " (", 
                         round(quantile(myResults$sims.list$betaWCountry[,1],probs=c(0.025,0.975),na.rm=T),digits = 2)[1],
                         ";",
                         round(quantile(myResults$sims.list$betaWCountry[,1],probs=c(0.025,0.975),na.rm=T),digits = 2)[2],
                         ")", sep=""
)



tableS4[4,3] <-    paste(round(mean(myResults$sims.list$betaWCountry[,2],na.rm = T),digits = 2), 
                         " (", 
                         round(quantile(myResults$sims.list$betaWCountry[,2],probs=c(0.025,0.975),na.rm=T),digits = 2)[1],
                         ";",
                         round(quantile(myResults$sims.list$betaWCountry[,2],probs=c(0.025,0.975),na.rm=T),digits = 2)[2],
                         ")", sep=""
)
tableS4

### ====  2.4 TABLE.S5  ====
tableS5 <- matrix(0, nrow=4, ncol=(nYears-1)*4)
colnames(tableS5) <- rep(years[1:(nYears-1)],each=4)
tableS5[1,] <- rep(c("M","M","F","F"),nYears-1)
tableS5[2,] <- rep(c("1","2"),nYears-1)

sex <- c("F","M")
for(t in 1:(nYears-1)){
  for(s in 1:2){
    for(st in 1:2){
      wh <-  which(tableS5[1,]%in% sex[s] & colnames(tableS5) %in% as.character(years[t]) & tableS5[2,] %in% as.character(st))
      tableS5[3,wh] <- 
        paste(round(median(myResults$sims.list$mhH[,st,s,t],na.rm = T),digits = 1), 
              " (", 
              round(quantile(myResults$sims.list$mhH[,st,s,t],probs=c(0.025,0.975),na.rm=T),digits = 1)[1],
              ";",
              round(quantile(myResults$sims.list$mhH[,st,s,t],probs=c(0.025,0.975),na.rm=T),digits = 1)[2],
              ")", sep=""
        )
      
      tableS5[4,wh] <- 
        paste(round(median(myResults$sims.list$mhW[,st,s,t],na.rm = T),digits = 1), 
              " (", 
              round(quantile(myResults$sims.list$mhW[,st,s,t],probs=c(0.025,0.975),na.rm=T),digits = 1)[1],
              ";",
              round(quantile(myResults$sims.list$mhW[,st,s,t],probs=c(0.025,0.975),na.rm=T),digits = 1)[2],
              ")", sep=""
        ) 
    }
  }
}

row.names(tableS5) <- c("sex","State","h0","w0")

tableS5

### ====  2.5 TABLE.S6  ====
tableS6 <- matrix(0, nrow=2, ncol=(nYears-1))
colnames(tableS6) <- rep(years[1:(nYears-1)],1)
row.names(tableS6) <- c("rho","psi")
for(t in 1:(nYears-1)){
  
  n.recruit <- apply(myResults$sims.list$z[,,c(t,t+1)], 1, function(x) sum( x[,1] %in% c(1) & x[,2] %in% c(2) ))
  # NUMBER of individuals with state 1 at t-1 and 2 at t==> number of recruits
  alivetminus1 <- apply(myResults$sims.list$z[,,t], 1, function(x)sum(x %in% c(2,3)))
  temp <- n.recruit/alivetminus1
  
  
  tableS6[1,t] <-     paste(round(median(temp,na.rm = T),digits = 2), 
                      " (", 
                      round(quantile(temp,probs=c(0.025,0.975),na.rm=T),digits = 2)[1],
                      ";",
                      round(quantile(temp,probs=c(0.025,0.975),na.rm=T),digits = 2)[2],
                      ")", sep=""
  ) 
  
  
  tableS6[2,t] <-     paste(round(median(myResults$sims.list$psi[,t],na.rm = T),digits = 1), 
                      " (", 
                      round(quantile(myResults$sims.list$psi[,t],probs=c(0.025,0.975),na.rm=T),digits = 2)[1],
                      ";",
                      round(quantile(myResults$sims.list$psi[,t],probs=c(0.025,0.975),na.rm=T),digits = 2)[2],
                      ")", sep=""
  ) 
  
}

tableS6

### ====  2.6 TABLE.S7  ====
# habitat resolution needed to rescale units 
habitatResolution <- 20000
tableS7 <- matrix(0, nrow=3, ncol=(nYears)*4)
colnames(tableS7) <- rep(years[1:(nYears)],each=4)
tableS7[1,] <- rep(c("M","M","F","F"),nYears)
tableS7[2,] <- rep(c("1","2"),nYears)

sex <- c("F","M")
for(t in 1:(nYears)){
  for(s in 1:2){
    for(st in 1:2){
      wh <-  which(tableS7[1,]%in% sex[s] & colnames(tableS7) %in% as.character(years[t]) & tableS7[2,] %in% as.character(st))
      tableS7[3,wh] <- 
        paste(round(median(myResults$sims.list$sigma[,st,s,t]*habitatResolution/1000,na.rm = T),digits = 1), 
              " (", 
              round(quantile(myResults$sims.list$sigma[,st,s,t]*habitatResolution/1000,probs=c(0.025,0.975),na.rm=T),digits = 1)[1],
              ";",
              round(quantile(myResults$sims.list$sigma[,st,s,t]*habitatResolution/1000,probs=c(0.025,0.975),na.rm=T),digits = 1)[2],
              ")", sep=""
        )
    }
  }
}

row.names(tableS7) <- c("sex","State","sigma")
tableS7


### ====  2.7 TABLE.S8  ====
tableS8 <- matrix(NA, nrow=3, ncol=2)
tableS8[1,] <- rev(sex)
row.names(tableS8) <- c("Sex",1,2)
for(st in 1:2){
  for(s in 1:2){
    tableS8[st+1,s] <- paste(round(median(myResults$sims.list$dmean[,st,s]*habitatResolution/1000,na.rm = T),digits = 2), 
                       " (", 
                       round(quantile(myResults$sims.list$dmean[,st,s]*habitatResolution/1000,probs=c(0.025,0.975),na.rm=T),digits = 2)[1],
                       ";",
                       round(quantile(myResults$sims.list$dmean[,st,s]*habitatResolution/1000,probs=c(0.025,0.975),na.rm=T),digits = 2)[2],
                       ")", sep=""
    )
  }
}

tableS8

### ====  2.8 TABLE.S9  ====
tableS9 <- matrix(0, nrow=8, ncol=(nYears)*2)
colnames(tableS9) <- rep(years[1:(nYears)],each=2)
tableS9[1,] <- rep(c("M","F"),nYears)

row.names(tableS9) <- c("Sex","BEffort","BSnow","BResponse", "BOthRoad","BOthObs","BOthSnow","BOthResponse")

sex <- c("F","M")
for(t in 1:(nYears)){
  for(s in 1:2){
    
    
    for(j in 1:2){
      wh <-  which(tableS9[1,]%in% sex[s] & colnames(tableS9) %in% as.character(years[t]))# & m[2,] %in% as.character(st))
      tableS9[j+1,wh] <- 
        paste(round(median(myResults$sims.list$trapBetas[,j,s,t],na.rm = T),digits = 2), 
              " (", 
              round(quantile(myResults$sims.list$trapBetas[,j,s,t],probs=c(0.025,0.975),na.rm=T),digits = 2)[1],
              ";",
              round(quantile(myResults$sims.list$trapBetas[,j,s,t],probs=c(0.025,0.975),na.rm=T),digits = 2)[2],
              ")", sep=""
        )
    }
    
    tableS9[4,wh] <- 
      paste(round(median(myResults$sims.list$betaResponse[,s,t],na.rm = T),digits = 2), 
            " (", 
            round(quantile(myResults$sims.list$betaResponse[,s,t],probs=c(0.025,0.975),na.rm=T),digits = 2)[1],
            ";",
            round(quantile(myResults$sims.list$betaResponse[,s,t],probs=c(0.025,0.975),na.rm=T),digits = 2)[2],
            ")", sep=""
      )
    
    
    
    for(j in 1:3){
      wh <-  which(tableS9[1,]%in% sex[s] & colnames(tableS9) %in% as.character(years[t]))# & m[2,] %in% as.character(st))
      tableS9[j+4,wh] <- 
        paste(round(median(myResults$sims.list$trapBetasOth[,j,s,t],na.rm = T),digits = 2), 
              " (", 
              round(quantile(myResults$sims.list$trapBetasOth[,j,s,t],probs=c(0.025,0.975),na.rm=T),digits = 2)[1],
              ";",
              round(quantile(myResults$sims.list$trapBetasOth[,j,s,t],probs=c(0.025,0.975),na.rm=T),digits = 2)[2],
              ")", sep=""
        )
    }
    
    
    tableS9[8,wh] <- 
      paste(round(median(myResults$sims.list$betaResponseOth[,s,t],na.rm = T),digits = 2), 
            " (", 
            round(quantile(myResults$sims.list$betaResponseOth[,s,t],probs=c(0.025,0.975),na.rm=T),digits = 2)[1],
            ";",
            round(quantile(myResults$sims.list$betaResponseOth[,s,t],probs=c(0.025,0.975),na.rm=T),digits = 2)[2],
            ")", sep=""
      )
  }
}

tableS9
### ====  2.9 Result section: different number and % of mortality events ====
## number of individuals that were found dead. 
## state set to 4 and 5 when found dead. 
ZDEAD <- ZDEAD1 <- nimData$z
ZDEAD[] <- ZDEAD1[] <- 0
for(i in 1:nrow(ZDEAD)){
  if(sum(nimData$z[i,] %in% c(4,5))>0){
    firstDead <- which(nimData$z[i,] %in% c(4,5))[1]
    ZDEAD[i,firstDead] <- 1
  }
}
### sum of ids found dead ###
colSums(ZDEAD)

#GET POSTERIOR OF NUMBER OF INDIVIDUALS ESTIMATED TO BE DEAD #
NDEADPerYEAR <- list() 
for(ite in 1:dim(myResults$sims.list$z)[1]){
  ZDEAD1[] <- 0
  for(i in 1:nrow(ZDEAD)){
    
    tmpz <- myResults$sims.list$z[ite,i,] %in% c(4,5)
    if(sum(tmpz)>0){
      firstDead <- which(tmpz)[1]
      ZDEAD1[i,firstDead] <- 1
    }
  }
  NDEADPerYEAR[[ite]] <- colSums(ZDEAD1)
}
#GET AVERAGE NUMBER OF INDIVIDUALS DEAD NOT DETECTED
NDEADPerYEARMatrix <- do.call(rbind, NDEADPerYEAR)
averageNumberNotDetected <- list()
for(t in 2:nYears ){
  averageNumberNotDetected[[t]]  <- colSums(ZDEAD)[t]/NDEADPerYEARMatrix[,t]
}
## get yearly mean
lapply(averageNumberNotDetected, mean)
## get yearly 95CrI
lapply(averageNumberNotDetected, function(x){quantile(x, probs=c(2.5,97.5)/100)})

##estimated number of wolves that diead yearly
colMeans(NDEADPerYEARMatrix)
apply(NDEADPerYEARMatrix,2,function(x) quantile(x,probs=c(0.025,0.975)))

###annual proportion of wolves that died yearly 
colMeans(NDEADPerYEARMatrix)[2:(nYears)]/colMeans(myResults$sims.list$N)[1:(nYears-1)]
#CRI
apply(NDEADPerYEARMatrix,2,function(x) quantile(x,probs=c(0.025,0.975)))[,2:(nYears)]/
  apply(myResults$sims.list$N,2,function(x) quantile(x,probs=c(0.025,0.975)))[,1:(nYears-1)]

#overall proportion of wolves that died during the entire period
mean(colMeans(NDEADPerYEARMatrix)[2:(nYears)]/colMeans(myResults$sims.list$N)[1:(nYears-1)])



#Average number of individuals that remained undetected
1-mean(do.call(c, averageNumberNotDetected))
1-quantile(do.call(c, averageNumberNotDetected),probs=c(0.025,0.975))



### ==== 3. MORTALITY MAPS ====
### ====  3.1 GET SPATIAL OBJECTS READY FOR PREDICTIONS ====
habCov.arr <- nimConstants$habCov
# OBTAIN MAPS AT 5KMS 
dev.off()
habbR <- raster::disaggregate(habitat.r, fact=2)
plot(habbR)

# CREATE COUNTRY OBJECTS FOR EXTRACTION 
#COUNTRIES <- aggregate(x = COMMUNES, by = "ISO")
COUNTRIES <- COMMUNESsf %>% group_by(ISO) %>% summarize()
COUNTRIES$NAME_1 <- COUNTRIES$ISO
plot(habitat.r)
plot(COUNTRIES$geometry, add=T)

## IDENTIFY COUNTIES
box <- st_as_sfc(st_bbox(habitat.poly))

## GET OBJECTS READY FOR THE DENSITY EXTRACTION ##
habIDCells.mx <- habbR
habIDCells.mx[] <- 1:ncell(habbR)
habIDCells.mx <- as.matrix(habIDCells.mx)

DensiCountriesCPP <- list()
regionID <- habbR[]
regionIDunique <- unique(regionID)
regionIDunique <- regionIDunique[!is.na(regionIDunique)]
regionIDmat <-do.call(rbind,lapply(regionIDunique,function(x)habbR[]== x  ))
regionIDmat[is.na(regionIDmat)] <- 0
row.names(regionIDmat) <- unique(regionIDunique)


### rescale sxy (ACs) coordinates to original scale
dimnames(myResults$sims.list$sxy)[[3]] <- c("x", "y")
sxyScaledOri  <- myResults$sims.list$sxy 
myResults$sims.list$sxy <- scaleCoordsToHabitatGrid(coordsData = myResults$sims.list$sxy,
                                                        coordsHabitatGridCenter = habitat.xy,
                                                        scaleToGrid = FALSE)$coordsDataScaled



xyhabNew <- coordinates(habbR)
colnames(xyhabNew) <-  c("x", "y")
rescaledSXYNew <- scaleCoordsToHabitatGrid(coordsData = myResults$sims.list$sxy,#rescaledSXY$coordsDataScaled,
                                           coordsHabitatGridCenter =  xyhabNew,
                                           scaleToGrid = T
)

### ====    3.3.1 GET DENSITY  ==== ####
# DensityCountriesRegions <- list()
# for(t in 1:nYears){
#   DensityCountriesRegions[[t]] <- GetDensity_PD(
#     sx = rescaledSXYNew$coordsDataScaled[,,1,t],
#     sy =  rescaledSXYNew$coordsDataScaled[,,2,t],
#     z = myResults$sims.list$z[,,t],
#     IDmx = habIDCells.mx,
#     aliveStates = c(2,3),
#     regionID = regionIDmat,
#     returnPosteriorCells = T)
# }
# 
### ====    3.3.2 PLOT DENSITY MAPS  ==== ####
# SpatialRasterDensity <- list()
# for(t in 1:(nYears-1)){
#   SpatialRasterDensity[[t]] <- habbR
#   SpatialRasterDensity[[t]][] <- DensityCountriesRegions[[t]]$MeanCell
#   plot(SpatialRasterDensity[[t]])
#   plot(habitat.poly$geometry,add=T)
# }


# chose iterations from which to calculate realized density
itera <- dim(myResults$sims.list$sxy)[1]

# compile calculate Density for faster computation
CcalculateDensity <- compileNimble(calculateDensity)

#initiate array
dens <- dens1 <- array(0,c(dim(myResults$sims.list$sxy)[1],
                           nimConstants$numHabWindows,
                           nimConstants$n.years))
ite=1
# USE THE COMPILED NIMBLE FUNCTION
for(t in 1:nYears){
  for(ite in 1:dim(dens)[1]){
    dens[ite,,t] <-  CcalculateDensity(s = sxyScaledOri[ite, ,1:2, t],
                                       habitatGrid = nimData$habitatGrid,
                                       indicator =  myResults$sims.list$z[ite, , t] %in% c(2,3),
                                       numWindows = nimConstants$numHabWindows,
                                       nIndividuals = nimConstants$n.individuals)
    dens1[ite,, t] <- log(dens[ite,,t] + 1) - nimConstants$dens.offset
    
  }
}

### ====    3.3.3 PREDICT MAPS  ====
### ====      3.3.3.1 CALCULATE ====
# t = 3
# ss=2
# s=1
# Initiate object for storage
H <- W <- PHI <- list()
HCIL <- WCIL <- PHICIL <- list()
HCIH <- WCIH <- PHICIH <- list()

mhH1 <- mhW1 <- phi <- h <- w <- phi <- array(0,c(itera,2,2,
                                                  nimConstants$numHabWindows,
                                                  nYears))
# get predictions  for each sex (s), state (ss), year (t), and iteration (itera).
for(s in 1:2){
  for(ss in 1:2){
    for(t in 1:(nYears-1)){
      for(ite in 1:itera){
        mhH1[ite,s,ss,,t] <- exp(myResults$sims.list$mhH[ite,s,ss,t] + 
                                   myResults$sims.list$betaHDens[ite] * dens1[ite,,t] +
                                   myResults$sims.list$betaHCountry[ite,1,t] * habCov.arr[,"WolfZone",t] +
                                   myResults$sims.list$betaHReindeer[ite] * habCov.arr[,"Reinder",t] +
                                   myResults$sims.list$betaHRoads[ite] * habCov.arr[,"Roads",t] +
                                   myResults$sims.list$betaHHistPres[ite] * habCov.arr[,"HistPres",t] +
                                   myResults$sims.list$betaHMoose[ite] * habCov.arr[,"Moose",t]+ 
                                   myResults$sims.list$betaHHuman[ite] * habCov.arr[,"HumanDens",t]+
                                   myResults$sims.list$betaHSnow[ite] * habCov.arr[,"Snow",t])
        
        mhW1[ite,s,ss,,t] <- exp(myResults$sims.list$mhW[ite,s,ss,t] + 
                                   myResults$sims.list$betaWDens[ite] * dens1[ite,,t] +
                                   myResults$sims.list$betaWCountry[ite,1] * habCov.arr[,"WolfZone",t] +
                                   myResults$sims.list$betaWCountry[ite,2] * habCov.arr[,"Reinder",t] +
                                   myResults$sims.list$betaWCountry[ite,3] * habCov.arr[,"Roads",t] +
                                   myResults$sims.list$betaWCountry[ite,4] * habCov.arr[,"HistPres",t] +
                                   myResults$sims.list$betaWCountry[ite,5] * habCov.arr[,"Moose",t]+ 
                                   myResults$sims.list$betaWCountry[ite,6] * habCov.arr[,"HumanDens",t]+
                                   myResults$sims.list$betaWCountry[ite,7] * habCov.arr[,"Snow",t])
        
        #print(ite)
        # transform hazard rates to probabilities 
        #legal mortality
        h[ite,s,ss,,t] <- (1-exp(-(mhH1[ite,s,ss,,t] + mhW1[ite,s,ss,,t]))) * (mhH1[ite,s,ss,,t]/(mhH1[ite,s,ss,,t]  + mhW1[ite,s,ss,,t]))
        # other causes
        w[ite,s,ss,,t] <- (1-exp(-(mhH1[ite,s,ss,,t] + mhW1[ite,s,ss,,t]))) * (mhW1[ite,s,ss,,t]/(mhH1[ite,s,ss,,t]  + mhW1[ite,s,ss,,t]))
        #overal mortality 
        phi[ite,s,ss, ,t] <- h[ite,s,ss,,t] + w[ite,s,ss,,t] 
      }
    }}}

### ====      3.3.3.2 PLOT ====
#plot the predicted maps 
r <- r1 <-  habitat.r
r[r1[]==1] <-  colMeans(dens[,,t])#habCov.arr[,"Reinder",t]#
r[r1[]==0] <- NA
st <- c("other","scent-marking")
sex <- c("F","M")


cuts <- seq(0,1,by=0.01)
col = terrain.colors(length(cuts))
par(mfrow=c(3,nYears-1),mar=c(1,1,1,3))
state <- c("other","scent-marking") 
for(s in 1:2){
  for(ss in 1:2){
    for(t in 1:(nYears-1)){
      #Legal Mortality
      H[[t]] <- W[[t]] <- PHI[[t]] <-  r1 <- habitat.r
      H[[t]][r1[]==1] <- colMeans(h[,s,ss, ,t])
      H[[t]][r1[]==0] <- NA
      plot(H[[t]], breaks=cuts, col = col, main="",legend=F,axes=F ,box=F) #p
      ####CI
      HCIH[[t]] <- WCIH[[t]] <- PHICIH[[t]] <-  r1 <- habitat.r
      HCIL[[t]] <- WCIL[[t]] <- PHICIL[[t]] <-  r1 <- habitat.r
      
      H[[t]][r1[]==1] <- colMeans(h[,s,ss, ,t])
      H[[t]][r1[]==0] <- NA
      #CI
      HCIH[[t]][r1[]==1] <- apply(h[,s,ss, ,t],2, function(x) quantile(x,prob=97.5/100,na.rm=T))
      HCIH[[t]][r1[]==0] <- NA
      
      HCIL[[t]][r1[]==1] <- apply(h[,s,ss, ,t],2, function(x) quantile(x,prob=2.5/100,na.rm=T))
      HCIL[[t]][r1[]==0] <- NA
      
      if(t==3 ){
        mtext(paste(sex[s],state[ss], sep= "-"),3)
      } 
      if(t==1 ){
        mtext("h",2)
      } 
    }
    # Other mortality
    for(t in 1:(nYears-1)){
      W[[t]][r1[]==1] <- colMeans(w[,s,ss, ,t])
      W[[t]][r1[]==0] <- NA
      #CI
      WCIH[[t]][r1[]==1] <- apply(w[,s,ss, ,t],2, function(x) quantile(x,prob=97.5/100,na.rm=T))
      WCIH[[t]][r1[]==0] <- NA
      #CI
      WCIL[[t]][r1[]==1] <- apply(w[,s,ss, ,t],2, function(x) quantile(x,prob=2.5/100,na.rm=T))
      WCIL[[t]][r1[]==0] <- NA
      
      plot(W[[t]], breaks=cuts, col = col, main="",legend=F,axes=F ,box=F) #p
      if(t==1 ){
        mtext("w",2)
      } 
    }
    
    # overallmortality# 
    for(t in 1:(nYears-1)){
      PHI[[t]][r1[]==1] <- exp(-(colMeans(mhH1[,s,ss, ,t])+colMeans(mhW1[,s,ss, ,t])))#colMeans(1-phi)
      PHI[[t]][r1[]==0] <- NA
      
      #CI
      PHICIH[[t]][r1[]==1] <- apply(exp(-(mhH1[,s,ss, ,t]+mhW1[,s,ss, ,t])) ,2, function(x) quantile(x,prob=97.5/100,na.rm=T))
      #apply(mhW1[,s,ss, ,t],2, function(x) quantile(x,prob=2.5/100,na.rm=T))))
      PHICIH[[t]][r1[]==0] <- NA
      #CI
      PHICIL[[t]][r1[]==1] <-  apply(exp(-(mhH1[,s,ss, ,t]+mhW1[,s,ss, ,t])) ,2, function(x) quantile(x,prob=2.5/100,na.rm=T))
      PHICIL[[t]][r1[]==0] <- NA
      
      plot(PHI[[t]], breaks=cuts, col = col, main="",legend=F ,axes=F ,box=F) #p
      if(t==1 ){
        mtext("phi",2)
      } 
      if(t ==(nYears-1)){
        plot(PHI[[t]], legend.only=TRUE,breaks=cuts, col=col,
             legend.width = 2,
             axis.args=list(at=round(seq(0, 1, length.out = 5),digits = 1),
                            labels=round(seq(0, 1, length.out = 5),digits = 1), 
                            cex.axis=0.6),
             legend.args=list(text='prob', side=3, font=2, line=2.5, cex=0.8))
      }
    }
  }
}
# dev.off()


### ====      3.3.3.3 PREDICT MAPS WITH ABUNDANCE  ====
tmpW <- tmpH <- tmpPhi <- list()
Htmp <- Wtmp <- PHItmp <- list()
PointsPhi <- PointsH <- PointsW <- SY <- SX <- list()

c <- NULL
STATE <- c(2,3)

for(t in 1:(nYears-1)){
  PointsPhi[[t]] <- PointsH[[t]] <- PointsW[[t]] <- SY[[t]] <- SX[[t]] <- list()
  for(s in 1:2){
    PointsPhi[[t]][[s]] <- PointsH[[t]][[s]] <- PointsW[[t]][[s]] <- SY[[t]][[s]] <- SX[[t]][[s]] <- list()
    for(ss in 1:2){
      PointsPhi[[t]][[s]][[ss]] <- PointsH[[t]][[s]][[ss]]<- PointsW[[t]][[s]][[ss]] <- SY[[t]][[s]][[ss]] <- SX[[t]][[s]][[ss]] <- list()
      ##extract density #
      tmpRaster <- habbR#disaggregate(habbR,1)
      tmp <- mask(tmpRaster, mask = rasterToPolygons(habitat.rWthBuffer, fun=function(x) x==1))
      tmpPol <- rasterToPolygons(tmp, fun = function(x) x>0)
      
      #Extract survival/mort and so on 
      Htmp[[t]] <- Wtmp[[t]] <- PHItmp[[t]] <- r1 <- habitat.r
      #H
      Htmp[[t]][r1[]==1] <- colMeans(h[,s,ss, ,t])
      Htmp[[t]][r1[]==0] <- NA
      Htmp[[t]] <- mask(Htmp[[t]],mask = rasterToPolygons(habitat.rWthBuffer,fun=function(x) x==1))
      
      #w
      Wtmp[[t]][r1[]==1] <- colMeans(w[,s,ss, ,t])
      Wtmp[[t]][r1[]==0] <- NA
      Wtmp[[t]] <- mask(Wtmp[[t]],mask = rasterToPolygons(habitat.rWthBuffer,fun=function(x) x==1))
      
      #phi
      PHItmp[[t]][r1[]==1] <- exp(-(colMeans(mhH1[,s,ss, ,t])+colMeans(mhW1[,s,ss, ,t])))
      PHItmp[[t]][r1[]==0] <- NA
      PHItmp[[t]] <- mask(PHItmp[[t]], mask = rasterToPolygons(habitat.rWthBuffer, fun=function(x) x==1))
      
      
      ### GET SOME SXY FROM SOME ITERATIONS  
      iterar <- seq(1, dim(myResults$sims.list$sxy)[1],by=1)
      
      sxList <- syList <- list()
      #Subset per sex 
      ZF <- myResults$sims.list$z[iterar,,t]
      whichSex <- myResults$sims.list$sex[iterar,] %in% (s-1)
      ZF[!whichSex] <- 5
      
      for(i in 1:length(iterar)){
        sx = myResults$sims.list$sxy[iterar[i],,1,t]
        sy =  myResults$sims.list$sxy[iterar[i],,2,t]
        #z = ZF[i,]
        #Subset per state
        whichAlive <- which(ZF[i,] %in% c(STATE[ss]))
        
        sxList[[i]] <- sx[whichAlive]
        syList[[i]] <- sy[whichAlive]
      }
      
      SX[[t]][[s]][[ss]] <- do.call(base::c, sxList)
      SY[[t]][[s]][[ss]] <- do.call(base::c, syList)
      
      cells <- cellFromXY(object = PHItmp[[t]], xy = cbind(SX[[t]][[s]][[ss]], SY[[t]][[s]][[ss]]) )
      
      PointsPhi[[t]][[s]][[ss]] <- PHItmp[[t]][cells]#,length(Points[[i]])
      PointsH[[t]][[s]][[ss]]  <- Htmp[[t]][cells]#,length(Points[[i]])
      PointsW[[t]][[s]][[ss]] <- Wtmp[[t]][cells]#,length(Points[[i]])
    }
  }
}


#get color scale
rbPal <- colorRampPalette(c("white", "lightslateblue", "yellow", "orange", "red", "red"))


## assign color to each realized AC location depending the predicted mortality value
tmpW <- tmpH <- tmpPhi <-list()
tmpWCIH <- tmpHCIH <- tmpPhiCIH <-list()
tmpWCIL <- tmpHCIL <- tmpPhiCIL <-list()

PointsPhi <- PointsH <- PointsW <- SY <- SX <- list()
PointsPhiCIL <- PointsHCIL <- PointsWCIL <- list()
PointsPhiCIH <- PointsHCIH <- PointsWCIH <- list()
PointsPhiInterval <- PointsHInterval <- PointsWInterval <- list()


for(t in 1:(nYears-1)){
  
  tmpRaster <- habbR#disaggregate(SpatialRasterDensity[[t]],1)
  tmp <- mask(tmpRaster,mask = rasterToPolygons(habitat.rWthBuffer,fun=function(x) x==1))
  tmpPol <- rasterToPolygons(tmp,fun = function(x) x>0)
  
  #Overall mortality   
  tmpRasterPhi <- raster::disaggregate(PHI[[t]],2)
  tmpPhi[[t]] <- mask(tmpRasterPhi,mask = rasterToPolygons(habitat.rWthBuffer, fun=function(x) x==1))
  #tmpPHI <- tmpPhi[tmp[] > 0 | !is.na(tmp[]) ]
  #CI 
  tmpRasterPhiCIH <- raster::disaggregate(PHICIH[[t]],2)
  tmpPhiCIH[[t]] <- mask(tmpRasterPhiCIH,mask = rasterToPolygons(habitat.rWthBuffer, fun=function(x) x==1))
  tmpRasterPhiCIL <- raster::disaggregate(PHICIL[[t]],2)
  tmpPhiCIL[[t]] <- mask(tmpRasterPhiCIL,mask = rasterToPolygons(habitat.rWthBuffer, fun=function(x) x==1))
  
  
  #H legal mortality 
  tmpRasterH <- disaggregate(H[[t]],2)
  tmpH[[t]] <- mask(tmpRasterH,mask = rasterToPolygons(habitat.rWthBuffer,fun=function(x) x==1))
  #CI 
  tmpRasterHCIH <- raster::disaggregate(HCIH[[t]],2)
  tmpHCIH[[t]] <- mask(tmpRasterHCIH,mask = rasterToPolygons(habitat.rWthBuffer, fun=function(x) x==1))
  tmpRasterHCIL <- raster::disaggregate(HCIL[[t]],2)
  tmpHCIL[[t]] <- mask(tmpRasterHCIL,mask = rasterToPolygons(habitat.rWthBuffer, fun=function(x) x==1))
  
  #W other mortality
  tmpRasterW <- disaggregate(W[[t]],2)
  tmpW[[t]] <- mask(tmpRasterW,mask = rasterToPolygons(habitat.rWthBuffer,fun=function(x) x==1))
  #CI 
  tmpRasterWCIH <- raster::disaggregate(WCIH[[t]],2)
  tmpWCIH[[t]] <- mask(tmpRasterWCIH,mask = rasterToPolygons(habitat.rWthBuffer, fun=function(x) x==1))
  tmpRasterWCIL <- raster::disaggregate(WCIL[[t]],2)
  tmpWCIL[[t]] <- mask(tmpRasterWCIL,mask = rasterToPolygons(habitat.rWthBuffer, fun=function(x) x==1))
  
  
# chose for which iteration we extract the color 
  iterar <- seq(1,dim(myResults$sims.list$sxy)[1],by=1)
  sxList <- syList <- list()
  for(i in 1:length(iterar)){
    sx = myResults$sims.list$sxy[iterar[i],,1,t]#x value  
    sy =  myResults$sims.list$sxy[iterar[i],,2,t]#y value 
    z = myResults$sims.list$z[iterar[i],,t]# z
    whichAlive <- which(z %in% c(2,3))
    sxList[[i]] <- sx[whichAlive]
    syList[[i]] <- sy[whichAlive]
  }
  
  SX[[t]] <- do.call(base::c,sxList)
  SY[[t]] <- do.call(base::c,syList)
  
  #extract predicted value 
  cells <- cellFromXY(object = tmpPhi[[t]], xy = cbind(SX[[t]],SY[[t]]) )
  
  #Mean
  PointsPhi[[t]] <- tmpPhi[[t]][cells]#,length(Points[[i]])
  PointsH[[t]] <- tmpH[[t]][cells]#,length(Points[[i]])
  PointsW[[t]] <- tmpW[[t]][cells]#,length(Points[[i]])
  
  #low CI
  PointsPhiCIL[[t]] <- tmpPhiCIL[[t]][cells]#,length(Points[[i]])
  PointsHCIL[[t]] <- tmpHCIL[[t]][cells]#,length(Points[[i]])
  PointsWCIL[[t]] <- tmpWCIL[[t]][cells]#,length(Points[[i]])
  
  #High CI
  PointsPhiCIH[[t]] <- tmpPhiCIH[[t]][cells]#,length(Points[[i]])
  PointsHCIH[[t]] <- tmpHCIH[[t]][cells]#,length(Points[[i]])
  PointsWCIH[[t]] <- tmpWCIH[[t]][cells]#,length(Points[[i]])
  
  #CI width 
  PointsPhiInterval[[t]] <- tmpPhiCIH[[t]][cells]-tmpPhiCIL[[t]][cells]#,length(Points[[i]])
  PointsHInterval[[t]] <- tmpHCIH[[t]][cells]-tmpHCIL[[t]][cells]#,length(Points[[i]])
  PointsWInterval[[t]] <- tmpWCIH[[t]][cells]-tmpWCIL[[t]][cells]#,length(Points[[i]])
}

#PLOT for 
## DEFINE MAX VALUES 
MAX <- max(c(
  unlist(lapply(PointsW,function(x) max(x,na.rm=T))),
  unlist(lapply(PointsH,function(x) max(x,na.rm=T))),
  unlist(lapply(PointsPhi,function(x) max(1-x,na.rm=T)))))

CUT <- cut(seq(0, MAX,length.out=200), breaks = 200)

rbPal <- colorRampPalette(c("white", "slateblue", "yellow", "orange", "red", "red"))



### ====  3.4.PLOT ====
### ====    3.4.1. SURVIVAL ====
### ====      3.4.1.1 MEAN ====
t <- nYears-1 # plot the last year
begin <- 0.1
end <- 1
fillCountry <- c(grey(0.25),grey(0.20))


layout(matrix(1:2,ncol=2), width = rep(c(2,0.55),1),height = c(1,1))
par(mar=c(0,0,0,0))
#OVERALL MORTALITY 
plot(habitat.poly$geometry,border=NA)
rbPal <- colorRampPalette(c("white", "slateblue", "yellow", "orange", "red", "red"))
MAXPhi <- max(c(unlist(lapply(PointsPhi[t],function(x) max(1-x,na.rm=T)))))

sxy <- data.frame(cbind(SX[[t]],SY[[t]]))
colnames(sxy) <- c("x", "y")
sxy <-  st_as_sf(sxy, coords = c("x", "y"))
st_crs(sxy) <- st_crs(COUNTRIES)
CellsInBuffer  <- st_intersects( sxy, COUNTRIES,sparse = F)
                                   
  
cellsIn <-  which(!rowSums(CellsInBuffer)%in% 0)
Col <- rbPal(200)[as.numeric( cut(1-PointsPhi[[t]][cellsIn], seq(0,MAXPhi,length.out=200)))]
#PLOT COUNTIES 
plot(COUNTRIES$geometry,border=fillCountry,add=T,col=fillCountry)
#PLOT REALIZED SXY 
points(SY[[t]][cellsIn]~SX[[t]][cellsIn], col=adjustcolor(Col, alpha.f = 0.05), pch=16, cex=0.01)
plot(Linesf$geometry,add=T,col=adjustcolor("white",alpha.f = 0.8),lwd=0.1)
  
#legend
legend_image <- as.raster(matrix(rev(rbPal(20)), ncol=1))
par(mar=c(0,0,0,0))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=1.05, y = seq(0.25,0.65,l=5), labels = round(seq(0,1,l=5),digits=1),col="black",cex=2.8)
rasterImage(legend_image, xleft=0.1, xright= 0.5, ybottom=0.25, ytop= 0.65)
  
### ====      3.4.1.2 UNCERTAINTY CI WIDTH ====
layout(matrix(1:2,ncol=2), width = rep(c(2,0.55),1),height = c(1,1))
par(mar=c(0,0,0,0))
plot(habitat.poly$geometry, border=NA)
rbPal <- colorRampPalette(c("white", "slateblue", "yellow", "orange", "red", "red"))
MAXPhi <- max(c(
    unlist(lapply(PointsPhiInterval[t],function(x) max(x,na.rm=T))),
    unlist(lapply(PointsWInterval[t],function(x) max(x,na.rm=T))),
    unlist(lapply(PointsHInterval[t],function(x) max(x,na.rm=T)))))
  
minPhi <- max(c(
    unlist(lapply(PointsPhiInterval[t],function(x) min(x,na.rm=T))),
    unlist(lapply(PointsWInterval[t],function(x) min(x,na.rm=T))),
    unlist(lapply(PointsHInterval[t],function(x) min(x,na.rm=T)))))
  
sxy <- data.frame(cbind(SX[[t]],SY[[t]]))
colnames(sxy) <- c("x", "y")
sxy <-  st_as_sf(sxy, coords = c("x", "y"))
st_crs(sxy) <- st_crs(COUNTRIES)
CellsInBuffer  <- st_intersects( sxy, COUNTRIES,sparse = F)
cellsIn <-  which(!rowSums(CellsInBuffer)%in% 0)

Col <- viridis(200,begin=0.2)[as.numeric( cut(PointsPhiInterval[[t]][cellsIn], seq(minPhi,MAXPhi,length.out=200)))]
  
plot(COUNTRIES$geometry,border=fillCountry,add=T,col=fillCountry)
points(SY[[t]][cellsIn]~SX[[t]][cellsIn], col=adjustcolor(Col, alpha.f = 0.05), pch=16, cex=0.01)
plot(Linesf$geometry,add=T,col=adjustcolor("white",alpha.f = 0.8),lwd=0.1)
#legend
legend_image <- as.raster(matrix(rev(rbPal(20)), ncol=1))
legend_image <- as.raster(matrix(rev(viridis(20,begin=0.2)), ncol=1))
par(mar=c(0,0,0,0))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=1.05, y = seq(0.25,0.65,l=5), labels = round(seq(minPhi,MAXPhi,l=5),digits=1),col="black",cex=2.8)
rasterImage(legend_image, xleft=0.1, xright= 0.5, ybottom=0.25, ytop= 0.65)
  
  
### ====    3.4.2. OTHER CAUSES  ====
### ====      3.4.2.1 MEAN ====
MAXHW <- max(c(
  unlist(lapply(PointsW[t],function(x) max(x,na.rm=T))),
  unlist(lapply(PointsH[t],function(x) max(x,na.rm=T)))))
layout(matrix(1:2,ncol=2), width = rep(c(2,0.55),1),height = c(1,1))
par(mar=c(0,0,0,0))
plot(habitat.poly$geometry,border=NA)
Col <- rbPal(200)[as.numeric( cut(PointsW[[t]][cellsIn], seq(MAXPhi,0,length.out=200)))]
plot(COUNTRIES$geometry,border=fillCountry,add=T,col=fillCountry)
points(SY[[t]][cellsIn]~SX[[t]][cellsIn],col=adjustcolor(Col,alpha.f = 0.05),pch=16,cex=0.01)
plot(Linesf$geometry,add=T,col=adjustcolor("white",alpha.f = 0.8),lwd=0.1)

#legend
#legend_image <- as.raster(matrix((viridis(20,begin = begin,end=end)), ncol=1))
legend_image <- as.raster(matrix(rev(rbPal(20)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=1.05, y = seq(0.25,0.65,l=5), labels = round(seq(0,1,l=5),digits=1),col="black",cex=2.8)
rasterImage(legend_image, xleft=0.1, xright= 0.5, ybottom=0.25, ytop= 0.65)



### ====      3.4.2.2 CI WIDTH ====
layout(matrix(1:2,ncol=2), width = rep(c(2,0.55),1),height = c(1,1))
par(mar=c(0,0,0,0))
plot(habitat.poly$geometry, border=NA)
rbPal <- colorRampPalette(c("white", "slateblue", "yellow", "orange", "red", "red"))
MAXPhi <- max(c(unlist(lapply(PointsPhiInterval[t],function(x) max(x,na.rm=T))),
                unlist(lapply(PointsWInterval[t],function(x) max(x,na.rm=T))),
                unlist(lapply(PointsHInterval[t],function(x) max(x,na.rm=T)))))

minPhi <- max(c(unlist(lapply(PointsPhiInterval[t],function(x) min(x,na.rm=T))),
                unlist(lapply(PointsWInterval[t],function(x) min(x,na.rm=T))),
                unlist(lapply(PointsHInterval[t],function(x) min(x,na.rm=T)))))

sxy <- data.frame(cbind(SX[[t]],SY[[t]]))
colnames(sxy) <- c("x", "y")
sxy <-  st_as_sf(sxy, coords = c("x", "y"))
st_crs(sxy) <- st_crs(COUNTRIES)
CellsInBuffer  <- st_intersects( sxy, COUNTRIES,sparse = F)
cellsIn <-  which(!rowSums(CellsInBuffer)%in% 0)

Col <- viridis(200,begin=0.2)[as.numeric( cut(PointsWInterval[[t]][cellsIn], seq(minPhi,MAXPhi,length.out=200)))]

plot(COUNTRIES$geometry,border=fillCountry,add=T,col=fillCountry)
points(SY[[t]][cellsIn]~SX[[t]][cellsIn], col=adjustcolor(Col, alpha.f = 0.05), pch=16, cex=0.01)
plot(Linesf$geometry,add=T,col=adjustcolor("white",alpha.f = 0.8),lwd=0.1)

#legend
legend_image <- as.raster(matrix(rev(viridis(20,begin=0.2)), ncol=1))

par(mar=c(0,0,0,0))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=1.05, y = seq(0.25,0.65,l=5), labels = round(seq(minPhi,MAXPhi,l=5),digits=1),col="black",cex=2.8)
rasterImage(legend_image, xleft=0.1, xright= 0.5, ybottom=0.25, ytop= 0.65)


  
### ====    3.4.3. LEGAL MORTALITY  ====
### ====      3.4.3.1 MEAN ====
layout(matrix(1:2,ncol=2), width = rep(c(2,0.55),1),height = c(1,1))
par(mar=c(0,0,0,0))
plot(habitat.poly$geometry,border=NA)
Col <- rbPal(200)[as.numeric( cut(PointsH[[t]][cellsIn], seq(MAXPhi,0,length.out=200)))]
plot(COUNTRIES$geometry,border=fillCountry,add=T,col=fillCountry)
points(SY[[t]][cellsIn]~SX[[t]][cellsIn],col=adjustcolor(Col,alpha.f = 0.05),pch=16,cex=0.01)
plot(Linesf$geometry,add=T,col=adjustcolor("white",alpha.f = 0.8),lwd=0.1)
  
#legend
legend_image <- as.raster(matrix(rev(rbPal(20)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=1.05, y = seq(0.25,0.65,l=5), labels = round(seq(0,1,l=5),digits=1),col="black",cex=2.8)
rasterImage(legend_image, xleft=0.1, xright= 0.5, ybottom=0.25, ytop= 0.65)
  

  
### ====      3.4.3.2 CI WIDTH ====
layout(matrix(1:2,ncol=2), width = rep(c(2,0.55),1),height = c(1,1))
par(mar=c(0,0,0,0))
plot(habitat.poly$geometry,border=NA)
rbPal <- colorRampPalette(c("white", "slateblue", "yellow", "orange", "red", "red"))
MAXPhi <- max(c(unlist(lapply(PointsPhiInterval[t],function(x) max(x,na.rm=T))),
                  unlist(lapply(PointsWInterval[t],function(x) max(x,na.rm=T))),
                  unlist(lapply(PointsHInterval[t],function(x) max(x,na.rm=T)))))
  
minPhi <- max(c(unlist(lapply(PointsPhiInterval[t],function(x) min(x,na.rm=T))),
                  unlist(lapply(PointsWInterval[t],function(x) min(x,na.rm=T))),
                  unlist(lapply(PointsHInterval[t],function(x) min(x,na.rm=T)))))
  
sxy <- data.frame(cbind(SX[[t]],SY[[t]]))
colnames(sxy) <- c("x", "y")
sxy <-  st_as_sf(sxy, coords = c("x", "y"))
st_crs(sxy) <- st_crs(COUNTRIES)
CellsInBuffer  <- st_intersects( sxy, COUNTRIES,sparse = F)
cellsIn <-  which(!rowSums(CellsInBuffer)%in% 0)
Col <- viridis(200,begin=0.2)[as.numeric( cut(PointsHInterval[[t]][cellsIn], seq(0,1,length.out=200)))]
plot(COUNTRIES$geometry,border=fillCountry,add=T,col=fillCountry)
points(SY[[t]][cellsIn]~SX[[t]][cellsIn], col=adjustcolor(Col, alpha.f = 0.05), pch=16, cex=0.01)
plot(Linesf$geometry,add=T,col=adjustcolor("white",alpha.f = 0.8),lwd=0.1)
#legend
legend_image <- as.raster(matrix(rev(viridis(20,begin=0.2)), ncol=1))
par(mar=c(0,0,0,0))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=1.05, y = seq(0.25,0.65,l=5), labels = round(seq(minPhi,MAXPhi,l=5),digits=1),col="black",cex=2.8)
rasterImage(legend_image, xleft=0.1, xright= 0.5, ybottom=0.25, ytop= 0.65)

