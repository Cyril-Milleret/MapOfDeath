
#### 1.Density function ####
dCatZwhBhPsiHRSxy1Dens <- nimbleFunction(run = function( x = double(0)
                                                     , zit = double(0)
                                                     , gamma = double(0)
                                                     , mhH = double(1)
                                                     , mhW = double(1)
                                                     , betaH = double(1)
                                                     , betaW = double(1)
                                                     , betaHDens = double(0)
                                                     , betaWDens = double(0)
                                                     , psi = double(0)
                                                     , sxy = double(1)
                                                     , habCovsH = double(2)
                                                     , habCovsW = double(2)
                                                     , dens = double(1)
                                                     , habIndex = double(2)
                                                     , ResizeFactor = double(0)
                                                     , log = integer(0, default = 0)){
  # Return type declaration
  returnType(double(0))
  
  
  
  if(zit == 1){
    logLikelihood <- dcat(x, prob = c(1 - gamma, gamma), log=1)
    if(log == 1){return(logLikelihood)}else{return(exp(logLikelihood))}
  }
  
  if(zit == 2){
    sxyID <- habIndex[trunc(sxy[2]/ResizeFactor)+1, trunc(sxy[1]/ResizeFactor)+1]
    IndCovH <- habCovsH[sxyID,]
    IndCovW <- habCovsW[sxyID,]    
    densiCov <- dens[sxyID]

    
    mhH1 <- exp(mhH[1] + inprod(betaH, IndCovH) + betaHDens * densiCov)
    mhW1 <- exp(mhW[1] + inprod(betaW, IndCovW) + betaWDens * densiCov)
    
    # phi <- exp(-(mhH1+mhW))
    # h <- (1-phi)* (mhH1/(mhH1+mhW))
    # w <- (1-phi)* (mhW/(mhH1+mhW))
    
    h <- (1-exp(-(mhH1+mhW1)))* (mhH1/(mhH1+mhW1))
    w <- (1-exp(-(mhH1+mhW1)))* (mhW1/(mhH1+mhW1))
    phi <- 1-h-w
    # h1 <-  ilogit(logit(h) + IndCov*betaH)
    # phi <- 1-h1-w
    # # Constrain so phi never gets negative 
    # if(phi<0){
    #   if(log == 0) return(0.0)
    #   else return(-Inf)}
    
    logLikelihood <- dcat(x, prob = c(0, phi*(1-psi), phi*(psi), h, w), log=1)
    if(log == 1){return(logLikelihood)}else{return(exp(logLikelihood))}
  }
  
  if(zit == 3){
    sxyID <- habIndex[trunc(sxy[2]/ResizeFactor)+1, trunc(sxy[1]/ResizeFactor)+1]
    IndCovH <- habCovsH[sxyID,]
    IndCovW <- habCovsW[sxyID,]
    densiCov <- dens[sxyID]
    
    mhH1 <- exp(mhH[2] + inprod(betaH, IndCovH) + betaHDens * densiCov)
    mhW1 <- exp(mhW[2] + inprod(betaW, IndCovW) + betaWDens * densiCov)
    
    # phi <- exp(-(mhH1+mhW))
    # h <- (1-phi)* (mhH1/(mhH1+mhW))
    # w <- (1-phi)* (mhW/(mhH1+mhW))
    
    h <- (1-exp(-(mhH1+mhW1)))* (mhH1/(mhH1+mhW1))
    w <- (1-exp(-(mhH1+mhW1)))* (mhW1/(mhH1+mhW1))
    phi <- 1-h-w
    # h1 <-  ilogit(logit(h) + IndCov*betaH)
    # phi <- 1-h1-w
    # # Constrain so phi never gets negative 
    # if(phi<0){
    #   if(log == 0) return(0.0)
    #   else return(-Inf)}
    
    logLikelihood <- dcat(x, prob = c(0, 0, phi, h, w), log=1)
    if(log == 1){return(logLikelihood)}else{return(exp(logLikelihood))}
  }
  
  if(zit == 4 | zit == 5 ){
    logLikelihood <- dcat(x, prob = c(0, 0, 0, 0, 1), log=1)
    if(log == 1){return(logLikelihood)}else{return(exp(logLikelihood))}
  }
  
  
})

#### 2.Sampling function ####
rCatZwhBhPsiHRSxy1Dens <- nimbleFunction(run = function( n = integer(0)
                                                     , zit = double(0)
                                                     , gamma = double(0)
                                                     , mhH = double(1)
                                                     , mhW = double(1)
                                                     , betaH = double(1)
                                                     , betaW = double(1)
                                                     , betaHDens = double(0)
                                                     , betaWDens = double(0)
                                                     , psi = double(0)
                                                     , sxy = double(1)
                                                     , habCovsH = double(2)
                                                     , habCovsW = double(2)
                                                     , dens = double(1)
                                                     , habIndex = double(2)
                                                     , ResizeFactor = double(0)
){
  # Return type declaration
  returnType(double(0))
  
  ## EXTRACT LOCATION OF THE ID
  
  
  
  ## EXTRACT LOCATION OF THE ID
  # phi1 <-  ilogit(logit(phi[index]) + IndCov*betaPhiInd)
  # w1 <-  ilogit(logit(w[index]) + IndCov*betawInd)
  # h1 <-  ilogit(logit(h[index]) + IndCov*betaHInd)
  
  if(zit == 1){
    state <- rcat(1, prob = c(1 - gamma, gamma))
    return(state)
  }
  
  if(zit == 2){
    sxyID <- habIndex[trunc(sxy[2]/ResizeFactor)+1, trunc(sxy[1]/ResizeFactor)+1]
    IndCovH <- habCovsH[sxyID,]
    IndCovW <- habCovsW[sxyID,]
    densiCov <- dens[sxyID]
    
    mhH1 <- exp(mhH[1] + inprod(betaH, IndCovH)+ betaHDens * densiCov)
    mhW1 <- exp(mhW[1] + inprod(betaW, IndCovW)+ betaWDens * densiCov)
    
    # phi <- exp(-(mhH1+mhW))
    # h <- (1-phi)* (mhH1/(mhH1+mhW))
    # w <- (1-phi)* (mhW/(mhH1+mhW))
    
    h <- (1-exp(-(mhH1+mhW1)))* (mhH1/(mhH1+mhW1))
    w <- (1-exp(-(mhH1+mhW1)))* (mhW1/(mhH1+mhW1))
    phi <- 1-h-w
    
    # h <- (1-exp(-(mhH1+mhW)))* (mhH1/(mhH1+mhW))
    # w <- (1-exp(-(mhH1+mhW)))* (mhW/(mhH1+mhW))
    
    ## EXTRACT LOCATION OF THE ID
    #phi <- 1-h-w
    state <- rcat(1, prob = c(0, phi*(1-psi), phi*(psi), h, w) )
    return(state)
  }
  
  if(zit == 3){
    sxyID <- habIndex[trunc(sxy[2]/ResizeFactor)+1, trunc(sxy[1]/ResizeFactor)+1]
    IndCovH <- habCovsH[sxyID,]
    IndCovW <- habCovsW[sxyID,]
    densiCov <- dens[sxyID]
    
    mhH1 <- exp(mhH[2] + inprod(betaH, IndCovH)+ betaHDens * densiCov)
    mhW1 <- exp(mhW[2] + inprod(betaW, IndCovW)+ betaWDens * densiCov)
    
    
    # phi <- exp(-(mhH1+mhW))
    # h <- (1-phi)* (mhH1/(mhH1+mhW))
    # w <- (1-phi)* (mhW/(mhH1+mhW))
    
    h <- (1-exp(-(mhH1+mhW1)))* (mhH1/(mhH1+mhW1))
    w <- (1-exp(-(mhH1+mhW1)))* (mhW1/(mhH1+mhW1))
    phi <- 1-h-w
    
    # h <- (1-exp(-(mhH1+mhW)))* (mhH1/(mhH1+mhW))
    # w <- (1-exp(-(mhH1+mhW)))* (mhW/(mhH1+mhW))
    
    ## EXTRACT LOCATION OF THE ID
    #phi <- 1-h-w
    state <- rcat(1, prob = c(0, 0, phi*(psi), h, w))
    return(state)
  }
  
  if(zit == 4 | zit == 5 ){
    state <- 5
    return(state)
  }
  
})


#### 3.Registration ####
registerDistributions(list(
  dCatZwhBhPsiHRSxy1Dens = list(
    BUGSdist = "dCatZwhBhPsiHRSxy1Dens(zit, gamma, mhH, mhW, betaH, betaW,betaHDens, betaWDens, psi, sxy, habCovsH, habCovsW, dens, habIndex, ResizeFactor)",
    types = c( "value = double(0)", "zit = double(0)", "gamma = double(0)", "mhH = double(1)"
               , "mhW = double(1)","betaH = double(1)", "betaW = double(1)","betaHDens = double(0)"
               , "betaWDens = double(0)",  "psi = double(0)",
               "sxy = double(1)", "habCovsH = double(2)", "dens = double(1)",
               "habCovsW = double(2)","habIndex = double(2)","ResizeFactor = double(0)"
    ),
    pqAvail = FALSE,
    discrete = TRUE
  )))

