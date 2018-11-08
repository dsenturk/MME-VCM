MMEVCM_estimation <- function(data,       # data.frame in long format with nine labeled columns (described below)
                                          # and row length equal to the length of the vectorized observations across all 
                                          # subjects and facilities
                                          # DATA.FRAME COLUMNS: 
                                          # fid: facility IDs (vector of length sum(Nij))
                                          # sid: subject IDs (vector of length sum(Nij))
                                          # y: hospitalization outcome data (vector of length sum(Nij))
                                          # t: follow-up time (vector of length sum(Nij)) 
                                          # z1: facility-level covariate (vector of length sum(Nij))
                                          # z2: facility-level covariate (vector of length sum(Nij))
                                          # x1: subject-level covariate (vector of length sum(Nij))
                                          # x2: subject-level covariate (vector of length sum(Nij))
                                          # x0: a vector of 1s to add the intercept term (vector of length sum(Nij)) 
                             
                             bwThetaBeta  # bandwidth for estimating theta(t) and beta(t) (scalar)
                             ){
  
  #############################################################################
  ## Description: Function for estimation of MME-VCM model parameters described in "A Multilevel Mixed Effects Varying Coeffcient Model 
  ##              with Multilevel Predictors and Random Effects for Modeling Hospitalization Risk in Patients on Dialysis", including estimation 
  ##              and inference of time-varying effects of multilevel risk factors and variance of multilevel random effects. 
  ## Args:        see above
  ## Returns:     list()
  ##              theta: estimated facility-level risk factor effects (matrix of dimension 2*20)
  ##              beta: estimated subject-level risk factor effects and the intercept term (matrix of dimension 3*20) 
  ##              sigma2b: estimated variance of subject specific random effects (scalar)
  ##              sigma2gamma: estimated variance of facility specific random effects (scalar)
  ##              thetaSE: standard error of estimated facility-level risk factor effects (matrix of dimension 2*20)
  ##              betaSE: standard error of estimated subject-level risk factor effects and the intercept term (matrix of dimension 3*20)
  ##              sigma2bSE: standard error of estimated variance of subject specific random effects (scalar)
  ##              sigma2gammaSE: standard error of estimated variance of facility specific random effects (scalar)
  ##              bijEst: estimated posterior mean of subject-level random effects (vector of length sum(Ni))
  ##              gammaiEst: estimated posterior mean of facility-level random effects (vector of length numF)
  ##              grid: grid points used to estimate the varying coefficient functions gammai(t), theta(t) and beta(t) (vector of length 20)
  #############################################################################
  
  # Install missing packages
  list.of.packages <- c("MASS", "statmod", "mvtnorm", "bisoreg", "lme4")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages) 
  
  # Load packages  
  library(MASS)
  library(statmod)
  library(mvtnorm)
  library(bisoreg) 
  library(lme4)

  # Define functions 
  invLogit <- function(x) {
    return(exp(x)/(1+exp(x)))
  }
  getArea <- function(x,y) {
    if(length(x) != length(y)) {
      stop("The number of elements in the domain and range must equal.")
    }
    if(sum(sort(x) != x) > 0) {
      stop("The endpoints of the trapezoids must be ordered.")
    }
    k <- length(x)
    return(1/2*sum((y[-k]+y[-1])*(x[-1]-x[-k])))
  }
  numuniq <- function(x){
    return(length(unique(x)))
  }
  
  
  # Create a list for returns as outlined above
  MMEVCMEst <- list(11)
  
  # Number of subject-level risk factors (including the intercept term)
  nbeta <- 3
  
  # Number of facility-level risk factors
  ntheta <- 2
  
  # Empty matrix for constructing Hessian matrix of theta(t) and beta(t) used for estimation step 4 in Section 2.2
  h <- diag(nbeta+ntheta)   
  
  # Define the grid points used to estimate the varying coefficient functions theta(t) and beta(t)
  gridPoints <- seq(0,1,1/19)
  ngrid <- length(gridPoints)
  
  # Maximum number of iterations allowed for convergence
  maxNRIter <- 50
  
  # Read data from input
  df <- data

  # Index for updating varying coefficient functions
  df$r <- round(19*df$t)+1 
  
  numF <- max(df$fid) # Number of facilities
  sumNi <- max(df$sid) # Number of subjects
  numOB <- as.numeric(table(df$sid)) # Number of observations per subject 
  numSubPF <- aggregate(df$sid,by=list(df$fid),FUN=numuniq)[,2]
  Ni <- rep(1:numF,numSubPF)
  numDisPF <- as.numeric(table(df$fid)) 
  Nij <- rep(numOB,numOB)  
  df$numOB <- Nij

  
  ###########################################################################
  # Implement the approximate EM algorithm as described in Section 2.2
  ###########################################################################
  
  # Find initial values for theta(t) and beta(t) as well as the variances of multilevel random effects
  # from the non-time-varying risk effects generalized linear model
  gm1b <- glmer(y ~ x1 + x2 + x0 + z1 + z2 + (1 | fid) + (1 | sid) - 1, data = df, family = binomial, nAGQ = 1)
  
  IniSig <- getME(gm1b,"theta")^2
  IniThetaBeta <- getME(gm1b,"beta") 
  thetaEsts <- IniThetaBeta[(nbeta+1):(ntheta+nbeta)]
  betaEsts <- IniThetaBeta[1:nbeta]
  # End of finding initial values for theta(t) and beta(t)
  
  # Implement estimation steps as described in Section 2.2
  # Step 1: set initial values for fitting MME-VCM
  betaVCMEst <- matrix(betaEsts,nbeta,ngrid)
  betaVCMEst1 <- matrix(0,nbeta,ngrid)
  thetaVCMEst <- matrix(thetaEsts,ntheta,ngrid)
  thetaVCMEst1 <- matrix(0,ntheta,ngrid) 
  sigma2bEst <- IniSig[1]
  sigmagammaEst <- IniSig[2]
  df$gammaEst <- 0
  df$bisEst <- 0
  
  numIter <- 0
  Diff <- 1
  while(Diff > 0.001 & numIter <= maxNRIter) {
    print(paste("Iter",numIter)) # Print number of iterations
    oldlpred <- rowSums(t(thetaVCMEst)[df$r,]*df[,5:(4+ntheta)]) + rowSums(t(betaVCMEst)[df$r,]*df[,(5+ntheta):(4+ntheta+nbeta)])
    oldObj <- df$gammaEst + df$bisEst + oldlpred
    oldPijks <- invLogit(oldObj) # P0ijk from the previous iteration
    
    # Step 2 (E-step): update posterior means and variances of multilevel random effects using fully exponential Laplace approximation
    meangamma.C <- c()
    vargamma.C <- c()
    meanb.C <- c()
    varb.C <- c()
    corbgamma <- c()
    Tol <- 1e-8
    maxLapIter <- 30
    for(fi in 1:numF){
      #Facility ID
      facid <- fi
      ##### Start fully exponential Laplce approximation #####
      FacData <- df[df$fid==facid,] # Data within a facility
      numOBFac <- aggregate(FacData$numOB, by=list(FacData$sid), FUN=mean)[,2]
      lpredFac1 <- oldlpred[df$fid==facid]
      lpredFac <- FacData$gammaEst + FacData$bisEst + lpredFac1
      pijksFac <- invLogit(lpredFac)
      qijksFac <- 1 - pijksFac 
      # Obtain the posterior mode of random effects in Laplace aproximation using Newton-Raphson
      # Objective function: need to minimize
      l.rand <- function(randF){
        randpred <- rep(randF[-length(randF)], numOBFac) + randF[length(randF)]
        l1 <- sum(randpred * FacData$y) - sum(log(1 + exp(lpredFac1 + randpred)))
        return(-l1 + sum(randF[-length(randF)]^2)/2/sigma2bEst + 
                 randF[length(randF)]^2/2/sigmagammaEst)
      }
      uij0 <- rep(0,numSubPF[facid])
      ui0 <- 0
      facdiff <- 1
      LapIter <- 0
      while(facdiff > Tol & LapIter <= maxLapIter){
        olduij <- uij0
        oldui <- ui0
        OldFunVal <- l.rand(c(olduij, oldui))
        #Gradient
        L <- aggregate((FacData$y-pijksFac), by=list(FacData$sid), FUN=sum)[,2]
        LL <- c((-L+uij0/sigma2bEst), -sum(L)+ui0/sigmagammaEst)
        #Hessian
        B <- aggregate((pijksFac*qijksFac), by=list(FacData$sid), FUN=sum)[,2]
        D <- sum(B) + 1/sigmagammaEst
        A <- B + 1/sigma2bEst
        E <- D - sum(B^2/A)
        FF <- B/A
        invHess <- FF %*% t(FF) / E + diag(1/A)
        invHess <- cbind(invHess, -FF/E)
        invHess <- rbind(invHess, c((-FF/E),1/E))
        updatedir <- invHess %*% LL
        #line search
        for(Laps in 1:10){
          s <- 0.5 ^ (Laps-1)
          uij0 <- olduij - s * updatedir[1:numSubPF[facid]]
          ui0 <- oldui - s * updatedir[numSubPF[facid]+1]
          if(l.rand(c(uij0, ui0)) < OldFunVal) break
        }
        facdiff <- max(abs(c(uij0-olduij,oldui-ui0)))
        # print(facdiff)
        lpredFac <- ui0 + rep(uij0,numOBFac) + lpredFac1
        pijksFac <- invLogit(lpredFac)
        qijksFac <- 1 - pijksFac 
        LapIter <- LapIter + 1
      }
      # Start computing correction terms in fully exponential Laplace
      # Non zero entries in Aij from Supplement page 2
      C.1 <- aggregate((pijksFac*qijksFac^2 - pijksFac^2*qijksFac), by=list(FacData$sid), FUN=sum)[,2]
      # Non zero entries in Bij from Supplement page 2
      C.2 <- aggregate((pijksFac*qijksFac^3 - 4*pijksFac^2*qijksFac^2 + pijksFac^3*qijksFac), by=list(FacData$sid), FUN=sum)[,2]
      # Combine the correction term for gamma
      C.end <- c(C.1,sum(C.1))
      C.end.TTT <- invHess %*% C.end
      C.TTT <- matrix(0, numSubPF[facid]+1,numSubPF[facid]+1)
      for(subj in 1:numSubPF[facid]){
        # Construct columns of Aij
        C.subvec <- rep(0, numSubPF[facid] + 1)
        C.subvec[numSubPF[facid] + 1] <- C.1[subj]
        C.subvec[subj] <- C.1[subj]
        # Elements in derivative of Sigma wrt c from Supplement (2)
        C.TTT[subj,] <- invHess %*% C.subvec
      }
      C.comb <- NULL
      C.var.comb <- NULL
      C.cov.comb <- NULL
      # Correction terms for gamma 
      subjj <- numSubPF[facid] + 1
      C.mat1 <- diag(C.TTT[,subjj]) # Derivative of Sigma wrt c from Supplement (2)
      C.mat1[numSubPF[facid]+1,] <- C.TTT[,subjj]
      C.mat1[,numSubPF[facid]+1] <- C.TTT[,subjj]
      C.mat1[numSubPF[facid]+1,numSubPF[facid]+1] <- C.end.TTT[subjj]
      # Correction term for gamma posterior mean, trace(V) in main paper (3)
      C.comb[subjj] <- sum(diag(invHess %*% C.mat1))
      # Correction term 1 for gamma posterior variance, trace(VV^T) in main paper (3)
      C.var.comb[subjj] <- sum(diag(invHess %*% C.mat1 %*% C.mat1 %*% invHess))
      C.matF <- C.mat1
      # Correction terms for bij
      for(subjj in 1:numSubPF[facid]){
        C.mat1 <- diag(C.TTT[,subjj]) # Derivative of Sigma wrt c from Supplement (2)
        C.mat1[numSubPF[facid]+1,] <- C.TTT[,subjj]
        C.mat1[,numSubPF[facid]+1] <- C.TTT[,subjj]
        C.mat1[numSubPF[facid]+1,numSubPF[facid]+1] <- C.end.TTT[subjj]
        #Correction term for bij posterior mean, trace(V) in main paper (3)
        C.comb[subjj] <- sum(diag(invHess %*% C.mat1))
        #Correction term 1 for bij posterior variance, trace(VV^T) in main paper (3)
        C.var.comb[subjj] <- sum(diag(invHess %*% C.mat1 %*% C.mat1 %*% invHess))
        #Correction term 1 for posterior covariance, trace(VV^T) in main paper (3)
        C.cov.comb[subjj] <- sum(diag(invHess %*% C.mat1 %*% C.matF %*% invHess))
      }
      #Correction term 2 for posterior variance and covariance
      C.var.mat <- list()
      C.var.mat.2 <- list()
      C.var.TTT.2 <- matrix(0, numSubPF[facid]+1,numSubPF[facid]+1)
      for(subj2 in 1:numSubPF[facid]){
        # Second derivative of Sigma wrt ui for computing the first term in the second derivative of Sigma wrt c from Supplement (2)
        C.var.TTT <- matrix(0, numSubPF[facid]+1,numSubPF[facid]+1)
        C.var.TTT[subj2, subj2] <- C.2[subj2]
        C.var.TTT[subj2, numSubPF[facid] + 1] <- C.2[subj2]
        C.var.TTT[numSubPF[facid] + 1, subj2] <- C.2[subj2]
        C.var.TTT[numSubPF[facid] + 1, numSubPF[facid] + 1] <- C.2[subj2]
        C.var.mat[[subj2]] <- C.var.TTT
        # First derivative of Sigma wrt ui for computing the second term in the second derivative of Sigma wrt c from Supplement (2)
        C.var.TTT[subj2, subj2] <- C.1[subj2]
        C.var.TTT[subj2, numSubPF[facid] + 1] <- C.1[subj2]
        C.var.TTT[numSubPF[facid] + 1, subj2] <- C.1[subj2]
        C.var.TTT[numSubPF[facid] + 1, numSubPF[facid] + 1] <- C.1[subj2]
        C.var.mat.2[[subj2]] <- C.var.TTT
        C.var.TTT.2 <- C.var.TTT.2 + C.var.TTT
      }
      C.var.mat.2[[numSubPF[facid] + 1]] <- C.var.TTT.2
      # Compute the first term in the second derivative of Sigma wrt c from Supplement (2)
      C.var.2 <- NULL 
      C.cov.2 <- NULL
      for(subj2 in 1:(numSubPF[facid] + 1)){
        ttt <- (invHess[subj2,] + invHess[subj2, numSubPF[facid]+1])^2
        C.var.mat2 <- matrix(0, numSubPF[facid]+1,numSubPF[facid]+1)
        C.cov.mat2 <- matrix(0, numSubPF[facid]+1,numSubPF[facid]+1)
        for(subjj2 in 1:numSubPF[facid]){
          C.var.mat2 <- C.var.mat2 + C.var.mat[[subjj2]] * ttt[subjj2]
          ttt2 <- (invHess[subj2,subjj2] + invHess[subj2, numSubPF[facid] + 1]) * (invHess[numSubPF[facid]+1, subjj] + invHess[numSubPF[facid]+1,numSubPF[facid]+1])
          C.cov.mat2 <- C.cov.mat2 + C.var.mat[[subjj2]] * ttt2
        }
        # First term in correction term 2 for posterior variance 
        C.var.2[subj2] <- sum(diag(invHess %*% C.var.mat2))
        C.cov.2[subj2] <- sum(diag(invHess %*% C.cov.mat2))
      }
      # Compute the second term in the second derivative of Sigma wrt c from Supplement (2)
      C.var.2.2 <- NULL
      C.cov.2.2 <- NULL
      for(subj2 in 1:(numSubPF[facid] + 1)){
        J.mat <- -invHess %*% C.var.mat.2[[subj2]] %*% invHess %*% invHess
        C.var.mat2.2 <- matrix(0, numSubPF[facid]+1,numSubPF[facid]+1)
        C.cov.mat2.2 <- matrix(0, numSubPF[facid]+1,numSubPF[facid]+1)
        for(subjj2 in 1:(numSubPF[facid]+1)){
          C.var.mat2.2 <- C.var.mat2.2 + C.var.mat.2[[subjj2]] * J.mat[subj2, subjj2]
          C.cov.mat2.2 <- C.cov.mat2.2 + C.var.mat.2[[subjj2]] * J.mat[numSubPF[facid] + 1, subjj2]
        }
        # Second term in correction term 2 for posterior variance 
        C.var.2.2[subj2] <- sum(diag(invHess %*% C.var.mat2.2))
        C.cov.2.2[subj2] <- sum(diag(invHess %*% C.cov.mat2.2))
      }
      # Combine the correction terms for posterior variance and covariance
      C.var.comb.2 <- C.var.2 + C.var.2.2
      C.cov.comb.2 <- (C.cov.2 + C.cov.2.2)[1:numSubPF[facid]]
      # Apply correction to posterior mean as in main paper (3)
      uij0.C <- uij0 - 0.5 * C.comb[1:numSubPF[facid]]
      ui0.C <- ui0 - 0.5 * C.comb[numSubPF[facid] + 1]
      # Apply correction to posterior variance and covariance as in main paper (3)
      # Posterior variance of bij
      vij0.C <- diag(invHess)[1:numSubPF[facid]] + 0.5 * C.var.comb[1:numSubPF[facid]] - 0.5 * C.var.comb.2[1:numSubPF[facid]]
      # Posterior variance of gamma
      vi0.C <- invHess[numSubPF[facid]+1,numSubPF[facid]+1] + 0.5 * C.var.comb[numSubPF[facid] + 1] - 0.5 * C.var.comb.2[numSubPF[facid] + 1]
      # Posterior covariance
      rij0.C <- invHess[1:numSubPF[facid], numSubPF[facid]+1] + 0.5 * C.cov.comb - 0.5 * C.cov.comb.2
      ##### End of fully exponential Laplce approximation #####
      meangamma.C <- c(meangamma.C, ui0.C)
      vargamma.C <- c(vargamma.C, vi0.C)
      meanb.C <- c(meanb.C, uij0.C)
      varb.C <- c(varb.C, vij0.C)    
      corbgamma <- c(corbgamma, rij0.C)
    }
    df$bisEst <- rep(meanb.C, numOB) # posterior mean of subject specific random effects
    df$gammaEst <- rep(meangamma.C, numDisPF) # posterior mean of facility specific random effects
    df$vb <- rep(varb.C, numOB) # posterior variance of subject specific random effects
    df$vgamma <- rep(vargamma.C, numDisPF) # posterior variance of facility specific random effects
    df$rbgamma <- rep(corbgamma, numOB) # posterior covariance between subject and facility random effects
    
    # Step 3 (M-step): update sigma2b and sigmagamma by maximizing the incomplete log-likelihood
    sigmagammaEst <- mean(meangamma.C^2 + vargamma.C)
    sigma2bEst <- mean(meanb.C^2 + varb.C)
    
    # Step 4 (M-step): update theta(t) and beta(t) by maximizing the approximate expected local log-likelihood
    Phi_OneStep <- matrix(0,ngrid,2*(ntheta+nbeta))
    for(j in 1:ngrid) {
      to <- gridPoints[j]
      indices1 <- (abs(df$t-to) < bwThetaBeta)
      dfTemp <- df[indices1,]  # Data used for local MLE
      tDiff <- dfTemp$t - to
      lpredTemp <- colSums(thetaVCMEst[,j]*t(dfTemp[,5:(4+ntheta)])) + colSums(betaVCMEst[,j]*t(dfTemp[,(5+ntheta):(4+ntheta+nbeta)])) + 
        colSums(thetaVCMEst1[,j]*t(dfTemp[,5:(4+ntheta)]*tDiff)) + colSums(betaVCMEst1[,j]*t(dfTemp[,(5+ntheta):(4+ntheta+nbeta)]*tDiff))
      obj <- dfTemp$gammaEst + dfTemp$bisEst + lpredTemp
      pijks <- invLogit(obj)
      qijks <- 1 - pijks
      epan <- tDiff/bwThetaBeta
      kernels <- 0.75*(1-epan^2)
      kernels <- (kernels > 0)*kernels/bwThetaBeta  # Kernel function
      a1 <- (dfTemp$vb + dfTemp$vgamma + 2 * dfTemp$rbgamma)/2*(pijks*qijks^2 - pijks^2*qijks)
      a2 <- pijks*qijks + (dfTemp$vb + dfTemp$vgamma + 2 * dfTemp$rbgamma)/2*(pijks*qijks^3 - 4*pijks^2*qijks^2 + pijks^3*qijks)
      # Construct Gradient
      zxMatrix <- data.matrix(cbind(dfTemp[,5:(4+ntheta+nbeta)],dfTemp[,5:(4+ntheta+nbeta)]*tDiff))
      zxVector <- t(zxMatrix) %*% ((dfTemp$y - pijks - a1)*kernels)
      # Construct Hessian
      zzxxMatrix <- t(as.matrix(zxMatrix))      
      for(mIndex in 1:(2*(ntheta+nbeta))) {
        zzxxMatrix[mIndex,] <- zzxxMatrix[mIndex,]*a2*kernels
      }
      zzxxMatrix <- zzxxMatrix %*% zxMatrix
      # Newton-Raphson updating direction
      invHessian <- solve(zzxxMatrix)
      Phi_OneStep[j,] <- invHessian %*% zxVector
    } 
    # Begin line search for theta(t) and beta(t) as described in Web Appendix B
    # Calculate global log-likelihood
    OldbetaVCMEst <- betaVCMEst
    OldbetaVCMEst1 <- betaVCMEst1
    OldthetaVCMEst <- thetaVCMEst
    OldthetaVCMEst1 <- thetaVCMEst1
    OldObj <-df$gammaEst + df$bisEst + oldlpred
    OldPijks <- invLogit(OldObj)
    OldQijks <- 1-OldPijks
    OldlogLikelihood <- sum(OldObj*df$y+log(OldQijks)-(df$vb + df$vgamma + 2 * df$rbgamma)*OldPijks*OldQijks/2)
    for(s in 1:10){
      #step size
      ss <- (1 / 2) ^ (s - 1)
      #same step size is used for all time points
      thetaVCMEst <- OldthetaVCMEst + ss*t(Phi_OneStep[,1:ntheta]) # Update theta(t)
      betaVCMEst <- OldbetaVCMEst + ss*t(Phi_OneStep[,(ntheta+1):(ntheta+nbeta)]) # Update beta(t)
      thetaVCMEst1 <- OldthetaVCMEst1 + ss*t(Phi_OneStep[,(ntheta+nbeta+1):(2*ntheta+nbeta)])
      betaVCMEst1 <- OldbetaVCMEst1 + ss*t(Phi_OneStep[,(2*ntheta+nbeta+1):(2*ntheta+2*nbeta)])
      # Calculate global log-likelihood with updated of theta(t) and beta(t)
      newlpred <- rowSums(t(thetaVCMEst)[df$r,]*df[,5:(4+ntheta)]) + rowSums(t(betaVCMEst)[df$r,]*df[,(5+ntheta):(4+ntheta+nbeta)])
      newObj <-df$gammaEst + df$bisEst + newlpred
      newPijks <- invLogit(newObj)
      newQijks <- 1-newPijks
      newlogLikelihood <- sum(newObj*df$y+log(newQijks)-(df$vb + df$vgamma + 2 * df$rbgamma)*newPijks*newQijks/2)
      if(newlogLikelihood > OldlogLikelihood){
        break
      }
    }
    # End line search for theta(t) and beta(t)
    Obj <- df$gammaEst + df$bisEst + newlpred
    newPijks <- invLogit(Obj) # P0ijk from the current iteration
    Diff <- max(newPijks - oldPijks) 
    if(numIter == maxNRIter) { stop("WARNING: NO CONVERGENCE") } 
    print(paste("Diff=",Diff))
    numIter <- numIter + 1   
  }
  
  # Compute standard error
  score.sigmab <- NA
  score.sigmagamma <- NA
  score.coeff <- c()
  tscore.coeff <- matrix(0, nrow = numF, ncol = 2*(ntheta + nbeta))
  std.score.local <- matrix(0, nrow = nbeta + ntheta, ncol = ngrid)
  #Standard error
  for(j in 1:ngrid) {
    to <- gridPoints[j]
    indices1 <- (abs(df$t-to) < bwThetaBeta)
    dfTemp <- df[indices1,] 
    # Score from the jth time point
    tDiff <- dfTemp$t - to
    lpredTemp <- colSums(thetaVCMEst[,j]*t(dfTemp[,5:(4+ntheta)])) + colSums(betaVCMEst[,j]*t(dfTemp[,(5+ntheta):(4+ntheta+nbeta)])) + 
      colSums(thetaVCMEst1[,j]*t(dfTemp[,5:(4+ntheta)]*tDiff)) + colSums(betaVCMEst1[,j]*t(dfTemp[,(5+ntheta):(4+ntheta+nbeta)]*tDiff))
    obj <- dfTemp$gammaEst + dfTemp$bisEst + lpredTemp
    pijks <- invLogit(obj)
    qijks <- 1 - pijks
    epan <- tDiff/bwThetaBeta
    kernels <- 0.75*(1-epan^2)
    kernels <- (kernels > 0)*kernels/.75
    a1 <- (dfTemp$vb + dfTemp$vgamma + 2 * dfTemp$rbgamma)/2*(pijks*qijks^2 - pijks^2*qijks)
    a2 <- pijks*qijks + (dfTemp$vb + dfTemp$vgamma + 2 * dfTemp$rbgamma)/2*(pijks*qijks^3 - 4*pijks^2*qijks^2 + pijks^3*qijks)
    # Construct Gradient (Score)
    zxMatrix <- data.matrix(cbind(dfTemp[,5:(4+ntheta+nbeta)],dfTemp[,5:(4+ntheta+nbeta)]*tDiff))
    # Construct Fisher information matrix
    for(fi in 1:numF){
      zxMatrix.i <- zxMatrix[dfTemp$fid==fi,]
      if(is.null(dim(zxMatrix.i))){
        zxVector.i <- t(zxMatrix.i) * ((dfTemp$y - pijks-a1)*kernels)[dfTemp$fid==fi]
      }else{
        zxVector.i <- t(zxMatrix.i) %*% ((dfTemp$y - pijks-a1)*kernels)[dfTemp$fid==fi]
      }
      tscore.coeff[fi,] <- zxVector.i
    }
    std.score.local[,j] <-   sqrt(diag(solve(t(tscore.coeff) %*% tscore.coeff)))[1:(ntheta + nbeta)]
  }
  
  # Standard errors for varying coefficient functions
  stdtheta <- std.score.local[1:ntheta,] 
  stdbeta <- std.score.local[(ntheta+1):(nbeta + ntheta), ]

  # Score for variance components
  score.sigmagamma <- .5 * ((meangamma.C^2 + vargamma.C)/sigmagammaEst^2 - 1/sigmagammaEst)
  score.sigmab <- aggregate(((meanb.C^2 + varb.C)/sigma2bEst^2 - 1/sigma2bEst), by = list(Ni), FUN = sum)[,-1]
  score.sigmab <- .5 * score.sigmab
  score.sig <- cbind(score.sigmab,score.sigmagamma)
  inf.sig <- t(score.sig) %*% score.sig
  
  # Standard errors for variance components
  sigmagammaSD <- sqrt(solve(inf.sig)[2,2])
  sigma2bSD <- sqrt(solve(inf.sig)[1,1])
  
  # Store the estimates of multilevel random effects
  bijEst <- meanb.C
  gammaiEst <- meangamma.C
  
  # Output estimation results
  MMEVCMEst[[1]] <- thetaVCMEst
  MMEVCMEst[[2]] <- betaVCMEst
  MMEVCMEst[[3]] <- sigma2bEst
  MMEVCMEst[[4]] <- sigmagammaEst
  MMEVCMEst[[5]] <- stdtheta
  MMEVCMEst[[6]] <- stdbeta
  MMEVCMEst[[7]] <- sigma2bSD
  MMEVCMEst[[8]] <- sigmagammaSD
  MMEVCMEst[[9]] <- bijEst 
  MMEVCMEst[[10]] <- gammaiEst
  MMEVCMEst[[11]] <- gridPoints
  names(MMEVCMEst) <- c("theta", "beta", "sigma2b", "sigma2gamma","thetaSE", "betaSE", "sigma2bSE", "sigma2gammaSE", 
                       "bijEst", "gammaiEst", "grid")
  return(MMEVCMEst)
}
 