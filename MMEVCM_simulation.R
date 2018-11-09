MMEVCM_simulation <- function(numF # number of facilities (scalar)
                             ){
  
  ##############################################################################
  ## Description: Function for simulating one data set under the simulation design described
  ##              in Section 4.1.
  ## Args: see above
  ## Returns: data, data.frame with columns c("fid", "sid", "y", "t", "z1", "z2", "x1", "x2", "x0")
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
  ############################################################################## 
  
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
  
  invLogit <- function(x) {
    return(exp(x)/(1+exp(x)))
  }
  
  Selection <- function(x,ngrid) {
    OBpoints <- 1:ngrid
    return(OBpoints[OBpoints<=x])
  }
  
  # Function of subject-level effects
  beta1Fxn <- function(x) {
    return(cos(pi*x))
  }
  
  beta2Fxn <- function(x) {
    return(-cos(pi*x))
  }
  
  # Function of facility-level effects
  theta1Fxn <- function(x) {
    return(sin(pi*x))
  }
  
  theta2Fxn <- function(x) {
    return(-sin(pi*x))
  }
  
  # Function of the intercept term
  beta0Fxn <- function(x) {
    return(sqrt(x) - 1.8)
  }
  
  # Functions for facility-level covariates Zi(j)
  z1Fxn <- function(x){
    return(0.1*x)
  }
  
  z2Fxn <- function(x){
    return(-0.1*x)
  }
  
  # Variance of subject-specific random effects 
  sigma2b <- 1
  # Variance of facility-specific random effects
  sigmagamma <- 1
  
  # Number of subject-level risk factors (including the intercept term)
  nbeta <- 3
  
  # Number of facility-level risk factors
  ntheta <- 2
  
  # Define the grid points used for the varying coefficient functions theta(t) and beta(t)
  gridPoints <- seq(0,1,1/19)
  ngrid <- length(gridPoints)
  
  # Grid points of initiation time of dialysis
  etaGrid <- seq(0,1,1/19)
  nEtagrid <- length(etaGrid)
   
  
  #######################################################
  # Construct Data set
  #######################################################
  
  # Generate number of subjects per facility from a half normal distribution
  xx <- rnorm(numF * 4,sd=30)
  x <- abs(xx)[abs(xx)>=20]
  x <- round(x) 
  numSubPF <- x[1:numF]
  sumNi <- sum(numSubPF) # Total number of subjects
  Ni <- rep(1:numF,numSubPF)
  numOB <- rep(20, sumNi) # Number of observations per subject
  sumNij <- sum(numOB)
  
  # Create data.frame to store dataset
  df <- data.frame(matrix(ncol = 9, nrow = sumNij))
  colnames(df) <- c("fid", "sid", "y", "t", "z1", "z2", "x1", "x2", "x0")
  df$fid <- rep(rep(1:numF,numSubPF),numOB)
  numDisPF <- as.numeric(table(df$fid))
  df$sid <- rep(1:sumNi,numOB)
  df$t <- gridPoints[unlist(apply(matrix(numOB),1,Selection,ngrid))]
  IniDiaTime <- runif(sumNi)
  rIniDiaTime <- round(19*IniDiaTime) + 1
  df$c <- rep(IniDiaTime,numOB) 
  df$r <- round(19*df$t)+1 
  df$cr <- rep(rIniDiaTime,numOB) 
  
  # Generate multilevel random effects
  df$bi <- rep(rnorm(sumNi, 0, sqrt(sigma2b)), numOB) 
  df$gamma <- rep(rnorm(numF, 0, sqrt(sigmagamma)), numDisPF)
  
  # Generate subject-level covariates Xij
  covariates1 <- mvrnorm(sumNi, c(0,0), matrix(c(.125,.0625,.0625,.125),2,2))
  x1 <- covariates1[,1]
  x2 <- covariates1[,2]
  df$x1 <- c(rep(covariates1[,1],each=ngrid))
  df$x2 <- c(rep(covariates1[,2],each=ngrid))
  df$x0 <- 1 # Adding the intercept
  
  # Generate facility-level covariates Zi(j)
  eps <- mvrnorm(numF*nEtagrid, c(0,0), matrix(c(.125,.0625,.0625,.125),2,2))
  eps1 <- eps[(Ni-1)*nEtagrid+rIniDiaTime, 1]
  eps2 <- eps[(Ni-1)*nEtagrid+rIniDiaTime, 2]
  df$z1 <- z1Fxn(df$c)+rep(eps1,numOB) # continuous facility-level covariates Zi(j)*
  df$z2 <- z2Fxn(df$c)+rep(eps2,numOB)
  DIniDiaTime <- floor((rIniDiaTime)/7) * 7
  DIniDiaTime[DIniDiaTime==0] <- 1
  Dc <- floor((df$cr)/7) * 7 
  Dc[Dc==0] <- 1
  Deps1 <- eps[(Ni-1)*nEtagrid+DIniDiaTime, 1]
  Deps2 <- eps[(Ni-1)*nEtagrid+DIniDiaTime, 2]
  df$Dz1 <- z1Fxn((Dc-1)/19)+rep(Deps1,numOB)
  df$Dz2 <- z2Fxn((Dc-1)/19)+rep(Deps2,numOB)
  df$z1 <- df$Dz1 # discrete facility-level covariates Zi(j)
  df$z2 <- df$Dz2
  
  # Generate longitudinal effects of multilevel risk factors
  df$beta1 <- beta1Fxn(df$t)
  df$beta2 <- beta2Fxn(df$t) 
  df$beta0 <- beta0Fxn(df$t)
  df$theta1 <- theta1Fxn(df$t)
  df$theta2 <- theta2Fxn(df$t)
  
  # Generate Yijk* and Sij jointly as described in Section 4.1
  eta_new <- df$gamma + df$bi + df$beta1*df$x1 + df$beta2*df$x2 + df$beta0*df$x0 + df$theta1*df$z1 + df$theta2*df$z2
  prob_new <- invLogit(eta_new)
  survival_mu <- 0.8 
  survival_var <- 0.04
  covSY <- -0.01
  covYY <- 0
  p_alive <- pnorm(gridPoints, mean = survival_mu, sd = sqrt(survival_var), lower.tail = FALSE)
  joint_cov <- matrix(c(1,covSY,covSY,survival_var), nrow = 2)
  J_cov <- diag(c(rep(1,ngrid),survival_var), ngrid+1)
  J_cov[upper.tri(J_cov)] <- covYY
  J_cov[lower.tri(J_cov)] <- covYY
  J_cov[,ngrid+1] <- covSY
  J_cov[ngrid+1,] <- covSY
  J_cov[ngrid+1,ngrid+1] <- survival_var
  Uncond_mean <- rep(0,ngrid)
  y_pc <- c()
  df_index <- c()
  for(i in 1:sumNi){
    joint_p <- prob_new[df$sid==i] * p_alive
    ####start bisection####
    for(j in 1:ngrid){
      p <- joint_p[j]
      d <- 1
      lowerB <- -2
      upperB <- 2
      while(d > 0.0001){
        u <- 0.5 * (lowerB + upperB)
        p1 <- pmvnorm(lower = c(0,gridPoints[j]), upper = Inf, mean = c(u,survival_mu), sigma = joint_cov)
        if(p1[1] <= p){
          lowerB <- u
        } else{
          upperB <- u
        }
        d <- upperB - lowerB
      }
      Uncond_mean[j] <- u
    }
    Y <- rmvnorm(1, mean = c(Uncond_mean,survival_mu), sigma = J_cov)
    numOB[i] <- min(floor(Y[ngrid + 1] * 20), 20)
    if (numOB[i]<=0) numOB[i] <- 1
    df_index <- c(df_index, (1:numOB[i]) + 20*(i-1))
    y_pc <- c(y_pc, as.numeric((Y > 0)[1 : numOB[i]]))
  }
  
  df <- df[df_index,] # Truncation of outcome
  df$y <- y_pc # Patients' hospitalization outcome Yijk
  
  return(df[,1:9])
}
  