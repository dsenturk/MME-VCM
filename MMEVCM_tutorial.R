## MMEVCM_tutorial.R
#############################################################################
## Description: A step-by-step implementation of MME-VCM and the associated  
## procedures described in "A Multilevel Mixed Effects Varying Coeffcient Model 
## with Multilevel Predictors and Random Effects for Modeling Hospitalization 
## Risk in Patients on Dialysis".  
#############################################################################
## Functions implemented: 
## MMEVCM_simulation.R, MMEVCM_estimation.R
#############################################################################
## Tutorial Outline:
## 1. Simulate hospitalization outcome data (MMEVCM_simulation.R)
## 2. Perform MMEVCM estimation and inference on multilevel risk factors (MMEVCM_estimation.R)
## 3. Visualization of MMEVCM results
#############################################################################

# Install missing packages
list.of.packages <- c("MASS", "statmod", "mvtnorm","bisoreg","lme4")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages) 

# Load packages  
library(MASS)
library(statmod)
library(mvtnorm)
library(bisoreg) 
library(lme4)

#############################################################################
# 1. Simulate hospitalization outcome data
#############################################################################

# NOTE: Generating one dataset with 100 facilities will take approximately eight minutes.

# Simulate one dataset from the simulation design described in Section 4.1 with 100 facilities
data <- MMEVCM_simulation(numF = 100)  # MMEVCM_simulation.R

#############################################################################
# 2. Perform MMEVCM estimation and inference on multilevel risk factors 
#############################################################################

# NOTE: Performing MMEVCM estimation with 100 facilities will take approximately eight minutes.

MMEVCMEst <- MMEVCM_estimation(data = data, bwThetaBeta = .22)  # MMEVCM_estimation.R

 #############################################################################
# 3. Visualization of MMEVCM results
############################################################################# 

gridPoints <- seq(0,1,1/19) # Grid points for theta(t) and beta(t) 

# True varying coefficient functions
# Effects of subject-level risk factors
beta1Fxn <- function(x) {
  return(cos(pi*x))
}
beta2Fxn <- function(x) {
  return(-cos(pi*x))
}

# Effects of facility-level risk factors
theta1Fxn <- function(x) {
  return(sin(pi*x))
}
theta2Fxn <- function(x) {
  return(-sin(pi*x))
}

# The intercept term
beta0Fxn <- function(x) {
  return(sqrt(x) - 1.8)
}

invLogit <- function(x) {
  return(exp(x)/(1+exp(x)))
}

theta1F <- theta1Fxn(gridPoints)
theta2F <- theta2Fxn(gridPoints)
beta1F <- beta1Fxn(gridPoints)
beta2F <- beta2Fxn(gridPoints) 
beta0F <- beta0Fxn(gridPoints)
ntheta <- 2
nbeta <- 3

thetaVCMEst <- MMEVCMEst$theta
betaVCMEst <- MMEVCMEst$beta
sigma2bEst <- MMEVCMEst$sigma2b
sigma2gammaEst <- MMEVCMEst$sigma2gamma
stdbeta <- MMEVCMEst$betaSE
stdtheta <- MMEVCMEst$thetaSE
sigma2bSE <- MMEVCMEst$sigma2bSE
sigma2gammaSE <- MMEVCMEst$sigma2gammaSE

# Plot estimates of varying coefficient functions and +-2 standard errors
par(mfrow=c(3,2))

# Plot beta0(t)
beta0U <- betaVCMEst[3,] + 2 * stdbeta[3,]
beta0L <- betaVCMEst[3,] - 2 * stdbeta[3,]
u <- max(beta0U)
l <- min(beta0L)
plot(gridPoints,betaVCMEst[3,],'l',col="black",lwd=2,main="(d)", ylim = c(l,u),xlim=c(0,1),xaxs = "i",cex.main=1.8,xlab = "",ylab = "",cex.axis=1.4) 
title(xlab = "t",ylab=expression(widehat(beta)[0](t)), line=2, cex.lab=1.5)
lines(gridPoints,beta0F,lwd=2,lty=2)
polygon(c(gridPoints, rev(gridPoints)), c(beta0L, rev(beta0U)), col = rgb(.25, .25, .25, 0.5), border = NA)

# Plot beta1(t)
beta1U <- betaVCMEst[1,] + 2 * stdbeta[1,]
beta1L <- betaVCMEst[1,] - 2 * stdbeta[1,]
u <- max(beta1U)
l <- min(beta1L)
plot(gridPoints,betaVCMEst[1,],'l',col="black",lwd=2,main="(c)", ylim = c(l,u),xlim=c(0,1),xaxs = "i",cex.main=1.8,xlab = "",ylab = "",cex.axis=1.4) 
title(xlab = "t",ylab=expression(widehat(beta)[1](t)), line=2, cex.lab=1.5)
lines(gridPoints,beta1F,lwd=2,lty=2)
polygon(c(gridPoints, rev(gridPoints)), c(beta1L, rev(beta1U)), col = rgb(.25, .25, .25, 0.5), border = NA)

# Plot beta2(t)
beta2U <- betaVCMEst[2,] + 2 * stdbeta[2,]
beta2L <- betaVCMEst[2,] - 2 * stdbeta[2,]
u <- max(beta2U)
l <- min(beta2L)
plot(gridPoints,betaVCMEst[2,],'l',col="black",lwd=2,main="(d)", ylim = c(l,u),xlim=c(0,1),xaxs = "i",cex.main=1.8,xlab = "",ylab = "",cex.axis=1.4) 
title(xlab = "t",ylab=expression(widehat(beta)[2](t)), line=2, cex.lab=1.5)
lines(gridPoints,beta2F,lwd=2,lty=2)
polygon(c(gridPoints, rev(gridPoints)), c(beta2L, rev(beta2U)), col = rgb(.25, .25, .25, 0.5), border = NA)



# Plot theta1(t)
theta1U <- thetaVCMEst[1,] + 2 * stdtheta[1,]
theta1L <- thetaVCMEst[1,] - 2 * stdtheta[1,]
u <- max(theta1U)
l <- min(theta1L)
plot(gridPoints,thetaVCMEst[1,],'l',col="black",lwd=2,main="(a)", ylim = c(l,u),xlim=c(0,1),xaxs = "i",cex.main=1.8,xlab = "",ylab = "",cex.axis=1.4) 
title(xlab = "t",ylab=expression(widehat(theta)[1](t)), line=2, cex.lab=1.5)
lines(gridPoints,theta1F,lwd=2,lty=2)
polygon(c(gridPoints, rev(gridPoints)), c(theta1L, rev(theta1U)), col = rgb(.25, .25, .25, 0.5), border = NA)

# Plot theta2(t)
theta2U <- thetaVCMEst[2,] + 2 * stdtheta[2,]
theta2L <- thetaVCMEst[2,] - 2 * stdtheta[2,]
u <- max(theta2U)
l <- min(theta2L)
plot(gridPoints,thetaVCMEst[2,],'l',col="black",lwd=2,main="(b)", ylim = c(l,u),xlim=c(0,1),xaxs = "i",cex.main=1.8,xlab = "",ylab = "",cex.axis=1.4) 
title(xlab = "t",ylab=expression(widehat(theta)[2](t)), line=2, cex.lab=1.5)
lines(gridPoints,theta2F,lwd=2,lty=2)
polygon(c(gridPoints, rev(gridPoints)), c(theta2L, rev(theta2U)), col = rgb(.25, .25, .25, 0.5), border = NA)


# Print the estimated variance components
cat(paste("Estimated variance of subject-specific random effects:", sigma2bEst),"\n", paste("Standard error:", sigma2bSE), "\n", paste("True variance of subject-specific random effects:", 1),
    "\n", paste("Estimated variance of facility-specific random effects:", sigma2gammaEst),"\n", paste("Standard error:", sigma2gammaSE), "\n", 
    paste("True variance of facility-specific random effects:", 1))

dev.off()