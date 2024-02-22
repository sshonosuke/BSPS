###----------------------------------------------### 
###        Oneshot simulation study 2            ###
###----------------------------------------------###
rm(list=ls())
set.seed(2)

## load R codes and packages 
source("BSPS.R")
library(spgwr)   
library(MASS)
library(MCMCpack)
library(RandomForestsGLS)
library(gam)
library(spNNGP)


## settings
scenario <- 2   # DGP scenario (1 or 2)
p <- 5   # number of covariates (5 or 15)
add_p <- p - 5


## definition of function 
CI <- function(x){ quantile(x, prob=c(0.025, 0.975)) }


## sample size
n_train <- 300   # the number of locations (including non-sampled locations)
n_test <- 100   # the number of non-sampled locations 
n <- n_train + n_test


## spatial location
coords <- cbind(runif(n,-1,1), runif(n,-1,1))

##  covariates
phi <- 0.5    # range parameter for covariates
dd <- as.matrix(dist(coords))
mat <- exp(-dd/phi)
z1 <- mvrnorm(1, rep(0, n), mat)
z2 <- mvrnorm(1, rep(0, n), mat)
x1 <- z1
rr <- 0.2
x2 <- rr*z1 + sqrt(1-rr^2)*z2
x3 <- rnorm(n)
x4 <- rnorm(n)
x5 <- rnorm(n)
x <- cbind(x1, x2, x3, x4, x5)
if(add_p>0){
  x_add <- matrix(rnorm(n*add_p), n, add_p)
  x <- cbind(x, x_add)
}

## data generation 
sp_cov <- exp(-as.matrix(dist(coords))/0.3)
if(scenario==1){
  w <- as.vector(mvrnorm(1, rep(0,n), sp_cov))
  Mu <- w + x3^2*exp(-0.3*(coords[,1]^2+coords[,2]^2)) + sin(2*x2)*coords[,2]
}
if(scenario==2){
  w1 <- as.vector(mvrnorm(1, rep(0,n), sp_cov))
  w2 <- as.vector(mvrnorm(1, rep(0,n), sp_cov))
  w3 <- as.vector(mvrnorm(1, rep(0,n), sp_cov))
  Mu <- 2*w1 + (10*sin(pi*x1*x2) + 20*(x3-0.5)^2 + 10*x4 + 5*x5)/20
}
Sig <- 0.7
y <- rnorm(n, Mu, Sig)

## split data into train and test data
coords_train <- coords[1:n_train,]
y_train <- y[1:n_train]
x_train <- as.matrix(x[1:n_train,])

coords_test <- coords[(n_train+1):n,]
y_test <- y[(n_train+1):n]
x_test <- as.matrix(x[(n_train+1):n,])


### Three base models
## GWR
bandwidth.opt <- gwr.sel(y_train~x_train, coords=coords_train, verbose=F)
fit_gwr <- gwr(y_train~x_train, coords=coords_train, 
               bandwidth=bandwidth.opt, se.fit=T, hatmatrix=T, predictions=T)

mean_train_gwr <- fit_gwr$SDF$pred    # point prediction (used in BPS) 
sigma2 <- fit_gwr$results$sigma2.b
var_train_gwr <- rep(sigma2, n_train)    # prediction variance (used in BPS)

distance_test <- as.matrix(dist(coords))[(n_train+1):n,]
distance_test <- distance_test[,1:n_train]
grw_gauss <- gwr.Gauss(distance_test^2, bandwidth.opt)
x_train_const <- cbind(1, x_train)
gw_beta <- matrix(NA, n_test, ncol(x_train_const))
var_test_gwr <- c()
for(i in 1:n_test){
  mat <- solve(t(x_train_const)%*%(x_train_const*grw_gauss[i,]))
  gw_beta[i,] <- mat%*%t(x_train_const)%*%(y_train*grw_gauss[i,])
  var_test_gwr[i] <- sigma2
}
mean_test_gwr <- apply(cbind(1, x_test)*gw_beta, 1, sum)    # point prediction in non-sampled locations

## Additive model 
XX_train <- data.frame(x_train, coords_train)
XX_test <- data.frame(x_test, coords_test)
names(XX_train) <- names(XX_test) <- paste0("X", 1:p)
formula_gam <- as.formula( paste0("y_train~", paste0("s(X", 1:p, ")", collapse="+")) )
fit_gam <- gam(formula_gam, data=XX_train)
mean_train_gam <- predict(fit_gam)   # point prediction (used in BPS) 
var_train_gam <- summary(fit_gam)$dispersion   # prediction variance (used in BPS)
mean_test_gam <- predict(fit_gam, XX_test)  # point prediction in non-sampled locations
var_test_gam <- rep(summary(fit_gam)$dispersion, n_test)   # prediction variance (used in BPS)

## Spatial regression
priors <- list("phi.Unif"=c(0.01, 0.5), "sigma.sq.IG"=c(1,1), "tau.sq.IG"=c(1,1))
starting <- list("phi"=0.1, "sigma.sq"=1, "tau.sq"=1)
tuning <- list("phi"=0.1, "sigma.sq"=0.3, "tau.sq"=0.3)
NNGP_train <- spNNGP(y_train~x_train, coords=coords_train, method="latent", n.neighbors=10, 
                     cov.model="exponential", priors=priors, starting=starting, tuning=tuning, 
                     n.samples=2000, verbose=F)
XX_test <- as.matrix( cbind(1, x_test) )
class(XX_test) <- "numeric"
pred.NNGP <- predict(NNGP_train, X.0=XX_test, coords.0=coords_test)
mean_test_nngp <- apply(pred.NNGP$p.y.0, 1, mean)
var_test_nngp <- apply(pred.NNGP$p.y.0, 1, var)
XX_train <- as.matrix( cbind(1, x_train) )
mu_NNGP <- XX_train%*%t(NNGP_train$p.beta.samples) + NNGP_train$p.w.samples + t(NNGP_train$p.theta.samples[,2]*matrix(rnorm(2000*n_train), 2000, n_train))
mean_train_nngp <- apply(mu_NNGP, 1, mean)
var_train_nngp <- apply(mu_NNGP, 1, var)



## BSPS
# input parameters
mc <- 2000
bn <- 1000
bandwidth_grid.set <- seq(0.2, 0.8, by=0.05)[-1]
mean_models <- cbind(mean_train_gwr, mean_train_gam, mean_train_nngp)
var_models <- cbind(var_train_gwr, var_train_gam, var_train_nngp)
mean_models_test <- cbind(mean_test_gwr, mean_test_gam, mean_test_nngp)
var_models_test <- cbind(var_test_gwr, var_test_gam, var_test_nngp)

PI_l <- mean_models_test-1.96*sqrt(var_models_test)
PI_u <- mean_models_test+1.96*sqrt(var_models_test)

# BSPS with NNGP
m <- 10
fit_nngp <- BSPS_nngp(Y=y_train, SP=coords_train, mF=mean_models, vF=var_models, m=m, 
                      mcmc=mc, burn=bn, Phi.set=bandwidth_grid.set, print=T,
                      prediction=T, SP_test=coords_test, mF_test=mean_models_test, 
                      vF_test=var_models_test)
pred_BSPS_nngp <- apply(fit_nngp$pred, 2, mean)    # point prediction via BSPS with NNGP
PI_BSPS_nngp <- apply(fit_nngp$pred, 2, CI)    # prediction interval via BSPS with NNGP

# BSPS with variational approximation
fit_vb <- BSPS_vb(Y=y_train, SP=coords_train, mF=mean_models, vF=var_models, 
                  mc=mc-bn, Phi.set=bandwidth_grid.set, print=F,
                  prediction=T, SP_test=coords_test, mF_test=mean_models_test, 
                  vF_test=var_models_test)
pred_BSPS_vb <- fit_vb$pred     # point prediction via BPS-VB



## SpatialRF
estimation_result <- RFGLS_estimate_spatial(coords_train, y_train, x_train, ntree=50)
prediction_result <- RFGLS_predict_spatial(estimation_result, coords_test, x_test)
sig2 <- mean((y_train-estimation_result$predicted)^2)

mean_train_sprf <- estimation_result$predicted
var_train_sprf <- rep(sig2, n_train)
mean_test_sprf <- prediction_result$prediction
var_test_sprf <- rep(sig2, n_test)


## Bayesian Model Averaging
b1 <- (-0.5)*fit_gwr$results$AICc
b2 <- (-0.5)*AIC(fit_gam)
b3 <- (-0.5)*spDiag(NNGP_train)$WAIC[1,]
bb <- max(b1, b2, b3)
ww <- exp(c(b1, b2, b3)-bb)/(exp(b1-bb)+exp(b2-bb)+exp(b3-bb))
Pred_BMA <- apply(t(mean_models_test)*ww, 2, sum)


## simple average
Pred_SA <- apply(mean_models_test, 1, mean)


## Bates-Granger average
weight_test <- (1/var_models_test) / apply(1/var_models_test, 1, sum)
Pred_BG <- apply(mean_models_test*weight_test, 1, sum)

  
 
## performance evaluation 
# point prediction 
Pred <- cbind(pred_BSPS_nngp, pred_BSPS_vb, mean_models_test, mean_test_sprf, Pred_BMA, Pred_SA, Pred_BG)
colMeans((Pred-y_test)^2)

# prediction interval (coverage)
mean(PI_BSPS_nngp[1,]<y_test & PI_BSPS_nngp[2,]>y_test)
apply(PI_l<y_test & PI_u>y_test, 2, mean)

# prediction interval (length)
mean(PI_BSPS_nngp[2,]-PI_BSPS_nngp[1,])
apply(PI_u-PI_l, 2, mean)



  
  