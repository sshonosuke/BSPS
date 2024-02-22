###----------------------------------------------### 
###        Oneshot simulation study 1            ###
###----------------------------------------------###
rm(list=ls())
set.seed(1)

## load R codes and packages 
source("BSPS.R")
library(spgwr)   
library(MASS)
library(MCMCpack)
library(RandomForestsGLS)
library(gam)
library(spNNGP)

## settings 
n_train <- 300   # the number of locations (including non-sampled locations)
n_test <- 200   # the number of non-sampled locations 
n <- n_train + n_test


## generation of sampling locations
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


## data generation 
M1 <- x1 - 0.5*x2^2 
M2 <- x1^2 + x2^2
ID1 <- ifelse(coords[,1]<0, 1, 0)
ID2 <- ifelse(coords[,1]>0, 1, 0)
sp_cov <- exp(-as.matrix(dist(coords))/0.3)
w <- as.vector(mvrnorm(1, rep(0,n), sp_cov))
Mu <- ID1*M1 + ID2*M2 + 0.3*w 
Sig <- 1
y <- rnorm(n, Mu, Sig)
ID <- cbind(ID1, ID2)




## training & test data 
coords_train <- coords[1:n_train,]
y_train <- y[1:n_train]
x_train <- as.matrix(x[1:n_train,])
ID_train <- ID[1:n_train,]

coords_test <- coords[(n_train+1):n,]
y_test <- y[(n_train+1):n]
x_test <- as.matrix(x[(n_train+1):n,])
ID_test <- ID[(n_train+1):n,]



## Model 1: polynomial regression (left region)
EX_train <- cbind(x_train, x_train^2)
sub <- (ID_train[,1]==1)
fit <- lm(y_train[sub]~EX_train[sub,])
a1 <- as.vector( cbind(1,EX_train)%*%coef(fit) )   # point prediction (used in BPS) 
b1 <- rep(summary(fit)$sigma^2, n_train)   # prediction variance (used in BPS)

EX_test <- cbind(x_test, x_test^2)
fa1 <- as.vector( cbind(1,EX_test)%*%coef(fit) )
fb1 <- rep(summary(fit)$sigma^2, n_test)


## Model 2: polynomial regression (right region)
EX_train <- cbind(x_train, x_train^2)
sub <- (ID_train[,2]==1)
fit <- lm(y_train[sub]~EX_train[sub,])
a2 <- as.vector( cbind(1,EX_train)%*%coef(fit) )   # point prediction (used in BPS) 
b2 <- rep(summary(fit)$sigma^2, n_train)   # prediction variance (used in BPS)

EX_test <- cbind(x_test, x_test^2)
fa2 <- as.vector( cbind(1,EX_test)%*%coef(fit) )
fb2 <- rep(summary(fit)$sigma^2, n_test)


## Model 3: quadratic regression with spatial effects
EX_train <- cbind(x_train, x_train^2)
priors <- list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2,5), "tau.sq.IG"=c(2,5))
starting <- list("phi"=6, "sigma.sq"=5, "tau.sq"=5)
tuning <- list("phi"=0.5, "sigma.sq"=0.5, "tau.sq"=0.5)
NNGP_train <- spNNGP(y_train~EX_train, coords=coords_train, 
                     method="latent", n.neighbors=10, cov.model="exponential",
                     priors=priors, starting=starting, tuning=tuning, n.samples=2000)
pred.NNGP <- cbind(1,EX_train)%*%t(NNGP_train$p.beta.samples) + NNGP_train$p.w.samples + t(NNGP_train$p.theta.samples[,2]*matrix(rnorm(2000*n_train), 2000, n_train))
a3 <- apply(pred.NNGP, 1, mean)
b3 <- apply(pred.NNGP, 1, var)

EX_test <- cbind(x_test, x_test^2)
pred.NNGP <- predict(NNGP_train, X.0=cbind(1,EX_test), coords.0=coords_test)
fa3 <- apply(pred.NNGP$p.y.0, 1, mean)
fb3 <- apply(pred.NNGP$p.y.0, 1, var)


## Model 4: additive model 
gfit <- gam(y_train~s(x1)+s(x2)+s(x3)+s(x4)+s(x5), data=data.frame(x_train))
a4 <- predict(gfit)
b4 <- rep(sum(gfit$residuals^2)/gfit$df.residual, n_train)
fa4 <- predict(gfit, data.frame(x_test))
fb4 <- rep(sum(gfit$residuals^2)/gfit$df.residual, n_test)


## BSPS (setup)
m <- 10
mc <- 2000
bn <- 1000
bandwidth_grid.set <- seq(0.2, 1, by=0.05)[-1]

# BSPS1 (QR1+QR2)
mean_models <- cbind(a1, a2)
var_models <- cbind(b1, b2)
mean_models_test <- cbind(fa1, fa2)
var_models_test <- cbind(fb1, fb2)
fit1 <- BSPS_nngp(Y=y_train, SP=coords_train, mF=mean_models, vF=var_models, m=m, 
                 mcmc=mc, burn=bn, Phi.set=bandwidth_grid.set, print=T,
                 prediction=T, SP_test=coords_test, mF_test=mean_models_test, 
                 vF_test=var_models_test)
pred_BSPS1 <- apply(fit1$pred, 2, mean)    # point prediction via BPS
PI1 <- apply(fit1$pred, 2, quantile, prob=c(0.025, 0.975))    # prediction interval

# BSPS2 (QR1+QR2+SPR)
mean_models <- cbind(a1, a2, a3)
var_models <- cbind(b1, b2, b3)
mean_models_test <- cbind(fa1, fa2, fa3)
var_models_test <- cbind(fb1, fb2, fb3)
fit2 <- BSPS_nngp(Y=y_train, SP=coords_train, mF=mean_models, vF=var_models, m=m, 
                  mcmc=mc, burn=bn, Phi.set=bandwidth_grid.set, print=T,
                  prediction=T, SP_test=coords_test, mF_test=mean_models_test, 
                  vF_test=var_models_test)
pred_BSPS2 <- apply(fit2$pred, 2, mean)    # point prediction via BPS
PI2 <- apply(fit2$pred, 2, quantile, prob=c(0.025, 0.975))    # prediction interval

# BSPS3 (AM+SPR)
mean_models <- cbind(a3, a4)
var_models <- cbind(b3, b4)
mean_models_test <- cbind(fa3, fa4)
var_models_test <- cbind(fb3, fb4)
fit3 <- BSPS_nngp(Y=y_train, SP=coords_train, mF=mean_models, vF=var_models, m=m, 
                  mcmc=mc, burn=bn, Phi.set=bandwidth_grid.set, print=T,
                  prediction=T, SP_test=coords_test, mF_test=mean_models_test, 
                  vF_test=var_models_test)
pred_BSPS3 <- apply(fit3$pred, 2, mean)    # point prediction via BPS
PI3 <- apply(fit3$pred, 2, quantile, prob=c(0.025, 0.975))    # prediction interval


## simple average (SA)
pred_SA1 <- 0.5*fa1 + 0.5*fa2
pred_SA2 <- fa1/3 + fa2/3 + fa3/3
pred_SA3 <- 0.5*fa3 + 0.5*fa4


## Bayesian Model Averaging
BIC1 <- -2*sum(dnorm(y_train, a1, sqrt(b1), log=T)) + log(n_train)*(dim(x_train)[2]+2)
BIC2 <- -2*sum(dnorm(y_train, a2, sqrt(b2), log=T)) + log(n_train)*(dim(x_train)[2]+2)
w1 <- (-0.5)*BIC1
w2 <- (-0.5)*BIC2
w3 <- (-0.5)*spDiag(NNGP_train)$WAIC[1,]
w4 <- (-0.5)*AIC(gfit)

max_w <- max(w1, w2)
ww1 <- exp(c(w1, w2)-max_w)/(exp(w1-max_w)+exp(w2-max_w))
pred_BMA1 <- ww1[1]*fa1 + ww1[2]*fa2

max_w <- max(w1, w2, w3)
ww2 <- exp(c(w1, w2, w3)-max_w)/(exp(w1-max_w)+exp(w2-max_w)+exp(w3-max_w))
pred_BMA2 <- ww2[1]*fa1 + ww2[2]*fa2 + ww2[3]*fa3

max_w <- max(w3, w4)
ww3 <- exp(c(w3, w4)-max_w)/(exp(w3-max_w)+exp(w4-max_w))
pred_BMA3 <- ww3[1]*fa3 + ww3[2]*fa4


## evaluation (MSE) 
mean( (pred_BSPS1 - y_test)^2 )
mean( (pred_SA1 - y_test)^2 )
mean( (pred_BMA1 - y_test)^2 )

mean( (pred_BSPS2 - y_test)^2 )
mean( (pred_SA2 - y_test)^2 )
mean( (pred_BMA2 - y_test)^2 )

mean( (pred_BSPS3 - y_test)^2 )
mean( (pred_SA3 - y_test)^2 )
mean( (pred_BMA3 - y_test)^2 )

mean( (fa1 - y_test)^2 )
mean( (fa2 - y_test)^2 )
mean( (fa3 - y_test)^2 )
mean( (fa4 - y_test)^2 )


# coverage prob
mean(PI1[1,]<y_test & y_test<PI1[2,])
mean(PI2[1,]<y_test & y_test<PI2[2,])
mean(PI3[1,]<y_test & y_test<PI3[2,])



##  spatial plot function
SPlot <- function(Sp, value, ran, xlim=NULL, title=""){
  cs <- colorRamp( c("blue", "green", "yellow", "red"), space="rgb")
  value <- (value-ran[1])/(ran[2]-ran[1])
  cols <- rgb( cs(value), maxColorValue=256 )
  plot(Sp, col=cols, ylab="Latitude", xlab="Longitude", xlim=c(-1, 1.5), main=title, pch=20, cex=2)
  cs <- colorRamp( c("blue", "green", "yellow", "red"), space="rgb")
  cols <- rgb(cs(0:1000/1000), maxColorValue=256)
  rect(1.25, seq(-1,1,length=1001), 1.45, seq(-1,1,length=1001), col=cols, border=cols)
  tx <- round(seq(ran[1], ran[2], length=5), 2)
  text(x=1.5, y=seq(-1, 1, length=5), tx, cex=0.7)
  yy <- seq(0, 2, length=5)
  for (i in 1:5){
    segments(1.1, yy[i], 1.3, yy[i], col="white")
  }
}


##  plot (weight & interval prediction)
hBeta <- apply(fit1$beta, c(2,3), mean)
head(hBeta)
ab.hBeta <- abs(hBeta)
sel <- ab.hBeta[,2]/(ab.hBeta[,2]+ab.hBeta[,3])


par(mfcol=c(1,2))
# weight
SPlot(coords_train, sel, ran=c(0,1))
abline(v=0)
# interval
ran <- range(pred_BSPS1)
plot(y_test, pred_BSPS1, xlim=ran, ylim=ran, pch=20, cex=1.5, 
     xlab="true value", ylab="prediction")
abline(0, 1)
for(i in 1:n_test){
  lines(rep(y_test[i],2), PI1[,i], col=1)
}
points(y_test, pred_BMA1, col=2, pch=2)
points(y_test, pred_SA1, col=4, pch=8)
legend("bottomright", legend=c("point prediction (BMA)", "point prediction (SA)", "point prediction (BSPS)", "95% prediction interval"), 
       pch=c(2, 8, 20, NA), lty=c(NA, NA, NA, 1), col=c(2,4,1,1))
