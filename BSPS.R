###----------------------------------------------------### 
###    Bayesian spatial predictive synthesis (BSPS)    ###
###       with nearest-neighbor Gaussian process       ###
###----------------------------------------------------###

## INPUT
# Y: response vector 
# SP: location information
# mF: matrix of point prediction of agent (without intercept)
# vF: matrix of prediction variance of agent 
# m: number of neighbors in NNGP (default is 10)
# Remark: exponential covariance function is used.  

## OUTPUT
# posterior samples of beta, sigma, tau, phi (bandwidth)

BSPS_nngp <- function(Y, SP, mF, vF, m=10, mcmc=1000, burn=200, Phi.set=NULL, print=T,
                      prediction=F, SP_test=NULL, mF_test=NULL, vF_test=NULL, constrain=F){
  # preparation
  n <- length(Y)     # the number of sampled locations
  MF <- cbind(1, mF)
  VF <- cbind(1, vF)
  J <- dim(MF)[2] 
  prior_var <- 1     # hyperparameter in inverse gamma prior for variance parameters 
  
  # correlation function 
  CF <- function(dist.mat, phi){ exp(-dist.mat/phi) }   
  
  # neighbor location 
  dd <- as.matrix( dist(SP) )
  NN <- matrix(NA, n, m)
  for(i in 2:n){
    if(i<=m){  NN[i, 1:(i-1)] <- sort(order(dd[i,1:(i-1)]))  }
    if(i>m){   NN[i,] <- sort(order(dd[i,1:(i-1)])[1:m])  }
  }
  
  # set of children indices 
  In <- function(x, i){ i %in% x }
  UU <- list()
  for(i in 1:n){
    UU[[i]] <- (1:n)[apply(NN, 1, In, i=i)]
  }
  
  # set of spatial range parameters
  L <- length(Phi.set)
  if(is.null(Phi.set)){
    L <- 10
    Phi.range <- c(0, median(dist(SP)))
    Phi.set <- seq(Phi.range[1], Phi.range[2], length=L+1)[-1]
  }
  
  # set of prior parameters in NNGP for each candidate value of Phi
  BB.set <- list()    # coefficient for neighbor locations
  FF.set <- list()    # conditional variance 
  for(l in 1:L){
    BB <- matrix(NA, n, m)
    FF <- c()
    FF[1] <- 1
    for(i in 2:n){
      if(i<=m){
        mat <- CF(as.matrix( dist(SP[c(NN[i,1:(i-1)], i),])), Phi.set[l])
        C1 <- solve(mat[1:(i-1), 1:(i-1)])
        C2 <- mat[i, 1:(i-1)]
        FF[i] <- mat[i, i] - t(C2)%*%C1%*%C2   # conditional variance  
        BB[i,1:(i-1)] <- as.vector(C2%*%C1)   # coefficient for neighbour variables 
      }
      if(i>m){
        mat <- CF(as.matrix(dist(SP[c(NN[i,], i),])), Phi.set[l])
        C1 <- solve(mat[1:m, 1:m])
        C2 <- mat[m+1, 1:m]
        FF[i] <- mat[m+1, m+1] - t(C2)%*%C1%*%C2   # conditional variance  
        BB[i,] <- as.vector(C2%*%C1)   # coefficient for neighbour variables 
      }
    }
    BB.set[[l]] <- BB
    FF.set[[l]] <- FF
  }
  
  # matrices/arrays to store posterior samples 
  Beta.pos <- array(NA, c(mcmc, n, J))
  Phi.pos <- matrix(NA, mcmc, J)
  Tau.pos <- matrix(NA, mcmc, J)
  Sig.pos <- c()
  
  # initial values
  Beta <- cbind(0, matrix(1/(J-1), n, J-1))     # coefficients 
  Tau <- rep(0.7, J)     # scale parameters in NNGP
  Sig <- 0.2     # error variance 
  Phi.index <- rep(round(L/2), J)     # spatial ranges in NNGP
  TP <- MF     # true predictor 
  Beta_center <- c(0, rep(1/(J-1), J-1))   # fixed
  
  # MCMC iteration 
  BB <- array(NA, c(n, m, J))
  FF <- matrix(NA, n, J)
  
  for(r in 1:mcmc){
    # Beta (update centered model weight)
    sY <- Y - apply(t(TP)*Beta_center, 2, sum)
    sBeta <- t( t(Beta) - Beta_center )   # centered model weight 
    for(k in 1:J){
      BB[,,k] <- BB.set[[Phi.index[k]]]
      FF[,k] <- FF.set[[Phi.index[k]]]
    }
    tFF <- t( Tau^2*t(FF) )
    for(i in 1:n){
      Us <- UU[[i]]
      if(length(Us)==0){ 
        V <- TP[i,]%*%t(TP[i,])/Sig^2 + diag(1/tFF[i,])
        prior.mean <- apply(na.omit(BB[i,,]*sBeta[NN[i,],]), 2, sum) / tFF[i,] 
        pos.mean <- TP[i,]*sY[i]/Sig^2 + prior.mean
      }
      if(length(Us)>0){
        sNN <- matrix(NN[Us,], length(Us), m)     # set of parents ID including "i" (matrix object)
        prior.prec <- matrix(0, J, J)
        for(k in 1:J){
          prior.prec[k, k] <- 1/tFF[i, k] + sum( na.omit((BB[Us,,k]^2/tFF[Us, k])[sNN==i]) )
        }
        V <- TP[i,]%*%t(TP[i,])/Sig^2 + prior.prec
        length.Us <- length(Us)
        prior.mean <- c()
        for(k in 1:J){
          mBeta <- matrix(sBeta[as.vector(sNN), k], length.Us, m)
          mBeta[sNN==i] <- 0
          mBeta[is.na(mBeta)] <- 0
          sBB <- BB[Us,,k]
          sBB[is.na(sBB)] <- 0
          bb <- na.omit( (t(sBB))[t(sNN)==i] )
          sBB[sNN==i] <- 0
          a <- sBeta[Us, k] - apply(mBeta*sBB, 1, sum)
          prior.mean[k] <- sum(bb*a/tFF[Us,k])
        }
        prior.mean <- prior.mean + apply(na.omit(BB[i,,]*sBeta[NN[i,],]), 2, sum) / tFF[i,] 
        pos.mean <- TP[i,]*sY[i]/Sig^2 + prior.mean
      }
      IV <- solve(V)
      sBeta[i,] <- mvrnorm(1, IV%*%pos.mean, IV) 
    }
    
    # Beta
    Beta <- t( t(sBeta) + Beta_center )
    if(constrain){
      for(k in 2:J){
        Beta[,k] <- (Beta[,k]>0)*Beta[,k]   # positive constraint
      }
    }
    Beta.pos[r,,] <- Beta
    
    # TP: Predictor (latent factor)
    for(j in 2:J){
      A <- 1/(Beta[,j]^2/Sig^2 + 1/VF[,j])
      resid <- Y - apply(as.matrix(TP[,-j]*Beta[,-j]), 1, sum)
      B <- Beta[,j]*resid/Sig^2 + MF[,j]/VF[,j]
      TP[,j] <- rnorm(n, A*B, sqrt(A))
    }
    
    # Sig
    resid <- Y - apply(TP*Beta, 1, sum)
    Sig <- sqrt( rinvgamma(1, prior_var+n/2, prior_var+sum(resid^2)/2) )
    Sig.pos[r] <- Sig
    
    # Tau 
    sBeta <- t( t(Beta) - Beta_center )
    prior.sBeta <- matrix(NA, n, J)
    for(k in 1:J){
      prior <- BB[,,k]*matrix(sBeta[NN,k], n, m)
      prior[is.na(prior)] <- 0
      prior.sBeta[,k] <- apply(prior, 1, sum)
    }
    ss <- apply((sBeta - prior.sBeta)^2/FF, 2, sum)
    Tau <- sqrt( rinvgamma(J, prior_var+n/2, prior_var+ss/2) )
    Tau.pos[r,] <- Tau
    
    # Phi 
    for(k in 1:J){
      new.index <- Phi.index[k] + sample(c(1, -1), 1) 
      if(new.index<1){ new.index <- 1 }
      if(new.index>L){ new.index <- L }
      new.BB <- BB.set[[new.index]]
      new.FF <- FF.set[[new.index]]
      new.prior.sBeta <- new.BB*matrix(sBeta[NN,k], n, m)
      new.prior.sBeta[is.na(new.prior.sBeta)] <- 0
      new.prior.sBeta <- apply(new.prior.sBeta, 1, sum)
      L1 <- sum( dnorm(x=sBeta[,k], mean=prior.sBeta[,k], sd=Tau[k]*sqrt(FF[,k]), log=T) )
      L2 <- sum( dnorm(x=sBeta[,k], mean=new.prior.sBeta, sd=Tau[k]*sqrt(new.FF), log=T) )
      pp <- min(1, exp(L2-L1))
      if(runif(1)<pp){ Phi.index[k] <- new.index }
    }
    Phi.pos[r,] <- Phi.set[Phi.index]
    
    # print
    if(round(10*r/mcmc)==(10*r/mcmc) & print ){ 
      print( paste0("MCMC ", 100*r/mcmc, "% completed") )
    }
  }
  
  ## omit burn-in samples
  om <- 1:burn
  Beta.pos <- Beta.pos[-om,,]
  Sig.pos <- Sig.pos[-om]
  Phi.pos <- Phi.pos[-om,]
  Tau.pos <- Tau.pos[-om,]
  
  ## Prediction 
  if(prediction){
    MC <- mcmc - burn
    MF_test <- cbind(1, mF_test)
    VF_test <- cbind(NA, vF_test)
    tau2.pos <- Tau.pos^2
    n_test <- dim(SP_test)[1]
    Pred.pos <- matrix(NA, MC, n_test)
    for(i in 1:n_test){
      dd <- sqrt( apply((SP_test[i,] - t(SP))^2, 2, sum) )
      NN <- sort(order(dd)[1:m]) 
      dd2 <- as.matrix(dist(SP[NN,]))
      dd3 <- sqrt( apply((SP_test[i,] - t(SP[NN,]))^2, 2, sum) )
      fBeta <- matrix(NA, MC, J)
      for(k in 1:J){
        for(r in 1:MC){
          C1 <- solve( tau2.pos[r,k]*exp(-dd2/Phi.pos[r,k]) )
          C2 <- tau2.pos[r,k]*exp(-dd3/Phi.pos[r,k])
          fFF <- tau2.pos[r,k] - t(C2)%*%C1%*%C2   
          fBB <- as.vector(C2%*%C1) 
          mm <- sum( fBB*(Beta.pos[r, NN, k]-Beta_center[k]) )
          fBeta[r,k] <- Beta_center[k] + rnorm(1, mm, sqrt(fFF)) 
        }
      }
      FP_test <- rep(1, MC)
      for(j in 2:J){
        FP_test <- cbind(FP_test, rnorm(MC, MF_test[i,j], sqrt(VF_test[i,j])))
      }
      Pred.pos[,i] <- apply(fBeta*FP_test, 1, sum) + Sig.pos*rnorm(MC)
    }
  }
  
  ## Summary
  if(prediction==F){ Pred.pos <- NULL }
  result <- list(beta=Beta.pos, sig=Sig.pos, phi=Phi.pos, tau=Tau.pos, pred=Pred.pos)
  return(result)
}















###----------------------------------------------------### 
###    Bayesian spatial predictive synthesis (BSPS)    ###
###       with variational Bayes approximation         ###
###----------------------------------------------------###

## INPUT
# Y: response vector 
# X: covariate matrix 
# mF: matrix of point prediction of agent (without intercept)
# vF: matrix of prediction variance of agent 
# Remark: exponential covariance function is used.  

## OUTPUT
# posterior samples of beta, sigma, tau, phi (bandwidth)

BSPS_vb <- function(Y, SP, mF, vF, maxit=100, mc=1000, Phi.set=NULL, print=T,
                    prediction=F, SP_test=NULL, mF_test=NULL, vF_test=NULL){
  # preparation
  n <- length(Y)     # the number of sampled locations
  MF <- cbind(1, mF)
  VF <- cbind(1, vF)
  J <- dim(MF)[2] 
  prior_var <- 1     # hyperparameter in inverse gamma prior for variance parameters 
  mat <- as.matrix(dist(SP))    # distance matrix
  
  # correlation function 
  CF <- function(dist.mat, phi){ exp(-dist.mat/phi) } 
  
  # function of inverse matrix
  ep <- 10^(-20)    # cut-off point 
  Inv <- function(Mat){
    Mat <- (Mat + t(Mat))/2
    dec <- eigen(Mat)
    lam <- dec$values
    lam[lam<ep] <- ep
    return((dec$vectors)%*%diag(1/lam)%*%t(dec$vectors))
  }
  
  # trace function 
  Tr <- function(x){ sum(diag(x)) }
  
  # set of spatial range parameters
  L <- length(Phi.set)
  if(is.null(Phi.set)){
    L <- 10
    Phi.range <- c(0, median(dist(SP)))
    Phi.set <- seq(Phi.range[1], Phi.range[2], length=L+1)[-1]
  }
  
  # set of spatial covariance matrices 
  L <- length(Phi.set)
  IH.set <- array(NA, c(n, n, L))
  bb.set <- c()
  for(l in 1:L){
    H <- exp(-mat/Phi.set[l])
    IH.set[,,l] <- Inv( H )
    EIG <- eigen(IH.set[,,l])$values
    EIG[EIG<ep] <- ep
    bb.set[l] <- sum( log(EIG) )
  }
  
  # variational parameters (initial values)
  Asig <- 1
  Bsig <- 1
  Atau <- rep(1, J)
  Btau <- rep(1, J)
  Mut <- matrix(0, n, J)
  Sit <- array(NA, c(n, n, J))
  PP <- matrix(1/L, L, J)
  MFt <- cbind(1, mF)
  VFt <- cbind(0, vF)
  dd <- 1    # initial value of convergence criteria 
  
  ## iteration 
  for(k in 1:maxit){
    # current values
    old.Mut <- Mut
    # mIH
    mIH <- array(0, c(n, n, J))
    for(j in 1:J){
      for(l in 1:L){
        mIH[,,j] <- mIH[,,j] + PP[l,j]*IH.set[,,l]
      }
    }
    # Beta
    for(j in 1:J){
      Om <- diag(MFt[,j]^2+VFt[,j]) 
      A <- solve( Om*Asig/Bsig + mIH[,,j]*Atau[j]/Btau[j] )
      B <- Asig/Bsig * MFt[,j] * (Y - apply(MFt[,-j]*Mut[,-j], 1, sum)) 
      Mut[,j] <- as.vector( A%*%B )
      Sit[,,j] <- A
    }
    # tau 
    for(j in 1:J){
      Atau[j] <- prior_var + n/2
      Btau[j] <- prior_var + 0.5*Tr( (Mut[,j]%*%t(Mut[,j]) + Sit[,,j])%*%mIH[,,j] )
    }
    
    # spatial range parameter
    for(j in 1:J){
      log.prob <- c()
      for(l in 1:L){
        log.prob[l] <- 0.5*bb.set[l] - 0.5*Atau[j]/Btau[j] * Tr( (Mut[,j]%*%t(Mut[,j])+Sit[,,j])%*%IH.set[,,l] )
      }
      prob <- exp(log.prob - max(log.prob))
      PP[,j] <- prob/sum(prob)
    }
    
    # sigma
    Asig <- prior_var + n/2
    resid <- Y - apply(MFt*Mut, 1, sum)
    vv <- c()
    for(j in 2:J){
      vv[j] <- Tr( (Mut[,j]%*%t(Mut[,j])+Sit[,,j]) * (MFt[,j]%*%t(MFt[,j])+diag(VFt[,j])) ) - sum(Mut[,j]^2*MFt[,j]^2)
    }
    Bsig <- prior_var + (sum(resid^2) + Tr(Sit[,,1]) + sum(vv[-1]))/2
    # F
    for(j in 2:J){
      A <- ( 1/VF[,j] + (Mut[,j]^2+diag(Sit[,,j]))*Asig/Bsig )^(-1)
      B <- MF[,j]/VF[,j] + (Mut[,j]^2+diag(Sit[,,j]))*Asig/Bsig
      MFt[,j] <-  A*B
      VFt[,j] <-  A
    }
    
    # convergence monitor 
    dd <- sum(abs(old.Mut-Mut))/sum(abs(old.Mut)+0.001)
    if(print){ print(k) }
    if(dd<0.005){ break }
  }
  
  # posterior samples
  Beta.pos <- array(NA, c(mc, n, J))
  Tau.pos <- matrix(NA, mc, J)
  Phi.pos <- matrix(NA, mc, J)
  for(j in 1:J){
    Beta.pos[,,j] <- mvrnorm(mc, Mut[,j], Sit[,,j]) 
    Tau.pos[,j] <- rinvgamma(mc, Atau[j], Btau[j])
    Phi.pos[,j] <- Phi.set[sample(1:L, mc, PP[,j], replace=T)]
  }
  Sig.pos <- rinvgamma(mc, Asig, Bsig)
  
  ## Prediction 
  if(prediction){
    MF_test <- cbind(1, mF_test)
    VF_test <- cbind(NA, vF_test)
    n_test <- dim(SP_test)[1]
    n_all <- n + n_test
    distance1 <- as.matrix(dist(SP))
    distance2 <- as.matrix(dist(rbind(SP, SP_test)))[(n+1):n_all,]
    distance2 <- distance2[,1:n]
    hPhi <- apply(Phi.pos, 2, mean)
    hBeta <- apply(Beta.pos, c(2,3), mean)
    fBeta <- matrix(NA, n_test, J)
    for(j in 1:J){
      IH <- Inv( exp(-distance1/hPhi[j]) )
      W <- exp(-distance2/hPhi[j])
      fBeta[,j] <- as.vector( W%*%IH%*%hBeta[,j] )
    }
    Pred <- apply(MF_test*fBeta, 1, sum)   # point prediction via BPS-VB
  }
  
  ## Summary
  if(prediction==F){ Pred.pos <- NULL }
  result <- list(beta=Beta.pos, sig=Sig.pos, phi=Phi.pos, tau=Tau.pos, pred=Pred)
  return(result)
}





