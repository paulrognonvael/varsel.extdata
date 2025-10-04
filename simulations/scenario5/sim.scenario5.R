library(mvtnorm)
library(tidyverse)
library(ggplot2)
library(mombf)
library(stringr)
library(glmnet)
library(ncvreg)
library(parallel)

# Working directory
setwd(dirname(utils::getSrcFilename(function(){}, full.names = TRUE)))

# Helper functions
source(paste0(path,"routines.R"))

# Scenario

nb.inactive2 <-function(n){
  return(sqrt(n)/2)
}

nb.inactive3 <-function(n){
  return(n/2-sqrt(n)/2)
}

nb.inactive.tot<-function(n,nb.inactive.b0=nb.inactive3,nb.inactive.b1=nb.inactive2){
  nb.inactive.b0(n)+nb.inactive.b1(n)
}

nb.active<- function(n){
  return(2 * round(3*log(n)/2)) # rounding to the closest even number
}

bmin <- 0.2


# Load simulated data
X <- as.matrix(read.csv(paste0(dirname(getwd()),"/X.designp2000n700.corr.csv")))
epsilon <- as.matrix(read.csv(paste0(dirname(getwd()),"/epsilon.csv")))
betas <- as.matrix(read.csv(paste0(dirname(getwd()),"/betas.csv")))[,1]

# Lists of methods
l0method.vec <- c("kappa.o", "EBIC", "S.EB", "S.A", "S.EB.b",
                  "S.A.b")

method.vec <- c(l0method.vec,"lasso.cv", "scad.cv")

###################################################

values.n <- c(20,40,60,80,seq(100,700,100))

for (n  in values.n){
  t0 <- Sys.time()
  betamin <- bmin
  beta_star0 <- c(0.8, betas[1:(nb.active(n)/2-1)], rep(0,nb.inactive3(n)))
  beta_star1 <- c(betamin, betas[(nb.active(n)/2):(nb.active(n)-2)], rep(0,nb.inactive2(n)))
  beta_star <- c(beta_star0, beta_star1)
  block0 <- c(1:length(beta_star0))
  block1 <- c((length(beta_star0)+1):(length(beta_star0)+length(beta_star1)))
  
  X.design <- X[1:n,1:length(beta_star)]
  
  # Simulation
  
  m <- 100 #number of simulations 
  # K <- 10 #number of folds for CV
  
  y.df <- replicate(m,(X.design%*%(beta_star))[,1])+epsilon[1:(n*m),1]
  
  sim.result.sel <- data.frame()
  sim.result.est.mse <- data.frame()
  #sim.result.cvmse <- data.frame()
  for(i in 1:m){
    res.sel.l0 <- selectionl0.comp(y.df[,i],X.design=X.design, block0=block0, block1=block1, beta_star= beta_star)
    res.est.mse.l0 <- sapply(res.sel.l0[,l0method.vec], mse.est, y=y.df[,i],X.design=X.design, beta_star= beta_star)
    res.est.mse.l0 <- data.frame(t(res.est.mse.l0))
    colnames(res.est.mse.l0) <- l0method.vec
    res.lasso.scad <- sel.lasso.scad.cv(y.df[,i],X.design=X.design, beta_star= beta_star)
    
    res.sel <- data.frame(res.sel.l0,res.lasso.scad$sel,sim=i,n=n,betamin=betamin)
    sim.result.sel <- bind_rows(sim.result.sel,res.sel)
    
    res.est.mse <- data.frame(res.est.mse.l0,res.lasso.scad$est.mse,sim=i,n=n,betamin=betamin)
    sim.result.est.mse <- bind_rows(sim.result.est.mse,res.est.mse)
    
    
    # res.mse.l0 <- cvmse.l0.comp(y.df[,i], X.design, block0, block1, K, mc.cores=4)
    # res.mse <- data.frame(res.mse.l0,res.lasso.scad$cvmse,sim=i,n=n,betamin=betamin)
    # sim.result.cvmse <- bind_rows(sim.result.cvmse,res.mse)
  }
  
  
  sim.result.df <- pivot_longer(sim.result.sel,cols=method.vec,names_to = 'method',values_to = 'sel.model')
  
  sim.result.est.mse.df <- pivot_longer(sim.result.est.mse,cols=method.vec,names_to = 'method',values_to = 'est.mse')
  sim.result.df <- merge(sim.result.df, sim.result.est.mse.df, by = c('method','sim','n','betamin'))
  
  # sim.result.cvmse.df <- pivot_longer(sim.result.cvmse,cols=method.vec,names_to = 'method',values_to = 'cv.mse')
  # sim.result.df <- merge(sim.result.df, sim.result.cvmse.df, by = c('method','sim','n','betamin'))
  
  sim.result.df <- postprocess.sim(sim.result.df,beta_star)
  
  write.csv(sim.result.df,paste0("sim.result.scenario5.n",n,".csv"),row.names = FALSE)
  
  t1 <- Sys.time() 
  cat('Time n=',n,':', round(difftime(t1, t0, units = "mins"),3),'minutes'); cat('\n');
  
}