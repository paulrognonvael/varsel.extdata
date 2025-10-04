library(mvtnorm)
library(tidyverse)
library(ggplot2)
library(mombf)
library(stringr)
library(glmnet)
library(future.apply)
library(ncvreg)
library(parallel)

# Working directory
setwd(dirname(utils::getSrcFilename(function(){}, full.names = TRUE)))

# Parallelizing
plan(multisession)


# Helper functions
source(paste0(dirname(getwd()),"/routines.R"))

# Scenario

nb.inactive2 <-function(n){
  return(sqrt(n))
}

nb.inactive3 <-function(n){
  return(1.5*n-sqrt(n))
}

nb.inactive.tot<-function(n,nb.inactive.b0=nb.inactive3,nb.inactive.b1=nb.inactive2){
  nb.inactive.b0(n)+nb.inactive.b1(n)
}

nb.active<- function(n){
  return(3*log(n))
}

bmin <- 0.33

# Load simulated data
X <- as.matrix(read.csv(paste0(dirname(getwd()),"/X.designp2000n700.corr.csv")))
epsilon <- as.matrix(read.csv(paste0(dirname(getwd()),"/epsilon.csv")))
betas <- as.matrix(read.csv(paste0(dirname(getwd()),"/betas.csv")))[,1]

###################################################

values.n <- c(20,40,60,80,seq(100,700,100))

iter <- 1
for (n  in values.n){
  t0 <- Sys.time() 
  niter <- length(values.n)
  betamin <- bmin
  beta_star0 <- c(0.8, betas[1:max(round(nb.active(n)/2-1),1)], rep(0,nb.inactive3(n)))
  beta_star1 <- c(betamin, betas[(max(round(nb.active(n)/2-1),1)+1):round(nb.active(n)-3)], rep(0,nb.inactive2(n)))
  beta_star <- c(beta_star0, beta_star1)
  block0 <- c(1:length(beta_star0))
  block1 <- c((length(beta_star0)+1):(length(beta_star0)+length(beta_star1)))
  
  X.design <- X[1:n,1:length(beta_star)]
  
  # Simulation
  
  m <- 100 #number of simulations 
  K <- 10 #number of folds for CV
  
  y.df <- replicate(m,(X.design%*%(beta_star))[,1])+epsilon[1:(n*m),1]
  
  # sim.result <- future_apply(y.df,2,selection.comp,X.design=X.design, block0=block0, block1=block1, beta_star= beta_star, future.seed=TRUE)
  # sim.result.df <- bind_rows(sim.result)
  sim.result.sel <- data.frame()
  sim.result.cvmse <- data.frame()
  
  res.sel.l0 <- future_apply(y.df,2,selectionl0.comp,X.design=X.design, block0=block0, block1=block1, beta_star= beta_star)
  res.mse.l0 <- future_apply(y.df,2, cvmse.l0.comp, X.design, block0, block1, K, mc.cores=4)
  res.lasso.scad <- future_apply(y.df,2,sel.lasso.scad.cv,X.design=X.design,future.seed=TRUE)
  res.lasso.scad.sel <- lapply(res.lasso.scad, function(l) l[[1]])
  res.lasso.scad.mse <- lapply(res.lasso.scad, function(l) l[[2]])
  
  sim.result.sel <- cbind(bind_rows(res.sel.l0),bind_rows(res.lasso.scad.sel))
  sim.result.cvmse <- cbind(bind_rows(res.mse.l0),bind_rows(res.lasso.scad.mse))
  
  sim.result.sel['sim'] <- 1:100
  sim.result.sel['n'] <- n
  sim.result.sel['betamin'] <- bmin
  
  sim.result.cvmse['sim'] <- 1:100
  sim.result.cvmse['n'] <- n
  sim.result.cvmse['betamin'] <- bmin
  
  method.vec <- c(
    #"kappa.o",
    "EBIC",
    "S.EB",
    "S.A",
    "S.EB.b",
    "S.A.b",
    "lasso.cv",
    "scad.cv"
  )
  
  sim.result.df <- pivot_longer(sim.result.sel,cols=method.vec,names_to = 'method',values_to = 'sel.model')
  
  sim.result.mse.df <- pivot_longer(sim.result.cvmse,cols=method.vec,names_to = 'method',values_to = 'cv.mse')
  
  sim.result.df <- merge(sim.result.df, sim.result.mse.df, by = c('method','sim','n','betamin'))
  
  sim.result.df['recovery']<-  sim.result.df$sel.model==paste(as.character(which(beta_star!=0)),collapse=",")
  sel.list <- lapply(sim.result.df$sel.model, model.char2vec)
  sim.result.df['nb.tI.b0'] <- sapply(sel.list, nb.error, true.beta=beta_star, type='I', block=block0)
  sim.result.df['nb.tI.b1'] <- sapply(sel.list, nb.error, true.beta=beta_star, type='I', block=block1)
  sim.result.df['nb.tII.b0'] <- sapply(sel.list, nb.error, true.beta=beta_star, type='II', block=block0)
  sim.result.df['nb.tII.b1'] <- sapply(sel.list, nb.error, true.beta=beta_star, type='II', block=block1)
  sim.result.df['est.size'] <- sapply(sel.list,length)
  
  write.csv(sim.result.df,paste0("sim.result.scenario1.n",n,"_future.csv"),row.names = FALSE)
  
  t1 <- Sys.time() 
  cat('Time n=',n,':', round(difftime(t1, t0, units = "mins"),3),'minutes'); cat('\n');
  
}
