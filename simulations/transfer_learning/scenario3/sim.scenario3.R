library(mvtnorm)
library(tidyverse)
library(ggplot2)
library(mombf)
library(stringr)
library(glmnet)
library(ncvreg)
library(parallel)

set.seed(1278)

# Working directory
path = "C:/Users/Rognon/Documents/GitHub/varsel.extdata/simulations/"
setwd(paste0(path,'transfer_learning/scenario3/'))

# Helper functions
source(paste0(path,"routines.R"))

# Scenario

nb.inactive2 <-function(n){
  return(sqrt(n))
}

nb.inactive3 <-function(n){
  return(1.5*n-sqrt(n))
}

nb.inactive.tot<-function(n,nb.inactive.b0=nb.inactive3,nb.inactive.b1=nb.inactive2){
  return(nb.inactive.b0(n)+nb.inactive.b1(n))
}

nb.active<- function(n){
  return(2 * round(3*log(n)/2)) # rounding to the closest even number
}

betamin <- 0.33

# Load simulated data
X <- as.matrix(read.csv(paste0(path,"X.designp2000n700.corr.csv")))
epsilon <- as.matrix(read.csv(paste0(path,"epsilon.csv")))
betas <- as.matrix(read.csv(paste0(path,"betas.csv")))[,1]

# Lists of methods
l0method.vec <- c("kappa.o", "EBIC", "S.EB", "S.A", "S.EB.b",
  "S.A.b")

method.vec <- c(l0method.vec,"lasso.cv", "scad.cv")

# Number of source datasets
K.source <- 2 


# Perturbation frequency
h <- 0.4

###################################################

values.n <- c(20,40,60,80,seq(100,700,100))


for (n  in values.n){
  cat('n=',n,'\n')
  t0 <- Sys.time() 
  beta_star <- c(betamin,betas[1:(nb.active(n)-1)], rep(0,nb.inactive.tot(n)))
  
  X.design <- X[1:n,1:length(beta_star)]
  
  # Simulation
  
  m <- 100 #number of simulations 
  
  y.df <- replicate(m,(X.design%*%(beta_star))[,1])+epsilon[1:(n*m),1]
  
  sim.result.sel <- data.frame()
  sim.result.est.mse <- data.frame()
  #sim.result.cvmse <- data.frame()
  
  ## Sources design matrix 
  X.design.source <- X[1:700,1:length(beta_star)]
  source_X <- replicate(K.source, X.design.source, simplify = FALSE)
  
  ## Number of perturbations
  nt1 <- nt2 <- 2*round(h*nb.active(n)/2)
  
  for(i in 1:m){
    cat('   iteration ',i,'\n')
    source_y <- vector(mode = "list", length = K.source)
    
    for(k in 1:K.source){
      #  Initialize beta k
      beta.k <- beta_star
      
      # Perturbation set
      H.t1 <- sample((nb.active(n)+1):(2*nb.active(n)), size=nt1)
      H.t2 <-sample(1:nb.active(n), size=nt2)
      
      # Generate beta k
      beta.k[H.t1] <- runif(nt1, 0.5, 1) 
      beta.k[H.t2] <- rep(0, nt2) 
      
      source_y[[k]] <- (X.design.source%*%(beta.k))[,1] + rnorm(700,0,1)
    }
    
    res.sel.l0 <- trans.s.ell0(target_y = y.df[,i], target_X = X.design, source_y=source_y, 
                               source_X = source_X, beta_star=beta_star)
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
  
  sim.result.df <- postprocess.sim(sim.result.df=sim.result.df, beta_star= beta_star)
  
  write.csv(sim.result.df,paste0("sim.result.scenario3.n",n,".csv"),row.names = FALSE)
  
  t1 <- Sys.time() 
  cat('Time n=',n,':', round(difftime(t1, t0, units = "mins"),3),'minutes'); cat('\n');
  
}