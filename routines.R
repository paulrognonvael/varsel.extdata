# Converts a string of variable indices separated by commas into a vector
model.char2vec <- function(model.char){
  return(as.numeric(strsplit(model.char,",")[[1]]))
}


## Returns the number of false positive and false negatives

nb.error<-function(sel.vect, true.beta, type, block){
  truth <- true.beta != 0
  sel.model <- rep(FALSE, length(true.beta))
  sel.model[sel.vect] <- TRUE
  
  if(type=='I'){
    tI <- truth==FALSE & sel.model==TRUE
    return(sum(tI[block]))
  } else if (type=='II'){
    tII <- truth==TRUE & sel.model==FALSE
    return(sum(tII[block]))
  } else{
    stop("Invalid type of error")
  }
  
}

# Unused
get.model <- function(glmnet.beta){
  return(paste(unname(which(glmnet.beta!=0)),collapse=","))
}

## Returns the l0 scores for a given model

l0score <- function(model, y, X.design, block0, block1, s0.hat, s1.hat, bmin0.hat, bmin1.hat,true.bmin){
  beta.hat <- rep(0,ncol(X.design))
  k <- 0
  k0 <-0
  k1 <- 0
  n <- length(y)
  if (sum(is.na(model))!=length(model)){
    model <- model[!is.na(model)]
    X.model <- X.design[,model]
    beta.hat[model] <- solve(t(X.model)%*%X.model)%*%t(X.model) %*% y
    k0 <- sum(beta.hat[block0] !=0)
    k1 <- sum(beta.hat[block1] !=0)
    k <- sum(beta.hat !=0)
  }
  
  s.hat <- s0.hat + s1.hat
  
  bmin.hat <- min(bmin0.hat,bmin1.hat)
  
  pred <- X.design%*%beta.hat
  
  # kappa^{0}
  kappa.o <- 0.5*n*log(sum((y-pred)^2)/n) + k*0.5*log(n) + k*log(ncol(X.design))
  
  # EBIC
  ebic <- 0.5*n*log(sum((y-pred)^2)/n) + 0.5*k*log(n) + log(choose(ncol(X.design),k))
  
  # S^EB, S^A
  S.EB <- 0.5*n*log(sum((y-pred)^2)/n) + 0.5*k*log(n) + k*log(ncol(X.design)/s.hat-1)
  S.A <- 0.5*n*log(sum((y-pred)^2)/n) + 0.5*k*log(n) + k*log(ncol(X.design)-s.hat)
  
  # S^EB,b
  S.EB.b <- 0.5*n*log(sum((y-pred)^2)/n) + 0.5*k*log(n) + k0*log(ncol(X.design[,block0])/s0.hat-1) + k1*log(ncol(X.design[,block1])/s1.hat-1)
  
  # S^A,b
  S.A.b <- 0.5*n*log(sum((y-pred)^2)/n) + 0.5*k*log(n) + k0*log(ncol(X.design[,block0])-s0.hat) + k1*log(ncol(X.design[,block1])-s1.hat)
  
  #oracle <- 0.5*n*log(sum((y-pred)^2)/n) + k*(0.5*true.bmin+log(ncol(X.design)/5-1)/(n*true.bmin))^2*(n/2)
  #block.oracle <- 0.5*n*log(sum((y-pred)^2)/n) + k0*(0.5*true.bmin+log(ncol(X.design[,block0])/4-1)/(n*true.bmin))^2*(n/2) + k1*(0.5*3+log(ncol(X.design[,block1])/1-1)/(n*3))^2*(n/2)
  
  models.val <- paste(model, collapse = "," )
  
  return(c(models.val,kappa.o, ebic, S.EB,
           S.A, S.EB.b, S.A.b))
}

## Returns the l0 scores for all models in a given list

compute.l0score<-function(models, y, X.design, block0, block1, s0.hat, s1.hat, bmin0.hat, bmin1.hat, true.bmin=NULL){
  l0score.df<- apply(models,1,l0score,y=y,X.design=X.design,
                     block0=block0, block1=block1, s0.hat=s0.hat, s1.hat=s1.hat,
                     bmin0.hat = bmin0.hat, bmin1.hat = bmin1.hat, true.bmin = true.bmin)
  l0score.df <-data.frame(t(l0score.df))
  l0score.df[,-1] <- sapply(l0score.df[, -1], as.numeric)
  colnames(l0score.df)<-c('model','kappa.o','EBIC','S.EB','S.A',
                          'S.EB.b', 'S.A.b')
  return(l0score.df)
}


## Returns model selected by S^EB, S^A, S^EB,b, S^A,b and EBIC 

selectionl0.comp <-function(y, X.design, block0, block1, beta_star=NULL){
  
  t0 <- Sys.time()
  
  # Get top models with kappa.o
  best.kappa.o.res <- mombf::bestIC(y,X.design, penalty= log(nrow(X.design))+2*log(ncol(X.design)), 
                              verbose=FALSE, maxvars= nrow(X.design)
                              #vmaxvars=nrow(X.design)*0.9,
                              # deltaini=ini,
  )
  sel.kappa.o <- best.kappa.o.res$models$modelid[1]
  top.kappa.o.list <- distinct(best.kappa.o.res$models[,1])
  
  # Add models visited by lasso
  glmnet.res<- cv.glmnet(X.design,y,intercept = FALSE)
  models.glmnet <- data.frame(modelid=apply(glmnet.res$glmnet.fit$beta,2,get.model))
  top.kappa.o.list<- distinct(rbind(top.kappa.o.list,models.glmnet))
  
  # Format list of models
  max.model.size <- max(sapply(str_split(top.kappa.o.list$modelid, ","),length))
  top.kappa.o.fmted <- str_split_fixed(top.kappa.o.list$modelid, ",",max.model.size)
  top.kappa.o.fmted[top.kappa.o.fmted == ""] <- NA
  class(top.kappa.o.fmted) <- "numeric"
  # sanity check can't have models with more than n variables
  top.kappa.o.fmted <- top.kappa.o.fmted[,1:min(nrow(X.design)-1,ncol(top.kappa.o.fmted))]
  #write.csv(topIC.lasso.fmted, 'topmodels.csv')
  
  # Estimate s, s0.hat, s1.hat
  margprob <- best.kappa.o.res$msfit$margpp
  s.hat <- sum(margprob,na.rm=TRUE)
  margprob0 <- margprob[block0]
  margprob1 <- margprob[block1]
  s0.hat <- sum(margprob0,na.rm=TRUE)
  s1.hat <- sum(margprob1,na.rm=TRUE)
  
  # Cleaning
  rm(best.kappa.o.res)
  
  # Estimate bmin0, bmin1 (not used)
  bmin0.hat <- NA
  bmin1.hat <- NA
  
  ## Sanity checks
  if(!is.null(beta_star)){
    truth <- beta_star!=0 
    truth.char <- paste(as.character(which(beta_star!=0)),collapse=",")
    is.true.model.intop.kappa.o <- truth.char %in% top.kappa.o.list$modelid
    is.true.model.intop.kappa.oplusLASSO <- truth.char %in% top.kappa.o.list$modelid
    is.s.hat.nan <- is.na(sum(margprob))
    true.bmin <- min(beta_star[beta_star!=0])
  } else {
    is.true.model.intop.kappa.o <- NA
    is.true.model.intop.kappa.oplusLASSO <- NA
    is.s.hat.nan <- NA
  }

  
  ## selected models
  
  l0score.df <- compute.l0score(top.kappa.o.fmted, y, X.design, block0, block1, s0.hat, s1.hat, bmin0.hat, bmin1.hat, true.bmin)
  
  t1 <- Sys.time()
  
  time.secs <- round(difftime(t1, t0, units = "secs"),3)
  
  sel.kappa.o <- l0score.df[which.min(l0score.df$kappa.o),"model"]
  sel.ebic <- l0score.df[which.min(l0score.df$EBIC),"model"]
  sel.S.EB <- l0score.df[which.min(l0score.df$S.EB),"model"]
  sel.S.A <- l0score.df[which.min(l0score.df$S.A),"model"]
  sel.S.EB.b <- l0score.df[which.min(l0score.df$S.EB.b),"model"]
  sel.S.A.b   <- l0score.df[which.min(l0score.df$S.A.b),"model"]
  
  
  ## formating results 
  method.vec <- c(
    "kappa.o",
    "EBIC",
    "S.EB",
    "S.A",
    "S.EB.b",
    "S.A.b"
  )
  
  sel.model.res <-c(
    sel.kappa.o,
    sel.ebic,
    sel.S.EB,
    sel.S.A,
    sel.S.EB.b,
    sel.S.A.b
  )
  
  res.all <- c(sel.model.res,
               is.true.model.intop.kappa.o,
               is.true.model.intop.kappa.oplusLASSO,
               is.s.hat.nan,
               s.hat,
               s0.hat,
               s1.hat,
               time.secs)
  
  df <- data.frame(t(res.all))
  colnames(df) <- c(method.vec,
                    'is.true.model.intop.kappa.o',
                    'is.true.model.intop.kappa.oplusLASSO',
                    'is.s.hat.nan',
                    's.hat',
                    's0.hat',
                    's1.hat',
                    'time.l0')
  
  gc()
  
  return(df)
}

## Returns MSE in coef estimates for models selected by l0 penalties
mse.est <- function(model,y, X.design, beta_star){
  
  beta.est <- rep(0,length(beta_star))
  X.design.df <- data.frame(X.design)
  model.vec <- model.char2vec(model)
  data<-cbind(data.frame(y=y),X.design.df[,model.vec,drop=FALSE])
  
  lin.reg <- lm(y~.+0,data=data)
  beta.est[model.vec] <- coef(lin.reg)
  
  mse <- sum((beta.est-beta_star)^2)/length(beta_star) 
  return(mse)
}

## Returns CV MSE with l0 penalties.

cvmse.l0.comp <- function(y, X.design, block0, block1, K, mc.cores=1){
  if (K > nrow(X.design)) stop("The number of folds cannot be larger than nrow(x)")
  subset <- rep(1:K,ceiling(nrow(X.design)/K))[1:nrow(X.design)]
  subset <- sample(subset,size=nrow(X.design),replace=FALSE)
  f <- function(k,...) {
    sel <- subset==k
    
    sel.l0 <- selectionl0.comp(y=y[!sel],X.design[!sel,,drop=FALSE], 
                               block0, block1)[,c("EBIC",
                                                   "S.EB",
                                                   "S.A",
                                                   "S.EB.b",
                                                   "S.A.b")]
    
    sel.list <- lapply(sel.l0, model.char2vec)
    
    pred.f <- function(model,y, X.design, sel){
      X.design.df <- data.frame(X.design)
      data<-cbind(data.frame(y=y[!sel]),X.design.df[!sel,model,drop=FALSE])
      lin.reg <- lm(y~.+0,data=data)
      pred <- predict(lin.reg, newdata= X.design.df[sel,model,drop=FALSE])
      return(pred)
    }
    
    pred.subset <- sapply(sel.list, pred.f, y, X.design, sel)
    
   }
  
  if (mc.cores > 1) {
    if ("parallel" %in% loadedNamespaces())  {
      allpred <- parallel::mclapply(1:K, f, mc.preschedule=FALSE)
    } else {
      stop("Did not find mclapply. Please load parallel package")
    }
  } else {
    allpred <- lapply(1:K, f)
  }
    
  pred= data.frame(matrix(NA,nrow=nrow(X.design), ncol= length(c("EBIC",
                                      "S.EB",
                                      "S.A",
                                      "S.EB.b",
                                      "S.A.b"))))
  for (k in 1:K) {
      pred[subset==k,] <- allpred[[k]]
  }
    
  res.cvmse <- colMeans((pred-y)^2)
  df.cvmse <- data.frame(t(res.cvmse))

  res.rsquared <- (apply(pred,2,cor,y))^2
  df.rsquared <- data.frame(t(res.rsquared))
  
  colnames(df.cvmse) <- c("EBIC", "S.EB", "S.A", "S.EB.b","S.A.b")
  colnames(df.rsquared) <- c("EBIC", "S.EB", "S.A", "S.EB.b","S.A.b")
  return(list(rsquared=df.rsquared,cvmse=df.cvmse))
  }


## Returns model selected, estimation MSE and CV MSE with CV LASSO and SCAD

sel.lasso.scad.cv <- function(y, X.design, beta_star=NULL){
  
  ## CV-LASSO 
  glmnet.res<- cv.glmnet(X.design,y,intercept = FALSE)
  sel.lasso.cv <- paste(which(coef(glmnet.res, s=glmnet.res$lambda.min)[-1]!=0),collapse=",")
  cv.mse.lasso <- glmnet.res$cvm[which(glmnet.res$lambda==glmnet.res$lambda.min)]
  
  ## SCAD
  cvscad <- cv.ncvreg(X=X.design,y=y,family="gaussian",penalty="SCAD",nfolds=10,dfmax=1000,max.iter=10^4)
  sel.scad.cv <- paste(which(cvscad$fit$beta[-1,cvscad$min]!=0),collapse=",")
  cv.mse.scad <- cvscad$cve[which(cvscad$lambda==cvscad$lambda.min)]
  
  df.sel <- data.frame(t(c(sel.lasso.cv, sel.scad.cv))) 
  df.cvmse <- data.frame(t(c(cv.mse.lasso, cv.mse.scad)))
  colnames(df.sel) <- c("lasso.cv","scad.cv")
  colnames(df.cvmse) <- c("lasso.cv","scad.cv")
  
  if(!is.null(beta_star)){
    est.mse.cvlasso <- sum((coef(glmnet.res, s=glmnet.res$lambda.min)[-1]-beta_star)^2)/length(beta_star)
    est.mse.cvscad <- sum((cvscad$fit$beta[-1,cvscad$min]-beta_star)^2)/length(beta_star)
    df.est.mse <- data.frame(t(c(est.mse.cvlasso, est.mse.cvscad)))
    colnames(df.est.mse) <- c("lasso.cv","scad.cv")
    return(list(sel=df.sel,cvmse=df.cvmse,est.mse = df.est.mse))
  }
  
  return(list(sel=df.sel,cvmse=df.cvmse))
}

postprocess.sim <- function(sim.result.df,beta_star){
  sim.result.df['recovery']<-  sim.result.df$sel.model==paste(as.character(which(beta_star!=0)),collapse=",")
  sel.list <- lapply(sim.result.df$sel.model, model.char2vec)
  sim.result.df['nb.tI.b0'] <- sapply(sel.list, nb.error, true.beta=beta_star, type='I', block=block0)
  sim.result.df['nb.tI.b1'] <- sapply(sel.list, nb.error, true.beta=beta_star, type='I', block=block1)
  sim.result.df['nb.tII.b0'] <- sapply(sel.list, nb.error, true.beta=beta_star, type='II', block=block0)
  sim.result.df['nb.tII.b1'] <- sapply(sel.list, nb.error, true.beta=beta_star, type='II', block=block1)
  sim.result.df['est.size'] <- sapply(sel.list,length)
  sim.result.df <- sim.result.df %>% mutate(FP=nb.tI.b0+nb.tI.b1, TP = sum(beta_star!=0)-(nb.tII.b0+nb.tII.b1), 
  FDR = FP/(FP+TP), power = TP/sum(beta_star!=0))
  return(sim.result.df)
}
