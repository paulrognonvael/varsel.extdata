library(mvtnorm)
library(tidyverse)

setwd("~/GitHub/varsel.extdata/simulations")

# Generate design matrixes 
n <- 700
p <- 2000
set.seed(27)
cov <- 0.5
Sigma <- matrix(rep(cov,p*p),ncol=p)
diag(Sigma) <- rep(1,p)
mu <- rep(0,p)

# Target data
X <- rmvnorm(n,mu,Sigma)
X.design <- X %*% diag(1/sqrt(diag(t(X) %*% X)/n))
write.csv(X.design,"X.designp2000n700.corr.csv",row.names = FALSE)

Sigma_source = matrix(NA, nrow=p, ncol=p)
for(i in 1:2000){
  for(j in 1:2000){
    Sigma_source[i,j] = 0.5^(abs(i-j))
  }
}

# Source data 1
X <- rmvt(n, sigma = Sigma_source, df= 8)
X.design <- X %*% diag(1/sqrt(diag(t(X) %*% X)/n))
write.csv(X.design,"X.source1.csv",row.names = FALSE)
# Source data 2
X <- rmvt(n, sigma = Sigma_source, df= 8)
X.design <- X %*% diag(1/sqrt(diag(t(X) %*% X)/n))
write.csv(X.design,"X.source2.csv",row.names = FALSE)


# Generate errors
epsilon <- rnorm(700*100,0,1)
write.csv(epsilon,"epsilon.csv",row.names = FALSE)


#Generate random beta values
betas <- runif(100,1/2,1)
write.csv(betas,"betas.csv",row.names = FALSE)

#Generate random beta values between 1/3 and 1
betas <- runif(100,1/3,1)
write.csv(betas,"betas2.csv",row.names = FALSE)