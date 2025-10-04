library(mvtnorm)
library(tidyverse)

# Generate design matrix
n <- 700
p <- 2000
set.seed(27)
cov <- 0.5
Sigma <- matrix(rep(cov,p*p),ncol=p)
diag(Sigma) <- rep(1,p)
mu <- rep(0,p)
X <- rmvnorm(n,mu,Sigma)
X.design <- X %*% diag(1/sqrt(diag(t(X) %*% X)/n))
write.csv(X.design,"X.designp2000n700.corr.csv",row.names = FALSE)


# Generate errors
epsilon <- rnorm(700*100,0,1)
write.csv(epsilon,"epsilon.csv",row.names = FALSE)


#Generate random beta values
set.seed(27)
betas <- runif(100,1,3)
write.csv(betas,"betas.csv",row.names = FALSE)