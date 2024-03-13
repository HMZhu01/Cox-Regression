## Projection of Gaussian × Gaussian-Bernoulli
GB_Proj <- function(theta, mu, sigma, m2, v2){
  
  
  et <- (1-theta)*(sqrt(v2/(sigma+v2)) * exp(-0.5*( (m2)^2 / v2 - (mu - m2)^2 / (sigma+v2) )))/theta
  Ex <- ((mu/sigma+m2/v2)/(1/sigma+1/v2)) * ( 1/ (1+et) )
  Ex2 <- (1/sigma+1/v2+(mu/sigma+m2/v2)^2) / (1/sigma+1/v2)^2 * ( 1/ (1+et) );
  
  #const <- theta*dnorm(0,mu-m2,sigma+v2)+(1-theta)*dnorm(0,m2,v2)
  #Ex <- theta*dnorm(0,mu-m2,sigma+v2)*(mu/sigma+m2/v2)/const/(1/sigma+1/v2)
  #Ex2 <- theta*dnorm(0,mu-m2,sigma+v2)*(1/sigma+1/v2+(mu/sigma+m2/v2)^2)/const/(1/sigma+1/v2)^2
  Dx <- Ex2-Ex^2
  
  ## Output
  output <- list(m=Ex, v=mean(Dx))
  return(output)
}

## Projection of Gaussian × Cox partial likelihood
Cox_Proj <- function(m, v, c, M, MaxIter=1000, esc=1e-20){
  ## Initialization
  hatz <- -1*rep(1, M)
  upM <- matrix(rep(1, M*M), M, M)
  lowM <- matrix(rep(1, M*M), M, M)
  upM[!upper.tri(upM, diag = TRUE)] <- 0
  lowM[!lower.tri(lowM, diag = TRUE)] <- 0
  hatz_old <- hatz

  ## Expectation
  for(i in 1:MaxIter){
    ld <- c-(upM %*% ((lowM %*% exp(hatz))^(-1)*c))*exp(hatz)
    ldd <- exp(2*hatz)*(upM %*% ((lowM %*% exp(hatz))^(-2)*c))  - exp(hatz)*(upM %*% ((lowM %*% exp(hatz))^(-1)*c))
    index <- (ldd==0)
    ldd[index] <- 0.00000000001
    if (i > 1){
      hatz_old <- hatz
    }
    zt <- hatz-ld/ldd
    hatz <-  zt+(zt-m)/(v*ldd-1)
    a <- sum((hatz_old-hatz)^2)/sum(hatz^2)
    if (a < esc){
      break
    }
  }
  
  ## Variance
  Q <- rep(0, M)
  z1 <- t(hatz)
  cpx1 <- c/t(exp(z1)%*%upM)
  cpx2 <- c/t((exp(z1)%*%upM)^2)
  for(k in 1:M){
    sum1 <- 0
    sum2 <- 0
    for(s in k:M){
      sum1 <- sum1 + cpx1[s]
      sum2 <- sum2 + cpx2[s]
    }
    Q[k] <-  (exp(z1[k])*sum1 + exp(2*z1[k])*sum2)
  }
  
  ## Output
  output <- list(m=hatz, v=mean((1/(1/v + Q))) )
  return(output)
}


## Projection of Gaussian × Laplace-Bernoulli
LB_Proj <- function(theta, lambda, m, v){
  eta <- pnorm((m-lambda*v)/sqrt(v))*exp(-lambda*m) + pnorm((-m-lambda*v)/sqrt(v))*exp(lambda*m) + ((1-theta)/(theta*lambda))*sqrt(2/(pi*v))*exp(-(m^2+lambda^2*v^2)/(2*v))   
  
  m1 <- ((m+v*lambda)*exp(m*lambda)*pnorm((-m-v*lambda)/sqrt(v)) +(m-v*lambda)*exp(-m*lambda)*pnorm((m-v*lambda)/sqrt(v))) / (eta)
  
  v1 <- (-2*v^2*lambda/sqrt(2*pi*v)*exp(-(m^2+lambda^2*v^2)/(2*v)) + (m-v*lambda)^2*pnorm((m-lambda*v)/sqrt(v))*exp(-2*lambda*m)+ (m+v*lambda)^2*pnorm((-m-lambda*v)/sqrt(v))*exp(2*lambda*m) + v*exp(-lambda*m)*pnorm((m-lambda*v)/sqrt(v)) + v*exp(lambda*m)*pnorm((-m-lambda*v)/sqrt(v))  ) /(eta) - m1^2

  ## Output
  output <- list(m=m1, v=mean(v1))
  return(output)
}