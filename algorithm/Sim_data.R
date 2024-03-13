generateData_laplacePiror <- function(mu,lambda,N){
  ## Random number sampled from Laplace distributions
  t <- runif(N,0,1)-0.5
  X <-  mu - (1/lambda)*sign(t)*log(1-2*abs(t))
  return(X)
}

SurvData <- function(M, N, lambda, snr, perc, rho, mu=0){
  ## Generate simulation data
  x <-  generateData_laplacePiror(mu,lambda,N)
  index <- runif(N, 0, 1)<rho
  x <- x*index
  
  H <- matrix(rnorm(M * N), nrow=M)
  st <- log( runif(M, 0, 1)) / -exp(H%*%x)
  ct <- Inf*rep(1, M)
  d <- runif(M, 0, 1)
  index <- d<perc
  ct[index] <- runif(M,0,st)[index]
  y <- pmin(st,ct)
  c <-  !index
  Output <- list(X=H,y=y,c=c,beta=x)
  return(Output)
}


