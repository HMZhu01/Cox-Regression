EM_GB <- function(theta,mu,sigma,m0_sub,v0_sub,N){
  ## Expectation Maximum
  pi <- 1/(1+((1-theta)/theta)*(dnorm(0,m0_sub,v0_sub)/dnorm(mu,m0_sub,sigma+v0_sub)))
  gamma <- ( (m0_sub/v0_sub) + (mu/sigma) )/( (1/v0_sub) + (1/sigma) )
  nu <- (v0_sub*sigma)/(sigma+v0_sub)
  
  theta_new <- sum(pi)/N
  mu_new <- sum(pi*gamma)/(theta_new*N)
  sigma_new <- sum( pi*( (abs(mu-gamma))^2 + nu ) )/(theta_new*N)
  
  output <- list(theta=theta_new, mu=mu_new, sigma=sigma_new)
  return(output)
}

EM_LB <- function(theta,lambda,m,v){
  ## Expectation Maximum
  tau1 <- (m-lambda*v)/sqrt(v)
  tau2 <- (-m-lambda*v)/sqrt(v)
  eta <- pnorm(tau1)*exp(-lambda*m) + pnorm(tau2)*exp(lambda*m) + ((1-theta)/(theta*lambda))*sqrt(2/(pi*v))*exp((m^2+lambda^2*v^2)/(-2*v))   
  
  return( mean( ((pnorm(tau1)*exp(-lambda*m) + pnorm(tau2)*exp(lambda*m) ) /eta) ))
}