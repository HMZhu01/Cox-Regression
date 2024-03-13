source("./algorithm/Pre_processing.R")
source("./algorithm/Projection.R")
source("./algorithm/Picd.R")
source("./algorithm/CI.R")
source("./algorithm/EM.R")

damping <- function(x, x_old, mes){
  X <- mes*x+(1-mes)*x_old
  X_old <- x
  out <- list(x=X, x_old=X_old)
  return(out)
}

EP <- function(X_pre, y_pre, c_pre, beta=NULL, MaxIter=1000, theta=0.01, mu=0, prior = "LB",
               sigma=1, lambda=1, em_on=TRUE, c_trick=FALSE, pic_on=TRUE, tor=1e-10, mes = 1)   {
  
  ## Data pre-processing
  N <- dim(X_pre)[2]
  M <- dim(X_pre)[1]
  flag <- is.null(beta)
  dataAfter <- Preprocessing(X_pre, y_pre, c_pre, c_trick)
  X <- dataAfter$X
  S <- dataAfter$S
  V <- dataAfter$V
  y <- dataAfter$y
  c <- dataAfter$c
  eps <- 0.00000000001

  ## Message initialization
  v1_plus <- 1
  m1_plus <- rep(0, M)
  v1_plus_old <- v1_plus
  m1_plus_old <- m1_plus
  v1_sub_old <- v1_plus
  m1_sub_old <- m1_plus
  v0_plus <- 1
  m0_plus <- rep(0, N)
  v0_plus_old <- v0_plus
  m0_plus_old <- m0_plus
  v0_sub_old <- v0_plus
  m0_sub_old <- m0_plus
  
  ## Evaluation initialization
  MSE_error <- rep(0, MaxIter)
  score <-  Inf*rep(1, MaxIter)
  ci <-  rep(0, MaxIter)
  beta_est <- matrix(rep(0,N*MaxIter), N, MaxIter)
  
  iter_times <- 0
  ## Main Loop
  for(i in 1:MaxIter){
    
    ## step one
    
    z_sub <- Cox_Proj(m1_plus, v1_plus, c, M)
    hatz_sub <- z_sub$m
    vz_sub <- z_sub$v
    
    vz_sub_old <- vz_sub
    v1_sub <- vz_sub/(1-vz_sub/v1_plus)
    m1_sub <- v1_sub*(hatz_sub/vz_sub-m1_plus/v1_plus)
    
    #damping
    mv <- damping(v1_sub,v1_sub_old, mes)
    v1_sub <- mv$x
    v1_sub_old <- mv$x_old
    mm <- damping(m1_sub,m1_sub_old, mes)
    m1_sub <- mm$x
    m1_sub_old <- mm$x_old
    
    ## step two
    
    Qx_sub <- V %*% diag(1/(S/v1_sub + 1/v0_plus)) %*% t(V)
    hatx_sub <- Qx_sub %*% (t(X) %*% m1_sub /v1_sub + m0_plus/v0_plus)
    vx_sub <- mean(diag(Qx_sub))
    
    v0_sub <- vx_sub/(1-vx_sub/v0_plus)
    m0_sub <- v0_sub*(hatx_sub/vx_sub - m0_plus/v0_plus)
    
    #damping
    mv <- damping(v0_sub,v0_sub_old, mes)
    v0_sub <- mv$x
    v0_sub_old <- mv$x_old
    mm <- damping(m0_sub,m0_sub_old, mes)
    m0_sub <- mm$x
    m0_sub_old <- mm$x_old
    
    ## step three
    if (prior == "GB") {
      x_plus <- GB_Proj(theta,mu,sigma,m0_sub,v0_sub)
    } else if(prior == "LB") {
      x_plus <- LB_Proj(theta,lambda,m0_sub,v0_sub)
    }
    hatx_plus <- x_plus$m
    vx_plus <- x_plus$v
    
    beta_est[,i] <- hatx_plus
    vx_plus_old <- vx_plus
    v0_plus <- 1/(1/vx_plus - 1/v0_sub)
    m0_plus <- v0_plus*(hatx_plus/vx_plus - m0_sub/v0_sub)
    
    #旧值替换
    if (v0_plus<eps) {
      v0_plus <- v0_plus_old
      m0_plus <- m0_plus_old
    }
    
    #damping
    mv <- damping(v0_plus,v0_plus_old, mes)
    v0_plus <- mv$x
    v0_plus_old <- mv$x_old
    mm <- damping(m0_plus,m0_plus_old, mes)
    m0_plus <- mm$x
    m0_plus_old <- mm$x_old
    
    ## step four
    Qx_plus <- V %*% diag(1/(S/v1_sub + 1/v0_plus)) %*% t(V)
    mx_plus <- Qx_plus %*% (t(X) %*% m1_sub /v1_sub + m0_plus/v0_plus)
    hatz_plus <- X %*% mx_plus
    vz_plus <- mean(diag(X %*% Qx_plus %*% t(X)))
    
    v1_plus <- 1/(1/vz_plus-1/v1_sub)
    m1_plus <- v1_plus*(hatz_plus/vz_plus-m1_sub/v1_sub)
    
    #damping
    mv <- damping(v1_plus,v1_plus_old, mes)
    v1_plus <- mv$x
    v1_plus_old <- mv$x_old
    mm <- damping(m1_plus,m1_plus_old, mes)
    m1_plus <- mm$x
    m1_plus_old <- mm$x_old
    
    ## EM
    if(em_on){
      if(prior == "LB") {
        theta <- EM_LB(theta,lambda,m0_sub,v0_sub)
      } else if(prior == "GB") {
        par <- EM_GB(theta,mu,sigma,m0_sub,v0_sub,N)
        theta <- par$theta
        mu <- par$mu
        sigma <- par$sigma
      }
    }
    
    
    ## Evaluation
    y_pred <-  exp(-X%*%hatx_plus)
    ci[i] <-  CI(y_pred,c,M)
    if (!flag) {
      MSE <- 10*log10(sum((hatx_plus-beta)^2)/sum(beta^2))
      MSE_error[i] <- MSE
    }
    score[i]  <- sum(abs(y_pred-y)/y_pred);
    
    ## Stop criteria
    if(i>1){
      if(score[i]>=score[i-1]){
          hatx_plus <- hatx_plus_old
        }
      if(   (sum((hatx_plus_old-hatx_plus)^2)/sum(hatx_plus^2)) < tor   ) {
        break
      }
    }
    hatx_plus_old <- hatx_plus
    iter_times <- i
  }
  
  ## Output
  index <- score[1:i]==min(score)
  pos <- which(index)
  pos <- iter_times
  f_ci <- ci[pos]
  if (!flag) {
    MSE_error <- MSE_error[1:pos]
  }
  ci <- ci[1:pos]
  
  if (!flag) {
    f_mse <- MSE_error[index]
    Output <- list(c_index=ci, MSE=MSE_error, beta_est=beta_est[,pos], pointed_ci=f_ci,
                   iter_times = iter_times, pointed_mse=f_mse, theta_est=theta)
  }else
    Output <- list(beta_est=beta_est[,pos], pointed_ci=f_ci, theta_est=theta, iter_times = iter_times)
  
  ## Plots
  if(pic_on){
    Picd(flag, pos, MSE_error, ci, beta, beta_est[,pos])
  }
  
  return(Output)
}

