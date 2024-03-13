## User should first using function setwd() to set the working directory!
setwd("C:/Users/Zhangshanshu/Desktop/Cox论文相关/Rscript")
source("./algorithm/EP.R")
source("./algorithm/Sim_data.R")

## Parameters Setting
M <- 400  # observations
N <- 200  # variables
lambda <- 1
perc <- 0.1   # censoring proportion
rho <- 0.1 # sparsity rate

## Survival data generation
data <-  SurvData(M, N, lambda, snr, perc, rho, 2)
X <- data$X
y <- data$y
c <- data$c
beta <- data$beta

## Estimation
a <- Sys.time()
result <- EP(X,y,c,beta,em_on = TRUE,mes = 1)
print(Sys.time()-a)