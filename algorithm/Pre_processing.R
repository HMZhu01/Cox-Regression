Preprocessing <- function(X, y, c, c_trick){
  ## Normalization
  X <- scale(X, center = TRUE)
  
  ## Sorting Data
  rankMatrix <- cbind(y,c,X)
  afterRank <- rankMatrix[order(rankMatrix[,1], decreasing=TRUE),]
  y_out <- afterRank[,1]
  if(!c_trick){
    c_out <- afterRank[,2]  
  }else{
    c_out <- rep(1,dim(X)[1]) 
  }

  X_out <- afterRank[,3:(dim(X)[2]+2)]

  ## SVD Decomposition
  out <- svd(X_out)
  S_out <- out$d^2
  V_out <- out$v
  
  ## Output
  output <- list(S=S_out, V=V_out, X=X_out, y=y_out, c=c_out)
  return(output)
}