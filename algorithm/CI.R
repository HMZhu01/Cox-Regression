CI <- function(y_pred,c,M){
  concordant_pair <- 0;
  discordant_pair <- 0;
  
  for (i in 1:(M-1)){
    for (j in (i+1):M){
 
      if (c[i]==TRUE && c[j]==TRUE){
      
        if (y_pred[i] >= y_pred[j]){
          concordant_pair <-  concordant_pair + 1
        } else{
          discordant_pair <-  discordant_pair + 1
        }
        
      } else if (c[i]==TRUE && c[j]==FALSE){
        
          if (y_pred[i]<y_pred[j]){
            discordant_pair <-  discordant_pair + 1
          }
          
      } else if (c[i]==FALSE && c[j]==TRUE){
        
          if (y_pred[i] >= y_pred[j]){
            concordant_pair <-  concordant_pair + 1
          }
      
      }
    }
  }
  
  CI <- concordant_pair/(concordant_pair+discordant_pair)
}

