Picd <- function(flag,pos,MSE,c_index,beta,beta_est){
  ## Clear all plots
  graphics.off()
  if(!flag){
    ## Plot
    par(mfrow=c(2,2)) 
    layout(matrix(c(1,2,1,3),2,2))
    
    # picture 1
    plot(beta, pch=21, main = "对比", type="h", col="black", ylab = "数值", xlab = "下标")
    points(beta_est, pch=8,col="red")
    legend(x = "top", 
           legend = c("真实值", "估计值"),
           pch = c(20,8),
           bty = 'n',
           horiz = T,
           col = c("black","red"))
    
    # picture 2
    plot(MSE, main = "均方误差", ylab = "均方误差", xlab = "迭代次数", type="l")
  }
  
  # picture 3
  plot(c_index, main = "一致性指数", ylab = "一致性指数", xlab = "迭代次数", type="l")
}
