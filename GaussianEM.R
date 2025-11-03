GaussianMixtureEM <- function(data,theta){
  
  #theta should be a vector of the form {mu1,mu2,v1,v2,p1,p2} 
  #representing the means, variances, and mixture parameters 
  #of the mixture components
  
  e <- 1e-8
  #initial guesses
  mu1 <- theta[1]
  mu2 <- theta[2]
  v1 <- theta[3]
  v2 <- theta[4]
  p1 <- theta[5]
  p2 <- theta[6]
  
  Q_t <- 0
  Q_t1 <-  sum(log(p1*((p1*dnorm(data,mu1,v1))) + (p2*(p2*dnorm(data,mu2,v2))))) 
  
  while(abs(Q_t1 - Q_t) >= e){
    ##----E-step----##
    #posterior distribution of the mixture parameters 
    tau1 <- p1*dnorm(data,mu1,v1)/( p1*dnorm(data,mu1,v1) + p2*dnorm(data,mu2,v2) )
    tau2 <- p2*dnorm(data,mu2,v2)/( p1*dnorm(data,mu1,v1) + p2*dnorm(data,mu2,v2) )
    
    #Expectation of the log likelihood with current parameters 
    Q_t <- sum(log(tau1*((p1*dnorm(data,mu1,v1))) + (tau2*(p2*dnorm(data,mu2,v2))))) 
    
    ##----M-step----##
    #new parameters to maximize log likelihood
    p1<-sum(tau1)/length(data)
    p2<-sum(tau2)/length(data)
    mu1<-sum(tau1*data)/sum(tau1)
    mu2<-sum(tau2*data)/sum(tau2)
    v1<-sqrt(sum(tau1*(data-mu1)^2)/sum(tau1))
    v2<-sqrt(sum(tau2*(data-mu2)^2)/sum(tau2))
    
    #updated log likelihood
    Q_t1 <-  sum(log(tau1*((p1*dnorm(data,mu1,v1))) + (tau2*(p2*dnorm(data,mu2,v2))))) 
  }
  theta_estimate <- c(mu1,mu2,v1,v2,p1,p2)
  return(theta_estimate)
}