#Gauss dis.
datagene_norm <- function(mu,Sigma,N){
  d <- dim(Sigma)[1]
  Z <- matrix(0,nrow=d,ncol=N)
  
  for (i in 1:N) {
    Z[,i] <- rnorm(d,0,1)
  }
  gamma<- matrix(0,d,d)
  diag(gamma) <- eigen(Sigma)$value^(1/2)
  H <- eigen(Sigma)$vector
  
  return(H%*%gamma%*%Z+mu)
}
#t dis.
datagene_t <- function(mu,Sigma,df,N){
  d <- dim(Sigma)[1]
  Sigma <- ((df-2)/df)*Sigma 
  X <- datagene_norm(rep(0,d),Sigma,N)
  Z <- datagene_norm(rep(0,df),diag(df),N) 
  xi <- apply(Z^2,2,sum)/df  
  Y <- X %*% diag(sqrt(1/xi))  
  t <- Y + mu 
  return(t)
}
#chisq dis.
datagene_chisq <- function(mu,Sigma,df,N){
  d <- dim(Sigma)[1]
  Z <-  matrix(0,d,as.integer(N))
  for (i in 1:N) {
    Z[,i] <- (rchisq(d,df) - df)/(2*df)^(1/2)
  }
  gamma<- matrix(0,d,d)
  diag(gamma) <- eigen(Sigma)$value^(1/2)
  H <- eigen(Sigma)$vector
  return(H%*%gamma%*%Z+mu)
}
