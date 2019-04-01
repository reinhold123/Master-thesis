small_geno_pred <- function(nmarker=10000, X, Y, A){
  s <- sample(1:ncol(X), nmarker)
  Xrandom <- X[,s]
  X <- Xrandom
  y <- Y[,2]
  if(any(is.na(y))){
    rms <- which(is.na(y))
    y <- y[-rms]
    X <- X[-rms,]
  }
  X <- scale(X,center=TRUE,scale=TRUE)
  G <- tcrossprod(X) / ncol(X)
  EVD <- eigen(G)
  ETA <- list(list(K=A, model='RKHS'), list(V=EVD$vectors,d=EVD$values, model='RKHS'))
  fmBRR <- BGLR(y=y,ETA=list(list(X=X,model='BRR')), nIter=6000, burnIn=1000,saveAt='brr_', verbose=F)
  accuracy <- cor(fmBRR$y, fmBRR$yHat)
  return(accuracy)
}
