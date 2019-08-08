small_geno_pred <- function(nmarker=10000, X, Y, ncv=50, tst_rate=0.2){
  #s <- sample(1:ncol(X), nmarker)
  #Xrandom <- X[,s]
  #X <- Xrandom
  y <- Y[,2]
  if(any(is.na(y))){
    rms <- which(is.na(y))
    y <- y[-rms]
    X <- X[-rms,]
  }
  X <- scale(X,center=TRUE,scale=TRUE)
  cvf <- matrix(0, nrow=length(y), ncol=as.numeric(ncv))
  for(j in 1:ncv){ 
    cvf[sample(1:nrow(cvf), round(nrow(cvf) * tst_rate, 0)), j] <- 1
  }
  yHat <- rep(NA, ncv)
  for(i in 1:ncv){
    cat("predicting cv-fold", i, "of", ncv)
    tst <- which(cvf[, i] == 1)
    yNA <- y
    yNA[tst] <- NA
    fm <- BGLR(y=yNA,ETA=list(list(X=X, model="BRR")), nIter=6000, burnIn=1000, saveAt='brr_', verbose=F)
    yHat[i] <- cor(fm$yHat[tst], y[tst])
  }
  accuracy <- mean(yHat)
  return(accuracy)
}
