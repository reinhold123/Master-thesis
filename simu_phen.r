###################### R script for simulating additive phenotypes from real genomic data of A. thaliana

simu_phen <- function(X, Y, K, herit=0.5, num_additive_markers=5, epistatic_order=2, epistatic_variance_explained=0.5, report=TRUE, model="additive"){
  if(model != "additive" & model != "epistatic" & model != "mixed"){
    cat("Mission aborted. Model has to be either mixed, additive or epistatic.\n")
    return(NULL)
  }
  X <- X
  Y_temp <- Y
  K <- K
  
  #latitude_value <- 0.5 + (((Y_temp$accession_latitude - min(na.omit(Y_temp$accession_latitude)))*0.5)/
  #                           (max(na.omit(Y_temp$accession_latitude)) - min(na.omit(Y_temp$accession_latitude))))
  latitude_value <- Y_temp$accession_latitude
  lat <- na.omit(latitude_value)
  for (k in 1:length(latitude_value)){
    if ((k %% 2) == 0){
      latitude_value[k] <- mean(lat)
    }
  }
  rm(lat)
  ecotype <- as.numeric(Y_temp$accession_id)
  rownames(X)<-ecotype
  
  herit <- herit
  num_causal_predictors <- num_additive_markers
  veb <- herit/num_causal_predictors
  effective_markers <- sample(ncol(X), size = num_causal_predictors)
  count <- 0
  prel_phenotype_value <- latitude_value
  
  if(model == "additive" | model == "mixed"){
    while(count < 5){
      genomic_value <- rep(0, length(rownames(X)))
      for (i in effective_markers){
        beta <- sqrt((veb*var(na.omit(prel_phenotype_value))) / ((var(X[,i])) - (veb * var(X[,i]))))
        for (j in 1:length(genomic_value)){
          genomic_value[j] <- genomic_value[j] + beta * X[j,i]}
      }
      prel_phenotype_value <- latitude_value + genomic_value
  
      noise_value <- var(genomic_value) * (1-herit)/herit
      prel_phenotype_value <- prel_phenotype_value + rnorm(length(prel_phenotype_value), sd=sqrt(noise_value))
      count <- count + 1
    }
  }
  
  if(model == "epistatic" | model == "mixed"){
    cat("Creating epistatic markers... this might take a while...\n")
    andmarkers <- as.data.frame(matrix(nrow=nrow(Y_temp), ncol=1000))
    ormarkers <- as.data.frame(matrix(nrow=nrow(Y_temp), ncol=1000))
    
    for (j in 1:1000){
      repeat {
        epistatic_markers <- sample(ncol(X), size=2)
        for (i in c(1,2)){
          assign(paste0("bar", i), sum(X[, epistatic_markers[i]])/nrow(X))
          if (get(paste0("bar", i)) > 1-get(paste0("bar", i))){
            assign(paste0("MAF", i), 1-get(paste0("bar", i)))
          }
          else{
            assign(paste0("MAF", i), get(paste0("bar", i)))
          }
        }
        if (MAF1 > 0.3 & MAF1 < 0.5 & MAF2 > 0.3 & MAF2 < 0.5) break
      }
      
      foo <- as.data.frame(matrix(nrow=nrow(Y_temp), ncol=4))
      colnames(foo) <- c("marker1", "marker2", "andmarker", "ormarker")
      foo$marker1 <- X[, epistatic_markers[1]]
      foo$marker2 <- X[, epistatic_markers[2]]
      
      for (i in 1:length(foo$marker1)){
        if (foo$marker1[i] == 1 & foo$marker2[i] == 1){
          foo$andmarker[i] <- 1
        }
        else{
          foo$andmarker[i] <- 0
        }
        #if (foo$marker1[i] == 1 | foo$marker2[i] == 1){
         # foo$ormarker[i] <- 1
        #}
        #else{
         # foo$ormarker[i] <- 0
        #}
      }
      andmarkers[, j] <- foo$andmarker
      #ormarkers [, j] <- foo$ormarker
    }
    
    #snp_mean <- apply(as.matrix(ormarkers), 2, mean)
    #snp_sd <- apply(as.matrix(ormarkers), 2, sd)
    #snpm <- matrix(nrow=nrow(ormarkers), ncol=ncol(ormarkers), data=snp_mean, byrow=T)
    #snpd <- matrix(nrow=nrow(ormarkers), ncol=ncol(ormarkers), data=snp_sd, byrow=T)
    #xx_gut <- (as.matrix(ormarkers)-snpm)/snpd
    
    snp_mean <- apply(X,2,mean)
    snp_sd <- apply(X,2,sd)
    snpm <- matrix(nrow=nrow(X), ncol=ncol(X), data=snp_mean, byrow=T)
    snpd <- matrix(nrow=nrow(X), ncol=ncol(X), data=snp_sd, byrow=T)
    xx_stand <- (X-snpm)/snpd
    
    #or_r2 <- (crossprod(xx_gut, xx_stand)/nrow(X))^2
    
    snp_mean <- apply(as.matrix(andmarkers), 2, mean)
    snp_sd <- apply(as.matrix(andmarkers), 2, sd)
    snpm <- matrix(nrow=nrow(andmarkers), ncol=ncol(andmarkers), data=snp_mean, byrow=T)
    snpd <- matrix(nrow=nrow(andmarkers), ncol=ncol(andmarkers), data=snp_sd, byrow=T)
    xx_gut <- (as.matrix(andmarkers)-snpm)/snpd
    
    and_r2 <- (crossprod(xx_gut, xx_stand)/nrow(X))^2
    
    rm(snp_mean, snp_sd, snpm, snpd, xx_stand, xx_gut)
    
    and_rms <- NULL
    #or_rms <- NULL
    
    for (k in 1:nrow(and_r2)){
      if (any(and_r2[k,] > 0.5)){
        and_rms[k] <- k
      }
      #if (any(or_r2[k,] > 0.5)){
      #  or_rms[k] <- k
      #}
    }
    
    andmarkers <- andmarkers[, -na.omit(and_rms)]
    #ormarkers <- ormarkers[, -na.omit(or_rms)]
    
    foo <- sample(ncol(andmarkers), size=5)
    #bar <- sample(ncol(ormarkers), size=5)
    epistatic_value <- andmarkers[,foo[1]] + andmarkers[,foo[2]] + andmarkers[,foo[3]] + andmarkers[,foo[4]] + andmarkers[,foo[5]] 
    #+ ormarkers[,bar[1]] + ormarkers[,bar[2]] + ormarkers[,bar[3]] + ormarkers[,bar[4]] + ormarkers[,bar[5]]
    
    if(model == "epistatic"){
      latitude_value <- 0.5 + (((latitude_value - min(na.omit(latitude_value)))*0.5/
                                  max(na.omit(latitude_value) - min(na.omit(latitude_value)))))
      prel_phenotype_value <- latitude_value + epistatic_value
      noise_value <- var(epistatic_value) * (1-herit)/herit
      prel_phenotype_value <- prel_phenotype_value + rnorm(length(prel_phenotype_value), sd=sqrt(noise_value))
    }
  }
  
  if(model == "mixed"){
    cat("Mixed model not yet implemented. Choose either additive or epistatic model.\n")
    return(NULL)
  }
  
  phenotype_value <- prel_phenotype_value
  Y <- data.frame(ecotype, phenotype_value)
  if(report == TRUE){
    cat("Brace yourself, results are coming...\n")
    h2test <- amm_gwas(Y, X, K, run=FALSE, report=FALSE)
    cat("Set heritability was", herit, "\t estimated heritability is", h2test[1], "\n")
    Y <- list(Y, h2test[1])
  }
  return(Y)
}
