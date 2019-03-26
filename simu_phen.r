###################### R script for simulating additive phenotypes from real genomic data of A. thaliana

simu_phen <- function(X, Y, K, herit=0.5, num_causal_predictors=5, report=TRUE){
  X <- X
  Y_temp <- Y
  K <- K
  
  latitude_value <- 0.5 + (((Y_temp$accession_latitude - min(na.omit(Y_temp$accession_latitude)))*0.5)/
                             (max(na.omit(Y_temp$accession_latitude)) - min(na.omit(Y_temp$accession_latitude))))
  ecotype <- as.numeric(Y_temp$accession_id)
  rownames(X)<-ecotype
  
  herit <- herit
  num_causal_predictors <- num_causal_predictors
  veb <- herit/num_causal_predictors
  effective_markers<- sample(ncol(X), size = num_causal_predictors)
  count <- 0
  prel_phenotype_value <- latitude_value
  
  while(count < 3){
    genomic_value <- rep(0, length(rownames(X)))
    for (i in effective_markers){
      beta <- sqrt((veb*var(na.omit(prel_phenotype_value))) / ((var(X[,i])) - (veb * var(X[,i]))))
      for (j in 1:length(genomic_value)){
        genomic_value[j] <- genomic_value[j] + beta * X[j,i]}
    }
    prel_phenotype_value <- latitude_value * genomic_value
    
    noise_value <- var(genomic_value) * (1-herit)/herit
    prel_phenotype_value <- prel_phenotype_value + rnorm(length(prel_phenotype_value), sd=sqrt(noise_value))
    count <- count + 1
  }
  
  phenotype_value <- prel_phenotype_value
  Y <- data.frame(ecotype, phenotype_value)
  if(report == TRUE){
    cat("Brace yourself, results are coming...\n")
    h2test <- amm_gwas(Y, X, K, run=FALSE, report=FALSE)
    cat("Set heritability is", herit, "\t estimated heritability is", h2test[1])
  }
  return(Y)
}