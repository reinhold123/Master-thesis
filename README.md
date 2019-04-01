# Master-thesis
Master thesis projects

simu_phen is a small R script designed to simulate phenotypes from real genomic data.
You can set the heritability you want the phenotype to have and the number of causative markers.

Required data are:
  - a genotype file in SNP-matrix format, where columns are SNPs (0/1) and rows are individuals
  - a phenotype file, which contains accession latitudes and IDs in Y$accession_latitude and Y$accession_id respectively
  - a kinship matrix, which you can calculate using the emma package (K <- emma.kinship(t(X)), where X is the genotype file)

The option "model" allows to switch between purely additive phenotypes, and phenotypes with epistatic effects. "epistatic_order" lets you set the number of SNPs which produce epistatic effects and "epistatic_variance_explained" sets the proportion of the epistatic effect in relation to the additive effect.
The heritability still has a rather large variance. Still needs some work, but is good enough to distinguish between e.g.
  0.5 and 0.8 heritability.
  
small_bglr is a script for quickly running genomic prediction using the BGLR package.
