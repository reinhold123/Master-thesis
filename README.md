# Master-thesis
Master thesis projects

simu_phen is a small R script designed to simulate phenotypes from real genomic data.
You can set the heritability you want the phenotype to have and the number of causative markers.

Required data are:
  - a genotype file in SNP-matrix format, where columns are SNPs (0/1) and rows are individuals
  - a phenotype file, which contains accession latitudes and IDs in Y$accession_latitude and Y$accession_id respectively
  - a kinship matrix, which you can calculate using the emma package (K <- emma.kinship(t(X)), where X is the genotype file)

As of now, the script only produces additive phenotypes. An implementation of epistatic effects is in work.
The heritability also still has a rather large variance. Still needs some work, but is good enough to distinguish between e.g.
  0.5 and 0.8 heritability.
