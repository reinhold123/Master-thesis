# Master-thesis
Master thesis projects

simu_phen is a small R script designed to simulate phenotypes from real genomic data.
You can set the heritability you want the phenotype to have and the number of causative markers.
Required data are:
  - a genotype file in SNP-matrix format, where columns are SNPs (0/1) and rows are individuals
  - a phenotype file, which contains accession latitudes and IDs in Y$accession_latitude and Y$accession_id respectively, as well       as phenotype values
  - a kinship matrix, which you can calculate using the emma package (K <- emma.kinship(t(X)), where X is the genotype file)
The option "model" allows to switch between purely additive phenotypes, and phenotypes with epistatic effects. "epistatic_order" lets you set the number of SNPs which produce epistatic effects.
The heritability will never be completely accurate to the set heritability value (due to some randnomness in noise generation), but will generally be close enough to distinguish between different values.
  
small_bglr is a script for quickly running genomic prediction via GBLUP using the BGLR package. Required data are X and Y (see above).

anngenopred is a Python script for building neural networks and running genomic prediction with them. It is based on the keras package (and therefore on tensorflow) and in addition requires pandas and numpy. The script consists of two functions:
  - build_network builds the neural network with specified paramters using the keras API and returns the model
  - genopred runs a genomic prediction with the created model and returns the prediction accuracies
In addition there are two pipelines in the script. One for a hyperparameter grid search, where a number of hyperparameters for the network can be set and will then be tested for their prediction accuracy. The other pipeline is designed for the final experiment of this thesis, where simply all simulated phenotypes are being predicted with set parameters (which can be taken from the grid search beforehand).
