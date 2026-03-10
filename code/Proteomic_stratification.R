### Proteomic_stratification based on NMF ###
library(NMF)
d=readRDS("protein_matrix_95sd.rds")
rank <- 3
seed <- 20220321
mut.nmf <- nmf(d, 
               rank = rank, 
               seed = seed, 
               nrun= 200,
               method = "brunet")