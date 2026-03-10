### immune cell infiltration analysis using xCell ###
library(xCell) # load xCell package

tpm = readRDS("rna_matrix.rds") # load TPM matrix, rows are genes, columns are samples
res<-xCellAnalysis(tpm) # run xCell analysis, the output is a matrix of cell type enrichment scores, rows are cell types, columns are samples

saveRDS(res,"xcell_res.rds") # save the xCell results as an RDS file for later use