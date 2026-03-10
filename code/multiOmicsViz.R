### multiOmicsViz ###
# load multiOmicsViz package 
library("multiOmicsViz")

# load multi-omics data
mut.df = readRDS('lcnec-data/mut.df.rds')
pro.df = readRDS('lcnec-data/pro.df.rds')
phs.df = readRDS('lcnec-data/phs.df.rds')
mrna.df = readRDS('lcnec-data/mrna.df.rds')
cnv.df = readRDS('lcnec-data/cnv.df.rds')

## multi-omics data visualization using multiOmicsViz package
sourceOmics<-cnv.df
targetOmicsList <- list()
targetOmicsList[[1]] <- mrna.df[,intersect(colnames(pro.df),colnames(mrna.df))]
targetOmicsList[[2]] <- as.data.frame(pro.df[,intersect(colnames(pro.df),colnames(mrna.df))])

outputfile <- "half_mrna-half_protein-0.05-2"

multiOmicsViz(sourceOmics,"CNA","All",targetOmicsList,c("half_mRNA","half_protein"),"All",0.05,outputfile,nThreads =2) 

phs.df  = data.frame(phs.df)
phs.df$gene = strsplit2(rownames(phs.df),'_')[,1]
phs.df = phs.df[!duplicated(phs.df$gene),]
rownames(phs.df) = phs.df$gene

sourceOmics<-cnv.df
targetOmicsList <- list()
targetOmicsList[[2]] <- mrna.df
targetOmicsList[[1]] <- phs.df

outputfile <- "half_phs-half_mrna-0.05-2"

multiOmicsViz(sourceOmics,"CNA","All",targetOmicsList,c("half_phs","half_mRNA"),"All",0.05,outputfile,nThreads =2) 
