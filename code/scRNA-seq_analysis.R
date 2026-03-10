### scRNA-seq analysis ###
library(Seurat) # load Seurat package for single-cell RNA-seq analysis
load("data_out/3_RNA_filtered.RData")  
RNAharm <- RunHarmony(RNA,"sample",plot_convergence = F, dims.use = 1:30) 
DimPlot(RNAharm, reduction = "harmony", label = F, group.by = "sample") 
ElbowPlot(RNAharm, reduction = "harmony", ndims = 30)

RNAharm <- RunUMAP(RNAharm, reduction = "harmony", dims = 1:30)
RNAharm <- FindNeighbors(object =RNAharm, reduction = "harmony", dims = 1:30)
RNAharm <- FindClusters(object =RNAharm, resolution = 0.4)

DimPlot(RNAharm, reduction = "umap",label = T) 
DimPlot(RNAharm, reduction = "umap",group.by = "sample")
FeaturePlot(object = RNAharm, features = c("Ptprc","Cd3d","Cd3e","Cd4","Cd8a",
                                           "Prf1","Nkg7","Cd79a","Ms4a1","Gzmb",
                                           "C1qc","S100a9","Csf3r","Lamp3","Clec10a",
                                           "Kitl","Cd37","Cst3","Il6st","Col3a1","Foxp3","Ifng","Pdcd1","Ccr7"), 
										   ncol=4,cols = c("grey", "red"), reduction = "umap")
