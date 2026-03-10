### immune cluster analysis using ConsensusClusterPlus ###
xcell = readRDS("xcell_res.rds") 

library(ConsensusClusterPlus) # load package
d<-as.matrix(xcell) # convert xCell results to a matrix format suitable for clustering analysis
maxK = 6 # maximum number of clusters to try

results = ConsensusClusterPlus(d,maxK=maxK,reps=200,pItem=0.8,pFeature=1,
                               title="xcellall_CCP",clusterAlg="km",
                               innerLinkage="complete",seed=1234,plot="pdf",
                               writeTable = TRUE) # run consensus clustering analysis using the ConsensusClusterPlus function, with specified parameters for clustering algorithm, distance metric, and number of repetitions. 