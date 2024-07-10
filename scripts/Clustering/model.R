# Generate a consensus clustering model with given arguments
CC_modelization <- function(datasets, title_model, cluster_alg = "hc", distance = "euclidean"){
  CC_model = ConsensusClusterPlus(datasets, maxK = 6, reps = 50, pItem = 0.8, pFeature = 1,
                                  clusterAlg = cluster_alg, distance = distance, 
                                  innerLinkage="ward.D2", finalLinkage="ward.D2",
                                  seed = 126, plot = "pdf", title = title_model)
  return(CC_model)
}