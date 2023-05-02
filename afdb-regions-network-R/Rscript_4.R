####################################################

###SCRIPT 4: hierarchical clustering of regions per protein

###This  script is ment to run in side a loop

i <- as.numeric(commandArgs(trailingOnly = TRUE))

pathID="cut_001eval/hierar/uniprotID_"
pathstart="cut_001eval/hierar/start_"
pathstop="cut_001eval/hierar/stop_"
pathcluster="cut_001eval/hierar/cluster_"

####Protein ID plus start stop positions
####
vector_guia_uniprot=readRDS(paste(pathID,i,sep=""))
vector_guia_from=readRDS(paste(pathstart,i,sep=""))
vector_guia_to=readRDS(paste(pathstop,i,sep=""))
vector_guia_cluster=rep(1,length(vector_guia_to))
protes=table(vector_guia_uniprot)


for (j in 1:length(protes)){
  
  temp=cbind(vector_guia_from[vector_guia_uniprot%in%names(protes)[j]],
             vector_guia_to[vector_guia_uniprot%in%names(protes)[j]])
  
  vector_guia_cluster[vector_guia_uniprot%in%names(protes)[j]]=cutree(hclust(dist(temp,method = "euclidean")),h= 250)
  
}

saveRDS(vector_guia_cluster,paste(pathcluster,i))


