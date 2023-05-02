

###SCRIPT2: cutting for evalue<=0.001 plus protein tables



###Argument to loop in a cluster system

i <- as.numeric(commandArgs(trailingOnly = TRUE))

path1="interactome_edges_temp/edges_" ###path with the 561 original file partitions

path2_edges="cut_001eval/interactome_edges_temp/edges_"###new location of filtered tables
path2_nodes="cut_001eval/interactome_nodes_temp/nodes_"
path2_nodes_nonred="cut_001eval/interactome_nodes_temp/nodes_prot"

library(igraph)

setwd("/hps/nobackup/beltrao/ibarrioh/AF-3D/071022/interactome_edges_temp")

###open original tables, filter names to keep uniprot IDs
pairs=readRDS(paste(path1,i,".rds",sep=""))
colnames(pairs)=c("Query","Target","Eval","score")
pairs=pairs[as.numeric(pairs[,"Eval"])<=0.001,]
pairs[,"Query"]=gsub("-F1-model_v3.cif","",pairs[,"Query"])
pairs[,"Query"]=gsub("AF-","",pairs[,"Query"])
pairs[,"Target"]=gsub("-F1-model_v3.cif","",pairs[,"Target"])
pairs[,"Target"]=gsub("AF-","",pairs[,"Target"])
gc()

###Use igraph to remove duplicities keeping max score and min evalue
net=graph_from_data_frame(d=pairs[,c("Query","Target")],directed=F)
E(net)$weight=as.numeric(pairs[,"Eval"])
E(net)$score=as.numeric(pairs[,"score"])
net=igraph::simplify(net,remove.loops = T,remove.multiple = T , edge.attr.comb = c(weight="min",score="max","ignore"))
Isolated = which(igraph::degree(net)==0)
length(Isolated)
net = delete.vertices(net, Isolated)

nodes=V(net)$name

###Saving protein sequence regions 
saveRDS(nodes,paste(path2_nodes,i,".rds",sep=""))
saveRDS(cbind(get.edgelist(net,names = T),
              E(net)$weight,
              E(net)$score),paste(path2_edges,i,".rds",sep=""))

###saving unique roteins the regions beong to
nodes=unlist(strsplit(nodes,"_"))[c(T,F,F)]
nodes=nodes[!duplicated(nodes)]
saveRDS(nodes,paste(path2_nodes_nonred,i,".rds",sep=""))
