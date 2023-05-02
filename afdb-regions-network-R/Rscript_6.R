




#####SCRIPT6: recoding the edges (similarity based on clustering results)

####Script meant to run inside a loop

library("igraph")

i <- as.numeric(commandArgs(trailingOnly = TRUE))

path=paste("cut_001eval/interactome_edges_temp/edges_",i,".rds",sep="")

####loading previous files
pairs=readRDS(path)
colnames(pairs)=c("Query","Target","Eval","score")
net=graph_from_data_frame(d=pairs[,c("Query","Target")],directed=F)

#loading mapping
nodes=readRDS("cut_001eval/mapping_chunks.rds")
nodes=nodes[names(nodes)%in%V(net)$name]
gc()

pairs=pairs[pairs[,"Query"]%in%names(nodes) & 
              pairs[,"Target"]%in%names(nodes),]

pairs[,"Query"]=nodes[pairs[,"Query"]]
pairs[,"Target"]=nodes[pairs[,"Target"]]

###recoding the network and eliminating duplicates

net=graph_from_data_frame(d=pairs[,c("Query","Target")],directed=F)
E(net)$weight=as.numeric(pairs[,"Eval"])
E(net)$score=as.numeric(pairs[,"score"])
net=igraph::simplify(net,remove.loops = T,remove.multiple = T , edge.attr.comb = c(weight="min",score="max","ignore"))
Isolated = which(igraph::degree(net)==0)
length(Isolated)
net = delete.vertices(net, Isolated)

###New edges

path_node=paste("cut_001eval/recode/node_",i,".rds",sep="")
path_edge=paste("cut_001eval/recode/edge_",i,".rds",sep="")

saveRDS(nodes,path_node)
saveRDS(cbind(get.edgelist(net,names = T),
              E(net)$weight,
              E(net)$score),path_edge)

##############















