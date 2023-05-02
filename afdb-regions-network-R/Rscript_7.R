############################################################################
############################################################################

###SCRIPT7: assembling the recoded tables, selection of connected components

#####Part1: assembly
####EDGES

path="cut_001eval/recode/edge_"


i=0

edges=readRDS(paste(path,i,".rds",sep=""))

for (i in 1:560){
  
  edges=rbind(edges,
              readRDS(paste(path,i,".rds",sep="")))
  
}

colnames(edges)=c("Query","Target","Eval","score")

###saving the assembly

saveRDS(edges,"cut_001eval/recode_edges_all.rds")

####NODES

path="cut_001eval/recode/node_"

i=0

nodes=unique(readRDS(paste(path,i,".rds",sep="")))

for (i in 1:560){
  
  nodes=c(nodes,
          readRDS(paste(path,i,".rds",sep="")))
  
}

###saving the assembly

saveRDS(nodes,"cut_001eval/recode_nodes_all.rds")

rm(edges)
rm(nodes)

###############################################################
###############################################################

###Part2: components

####We are going to selected connected components with at least 5 nodes 

library(igraph)

###################################################
####Volvemos localmente a ver que sale de todo esto

pairs=readRDS("cut_001eval/recode_edges_all.rds")

net=graph_from_data_frame(d=pairs[,c("Query","Target")],directed=F)
E(net)$weight=as.numeric(pairs[,"Eval"])
E(net)$score=as.numeric(pairs[,"score"])
net=igraph::simplify(net,remove.loops = T,remove.multiple = T , edge.attr.comb = c(weight="min",score="max","ignore"))
Isolated = which(igraph::degree(net)==0)
length(Isolated)
net = delete.vertices(net, Isolated)

##saving assmbled edges after removing dups

saveRDS(cbind(get.edgelist(net,names = T),
              E(net)$weight,
              E(net)$score),"cut_001eval/recode_edges_all_FILTER.rds")

###Components

subnet=components(net)

nodes_sel=subnet$membership[subnet$membership%in%names(table(subnet$membership)[table(subnet$membership)>=5])]

saveRDS(nodes_sel,"cut_001eval/recode_nodes_all_FILTER_components5.rds")
saveRDS(pairs[pairs[,1]%in%names(nodes_sel) & 
                pairs[,2]%in%names(nodes_sel)  ,],"cut_001eval/recode_edges_all_FILTER_components5.rds")











