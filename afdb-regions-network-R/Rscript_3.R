

####SCRIPT3: getting ready for hierarchical clustering of regions

####Here we generate vectors with uniprot identifiers, start and stop positions for each region

setwd("cut_001eval/interactome_nodes_temp")

files=list.files(pattern = "nodes_prot")
files1=list.files(pattern = "node")
files1=files1[!files1%in%files]

nodes=readRDS(files1[1])

for (i in 2:length(files1)){
  
  nodes=c(nodes,readRDS(files1[i]))
  
}

nodes=nodes[!duplicated(nodes)]
saveRDS(nodes,"cut_001eval/all_proteins.rds")


names_unlist=unlist(strsplit(nodes,"_"))[c(T,F,F)]
saveRDS(names_unlist,"cut_001eval/all_uniprotID.rds")
prots=table(names_unlist)

print(sum(prots==1))
print(length(prots))

rm(names_unlist)
gc()

start_unlist=unlist(strsplit(nodes,"_"))[c(F,T,F)]
saveRDS(start_unlist,"cut_001eval/all_start.rds")

rm(start_unlist)
gc()

stop_unlist= unlist(strsplit(nodes,"_"))[c(F,F,T)]
saveRDS(stop_unlist,"cut_001eval/all_stop.rds")


#######################################################
#####Generating 580 part files with protein ID, start and stop positions

####Lines for proteins with only one region (no clustering)

names_unlist=readRDS("cut_001eval/all_uniprotID.rds")
start_unlist=readRDS("cut_001eval/all_start.rds")
stop_unlist=readRDS("cut_001eval/all_stop.rds")

prots=table(names_unlist)

##clean out proteins with one region
prots=prots[prots>1]
start_unlist=start_unlist[names_unlist%in%names(prots)]
stop_unlist=stop_unlist[names_unlist%in%names(prots)]
names_unlist=names_unlist[names_unlist%in%names(prots)]

### regions cut

coor=cbind(seq(from=1, to=2004717, by=3452),
           c(seq(from=3452, to=2004717, by=3452),2004717))

pathID="cut_001eval/hierar/uniprotID_"
pathstart="cut_001eval/hierar/start_"
pathstop="cut_001eval/hierar/stop_"

for (i in 1:nrow(coor)){

cluster_temp=names_unlist[names_unlist%in%names(prots)[coor[i,1]:coor[i,2]]]
start_temp=start_unlist[names_unlist%in%cluster_temp]
stop_temp=stop_unlist[names_unlist%in%cluster_temp]

saveRDS(cluster_temp,paste(pathID,i,sep=""))
saveRDS(start_temp,paste(pathstart,i,sep=""))
saveRDS(stop_temp,paste(pathstop,i,sep=""))  
  
}

####################################################
####################################################


###Creasome seed

i <- as.numeric(commandArgs(trailingOnly = TRUE))

pathID="/hps/nobackup/beltrao/ibarrioh/AF-3D/181022/cut_001eval/hierar/uniprotID_"
pathstart="/hps/nobackup/beltrao/ibarrioh/AF-3D/181022/cut_001eval/hierar/start_"
pathstop="/hps/nobackup/beltrao/ibarrioh/AF-3D/181022/cut_001eval/hierar/stop_"
pathcluster="/hps/nobackup/beltrao/ibarrioh/AF-3D/181022/cut_001eval/hierar/cluster_"
####subimos los vectores 

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






























