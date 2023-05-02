


###SCRIPT1: partition of original table to handle size

######Due to file size, we split it in files containing 2,000,000 lines (561 files in total) 

####open connexion

raw_data_path <- "/nfs/research/beltrao/ibarrioh/AF-3D/071022/afdb50best_clu_nosingleton_repseq_ava.m8.gz"
raw_data_nrow <- length(count.fields(raw_data_path))

con = file(raw_data_path, "r")

###cutting loop

for (i in 0:560){
  
  pairs=read.delim(con, nrows=2000000,header=F)
  
  temp_node=unique(c(paste(as.character(pairs[,1]),
                           pairs[,7],
                           pairs[,8],sep="_"),
                     paste(as.character(pairs[,2]),
                           pairs[,9],
                           pairs[,10],sep="_")))
  
  temp_edge=cbind(paste(as.character(pairs[,1]),
                        pairs[,7],
                        pairs[,8],sep="_"),
                  paste(as.character(pairs[,2]),
                        pairs[,9],
                        pairs[,10],sep="_"),
                  pairs[,11],
                  pairs[,12])
  
  path_1=paste("interactome_nodes_temp/nodes_",i,".rds",sep="")
  path_2=paste("interactome_edges_temp/edges_",i,".rds",sep="")
  
  saveRDS(temp_node,path_1)
  saveRDS(temp_edge,path_2) 
  
}
close(con)


















#