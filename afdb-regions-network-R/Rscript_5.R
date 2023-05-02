

####SCRIPT5: assembling clustering results to one file

setwd("cut_001eval/hierar")

paths_clust=files=list.files(pattern="cluster")
paths_start=files=list.files(pattern="start")
paths_stop=files=list.files(pattern="stop")
paths_uniprot=files=list.files(pattern="uniprot")

paths_uniprot=paths_uniprot[!paths_uniprot%in%"uniprotID_0"]
paths_start=paths_start[!paths_start%in%"start_0"]
paths_stop=paths_stop[!paths_stop%in%"stop_0"]

####Empezamos  

temp=readRDS("uniprotID_0")

names_vect=paste(gsub("_1","",temp),
                 readRDS("start_0"),
                 readRDS("stop_0"),sep="_")

final_vect=temp
names(final_vect)=names_vect

for (i in 1:length(paths_uniprot)){
  
  temp_uni=readRDS(paths_uniprot[i])
  temp_from=readRDS(paths_start[i])
  temp_to=readRDS(paths_stop[i])
  temp_cl=readRDS(paths_clust[i])
  
  temp=paste(temp_uni,temp_cl,sep="_")
  names(temp)=paste(temp_uni,temp_from,temp_to,sep="_")
  
  final_vect=c(final_vect,temp)
  
}

saveRDS(final_vect,"cut_001eval/mapping_chunks.rds")










