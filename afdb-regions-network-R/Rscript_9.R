

#####SCRIPT9 : annotation of regions using pfam

####Lines to annotate a region using pfam domains,, 75% of the domain must be included in the region to be considered

###Ment to run inside a loop


i <- as.numeric(commandArgs(trailingOnly = TRUE))


path_node="cut_001eval/pfam_anot/recode_nodes_"
path_module="cut_001eval/pfam_anot/modules_"
path_pfam="cut_001eval/pfam_anot/pfam_"

nodes=readRDS(paste(path_node,i,sep=""))
modules=readRDS(paste(path_module,i,sep=""))
pfam=readRDS(paste(path_pfam,i,sep=""))

for (j in 1:nrow(modules)){
  
  temp=names(nodes)[nodes==rownames(modules)[j]]
  modules[j,"from"]=min(as.numeric(unlist(strsplit(temp,"_"))[c(F,T,F)]))
  modules[j,"to"]=max(as.numeric(unlist(strsplit(temp,"_"))[c(F,F,T)]))
  
}

modules=cbind(modules,"","")
colnames(modules)=c("subnet","uniprot","start","end","pfam_names","pfam_counts")

###Ahora anotacion de pfams 

for (j in c(1:nrow(modules))[modules[,"uniprot"]%in%pfam[,"uniprot_acc"]]){
  
  criteria=pfam[,"uniprot_acc"]%in%modules[j,"uniprot"] & (
    (as.numeric(pfam[,"seq_start"])>=as.numeric( modules[j,"start"]) & 
       as.numeric(pfam[,"seq_end"])<=as.numeric( modules[j,"end"])) |
      (as.numeric(pfam[,"seq_start"])<=as.numeric( modules[j,"start"]) & 
         as.numeric(pfam[,"seq_end"])>=as.numeric( modules[j,"start"])) |
      (as.numeric(pfam[,"seq_start"])<=as.numeric( modules[j,"end"]) & 
         as.numeric(pfam[,"seq_end"])>=as.numeric( modules[j,"end"]))) 
  temp=unique(pfam[criteria,"pfamA_acc"])
  modules[j,"pfam_names"]=paste(temp,collapse=";")
  modules[j,"pfam_counts"]=length(temp)
  
}

####Ahora vamos a buscar el porcentaje

modules=cbind(modules,"")
colnames(modules)[ncol(modules)]="pfam_chunk_rel"

for (j in c(1:nrow(modules))[modules[,"pfam_counts"]=="1"]){
  
  temp=as.numeric(pfam[pfam[,"uniprot_acc"]%in%modules[j,"uniprot"] & (
    (as.numeric(pfam[,"seq_start"])>=as.numeric( modules[j,"start"]) & 
       as.numeric(pfam[,"seq_end"])<=as.numeric( modules[j,"end"])) |
      (as.numeric(pfam[,"seq_start"])<=as.numeric( modules[j,"start"]) & 
         as.numeric(pfam[,"seq_end"])>=as.numeric( modules[j,"start"])) |
      (as.numeric(pfam[,"seq_start"])<=as.numeric( modules[j,"end"]) & 
         as.numeric(pfam[,"seq_end"])>=as.numeric( modules[j,"end"])))  ,c("seq_start","seq_end")])
  
  A=c(temp[1]:temp[2])
  B=c(as.numeric(modules[j,"start"]):as.numeric(modules[j,"end"]))
  
  modules[j,"pfam_chunk_rel"]=round(sum(A%in%B)/length(A),3)
  
}


modules=cbind(modules,"")
colnames(modules)[ncol(modules)]="pfam_all"

modules[modules[,"pfam_counts"]=="1","pfam_all"]=modules[modules[,"pfam_counts"]=="1","pfam_names"]


for (j in c(1:nrow(modules))[modules[,"pfam_counts"]%in%c("2","3","4","5","6","7","8","9")]){
  
  temp=cbind(pfam[pfam[,"uniprot_acc"]%in%modules[j,"uniprot"] & (
    (as.numeric(pfam[,"seq_start"])>=as.numeric( modules[j,"start"]) & 
       as.numeric(pfam[,"seq_end"])<=as.numeric( modules[j,"end"])) |
      (as.numeric(pfam[,"seq_start"])<=as.numeric( modules[j,"start"]) & 
         as.numeric(pfam[,"seq_end"])>=as.numeric( modules[j,"start"])) |
      (as.numeric(pfam[,"seq_start"])<=as.numeric( modules[j,"end"]) & 
         as.numeric(pfam[,"seq_end"])>=as.numeric( modules[j,"end"])))  ,c("pfamA_acc","seq_start","seq_end")],"")
  
  for (k in 1:nrow(temp)){
    
    A=c(as.numeric(temp[k,"seq_start"]):as.numeric(temp[k,"seq_end"]))
    B=c(as.numeric(modules[j,"start"]):as.numeric(modules[j,"end"]))
    temp[k,4]=round(sum(A%in%B)/length(A),3)
    
  }
  
  modules[j,"pfam_chunk_rel"]=paste( temp[,4],collapse=";")
  modules[j,"pfam_all"]=paste( temp[,1],collapse=";")
  
}


################################################
################################################
############Long version table (one line per pfam per domain)


mod_pfam_long=modules[modules[,"pfam_counts"]=="1",c("subnet",
                                                     "uniprot",
                                                     "start",
                                                     "end",
                                                     "pfam_names",
                                                     "pfam_counts",
                                                     "pfam_chunk_rel")]
mod_pfam_long=cbind(rownames(mod_pfam_long),mod_pfam_long)
colnames(mod_pfam_long)[1]="ID"

temp_pfam=modules[modules[,"pfam_counts"]%in%c("2","3","4","5","6","7","8","9"),]

temp=unlist(strsplit(temp_pfam[1,"pfam_all"],";"))
primer=cbind(rep(rownames(temp_pfam)[1],length(temp)),
             rep(temp_pfam[1,"subnet"],length(temp)),
             rep(temp_pfam[1,"uniprot"],length(temp)),
             rep(temp_pfam[1,"start"],length(temp)),
             rep(temp_pfam[1,"end"],length(temp)),
             temp,
             rep(temp_pfam[1,"pfam_counts"],length(temp)),
             unlist(strsplit(temp_pfam[1,"pfam_chunk_rel"],";")))

for (j in 2:nrow(temp_pfam)){
  
  temp=unlist(strsplit(temp_pfam[j,"pfam_all"],";"))
  temp1=cbind(rep(rownames(temp_pfam)[j],length(temp)),
              rep(temp_pfam[j,"subnet"],length(temp)),
              rep(temp_pfam[j,"uniprot"],length(temp)),
              rep(temp_pfam[j,"start"],length(temp)),
              rep(temp_pfam[j,"end"],length(temp)),
              temp,
              rep(temp_pfam[j,"pfam_counts"],length(temp)),
              unlist(strsplit(temp_pfam[j,"pfam_chunk_rel"],";")))
  
  primer=rbind(primer,temp1)
  
}

mod_pfam_long=rbind(mod_pfam_long,
                    primer)
rownames(mod_pfam_long)=NULL

############################################
############################################
############################################
####Saving the files 

path_result="cut_001eval/pfam_anot/result_modules_"###Table with annotation
path_result_long="cut_001eval/pfam_anot/result_long_"###long format
path_result_long_cut="cut_001eval/pfam_anot/result_longCut_"###Long format keeping only annotation

saveRDS(modules,paste(path_result,i,sep=""))
saveRDS(mod_pfam_long,paste(path_result_long,i,sep=""))
saveRDS(mod_pfam_long[as.numeric(mod_pfam_long[,"pfam_chunk_rel"])>=0.75,],paste(path_result_long_cut,i,sep=""))



############################################
############################################
####
####Assembling all together 

path_result="cut_001eval/pfam_anot/result_modules_"
path_result_long="cut_001eval/pfam_anot/result_long_"
path_result_long_cut="cut_001eval/pfam_anot/result_longCut_"

####Primero los resultados 

table1=readRDS(paste(path_result,1,sep=""))

for (i in 2:1000){
  
  table1=rbind(table1,
               readRDS(paste(path_result,i,sep=""))) 
  
}

saveRDS(table1,"cut_001eval/result_modules_all.rds")

########

table1=readRDS(paste(path_result_long,1,sep=""))

for (i in 2:1000){
  
  table1=rbind(table1,
               readRDS(paste(path_result_long,i,sep=""))) 
  
}

saveRDS(table1,"cut_001eval/result_long_all.rds")

########

table1=readRDS(paste(path_result_long_cut,1,sep=""))

for (i in 2:1000){
  
  table1=rbind(table1,
               readRDS(paste(path_result_long_cut,i,sep=""))) 
  
}

saveRDS(table1,"cut_001eval/result_longCut_all.rds")




















