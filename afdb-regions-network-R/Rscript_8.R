


####Script8: PFAM annotation

###Selecting pfam mapped to our sequences from the database 

genes=readRDS("cut_001eval/recode_nodes_all_FILTER_components5.rds")
genes=names(genes)
genes=unlist(strsplit(genes,"_"))[c(T,F)]
length(genes)

raw_data_path <- "Pfam-A.regions.uniprot.tsv.gz"
raw_data_nrow <- length(count.fields(raw_data_path))

###open connexion
con = file(raw_data_path, "r")
###Getting ready for filtering 
filtered_data=read.delim(con, nrows=raw_data_nrow-2.72e+08)
filtered_data=filtered_data[filtered_data[,1]%in%genes ,c("uniprot_acc","pfamA_acc","seq_start","seq_end")]
rownames(filtered_data)=NULL
gc()

for (i in 1:136){
  
  temp=read.delim(con, nrows=2000000,header=F)
  temp=temp[temp[,1]%in%genes ,c(1,5,6,7)]
  
  colnames(temp)= colnames(filtered_data)
  
  filtered_data=rbind(filtered_data,
                      temp)
  gc()
  print(i)
  
}
close(con)

saveRDS(filtered_data,
        "cut_eval001/filtered_pfam_cut001.rds")

##########################################################################################

####Splitting data to run annotation in paralel

###All positions
nodes=readRDS("cut_001eval/recode_nodes_all.rds")
###All components
modules=readRDS("cut_001eval/recode_nodes_all_FILTER_components5.rds")
modules=cbind(modules,
              unlist(strsplit(names(modules),"_"))[c(T,F)])
modules=cbind(modules,"","")
colnames(modules)=c("subnet","uniprot","from","to")
###pfam domains
pfam=as.matrix(readRDS("cut_001eval/filtered_pfam_cut001.rds"))
pfam=pfam[pfam[,"uniprot_acc"]%in%modules[,"uniprot"],]

####coordinates to cut in 1000 files 
coor=cbind(seq(from=1, to=3725915, by=3726),
           c(seq(from=3726, to=3725915, by=3726),3725915))

path_node="cut_001eval/pfam_anot/recode_nodes_"
path_module="cut_001eval/pfam_anot/modules_"
path_pfam="cut_001eval/pfam_anot/pfam_"

for (i in 1:nrow(coor)){
  
  temp_module=modules[coor[i,1]:coor[i,2],]
  saveRDS(nodes[nodes%in%rownames(temp_module)],paste(path_node,i,sep=""))
  saveRDS(pfam[pfam[,"uniprot_acc"]%in%temp_module[,"uniprot"],],paste(path_pfam,i,sep=""))
  saveRDS(temp_module,paste(path_module,i,sep=""))
}









