
####SCRIPT-10: trimming of the network, clustering of regions from connected components and connecting modules


library(igraph)

pairs=readRDS("cut_001eval/recode_edges_all_FILTER.rds")
modules=readRDS("cut_001eval/result_modules_all.rds")
colnames(pairs)=c("Query","Target","Eval","score")

pairs=pairs[as.numeric(pairs[,"Eval"])<=0.00001,]

modules=cbind(modules,"","")
colnames(modules)[(ncol(modules)-1):ncol(modules)]=c("subnet_trim","clust_trim")

####sizes of recoded regions 

length_structures=as.numeric(modules[,"end"])-as.numeric(modules[,"start"])

####stablishing filter to trim the network, length<350 and pval<=0.00001 -5

criteria=pairs[,"Query"]%in%rownames(modules)[length_structures<=350] & 
  pairs[,"Target"]%in%rownames(modules)[length_structures<=350]  

net=graph_from_data_frame(d=pairs[criteria,c("Query","Target")],directed=F)
E(net)$weight=as.numeric(pairs[criteria,"Eval"])
E(net)$score=as.numeric(pairs[criteria,"score"])
subnet=components(net)

modules[names(subnet$membership),"subnet_trim"]=subnet$membership

saveRDS(modules,"cut_001eval/result_modules_all_trim350_pval5.rds")

####Selection of connected components for walktrap clustering (at least 100 members)

sel=table(subnet$membership)
sel=sel[sel>=100]


criteria=pairs[,"Query"]%in%rownames(modules)[modules[,"subnet_trim"]%in%names(sel)] & 
  pairs[,"Target"]%in%rownames(modules)[modules[,"subnet_trim"]%in%names(sel)] 

pairs=pairs[criteria,]
rm(net)
gc()

for (i in 1:length(sel)){
  
  criteria=pairs[,"Query"]%in%rownames(modules)[modules[,"subnet_trim"]%in%names(sel)[i]] & 
    pairs[,"Target"]%in%rownames(modules)[modules[,"subnet_trim"]%in%names(sel)[i]] 
  
  
  net=graph_from_data_frame(d=pairs[criteria,c("Query","Target")],directed=F)
  E(net)$weight=as.numeric(pairs[criteria,"Eval"])
  E(net)$score=as.numeric(pairs[criteria,"score"])
  
  cwt=cluster_walktrap(	net, 
                        weights = E(net)$weight, 
                        steps = 6,
                        merges = TRUE, 
                        modularity = TRUE, 
                        membership = TRUE)
  
  modules[V(net)$name,"clust_trim"]=cwt$membership
  
}

###saving the tables with connected components and clustr information

saveRDS(modules,"cut_001eval/result_modules_all_trim350_pval5.rds")

################################################################################
################################################################################
################################################################################
####Table of modules (connected components plus clusters)

modules=readRDS("cut_001eval/result_modules_all_trim350_pval5.rds")
sum(modules[,"clust_trim"]!="")
sum(modules[,"subnet_trim"]!="")

modules[modules[,"clust_trim"]!="","subnet_trim"]=paste(modules[modules[,"clust_trim"]!="","subnet_trim"],
                                                        modules[modules[,"clust_trim"]!="","clust_trim"],
                                                        sep=";")

###Protein regions that are part of a community (connected component or cluster)

modules=modules[modules[,"subnet_trim"]!="",]

mod_pfam_long=readRDS("cut_001eval/result_longCut_all.rds")
mod_pfam_long=mod_pfam_long[mod_pfam_long[,"ID"]%in%rownames(modules),]
mod_pfam_long[,"subnet"]=modules[mod_pfam_long[,"ID"],"subnet_trim"]

#####Table of modules 

pfam_clusters=cbind(table(modules[,"subnet_trim"]),
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "")

colnames(pfam_clusters)=c("counts_chunks",
                          "counts_chunks_pfam",
                          "counts_chunks_pfam_075",
                          "pfam_first",
                          "pfam_first_count",
                          "pfam_first_rel",
                          "pfam_second",
                          "pfam_second_count",
                          "pfam_second_rel")

pfam_clusters=pfam_clusters[as.numeric(pfam_clusters[,1])>=5,]

####Annotating the modules with pfam from regions

modules=modules[modules[,"subnet_trim"]%in%rownames(pfam_clusters),]
mod_pfam_long=mod_pfam_long[mod_pfam_long[,"subnet"]%in%rownames(pfam_clusters),]

temp=table(modules[modules[,"pfam_counts"]!="","subnet_trim"])
temp=temp[names(temp)%in%rownames(pfam_clusters)]
pfam_clusters[names(temp),"counts_chunks_pfam"]=temp

###intermediate table

intermedio=mod_pfam_long[,c("subnet","pfam_names","ID")]
intermedio=intermedio[!duplicated(intermedio),]

for (i in 1:nrow(pfam_clusters)){
  
  temp=cbind(intermedio[intermedio[,"subnet"]%in%rownames(pfam_clusters)[i],"pfam_names"],
             intermedio[intermedio[,"subnet"]%in%rownames(pfam_clusters)[i],"ID"])
  
  if(length(temp)>2){
    
    temp1=table(temp[,1])
    temp1=temp1[order(temp1,decreasing=T)]
    
    if(length(temp1)>1){
      
      pfam_clusters[i,c("pfam_first","pfam_second")]=names(temp1)[1:2]
      pfam_clusters[i,c("pfam_first_count","pfam_second_count")]=temp1[1:2]
      pfam_clusters[i,c("counts_chunks_pfam_075")]=sum(!duplicated(temp[,2]))
      
    }else{
      
      pfam_clusters[i,c("pfam_first","pfam_second")]=c(names(temp1),"")
      pfam_clusters[i,c("pfam_first_count","pfam_second_count")]=c(temp1,"")
      pfam_clusters[i,c("counts_chunks_pfam_075")]=sum(!duplicated(temp[,2]))
    }
    
  }else{
    
    pfam_clusters[i,c("pfam_first","pfam_second")]=c(temp[1],"")
    pfam_clusters[i,c("pfam_first_count","pfam_second_count")]=c(1,"")
    pfam_clusters[i,c("counts_chunks_pfam_075")]=1
    
  }
  
}

####REcoding for missing

sum(is.na(pfam_clusters[,"pfam_first"]))
sum(is.na(pfam_clusters[,"pfam_second"]))

pfam_clusters[is.na(pfam_clusters[,"pfam_first"]),"counts_chunks_pfam_075"]="0"
pfam_clusters[is.na(pfam_clusters[,"pfam_first"]),"pfam_first_count"]="0"
pfam_clusters[is.na(pfam_clusters[,"pfam_first"]),"pfam_first"]=""
pfam_clusters[pfam_clusters[,"pfam_second_count"]=="","pfam_second_count"]="0"
pfam_clusters[pfam_clusters[,"pfam_first_count"]=="","pfam_first_count"]="0"
pfam_clusters[pfam_clusters[,"counts_chunks_pfam_075"]=="0","pfam_first_count"]="0"
pfam_clusters[pfam_clusters[,"counts_chunks_pfam_075"]=="0","pfam_second_count"]="0"

pfam_clusters[pfam_clusters[,"counts_chunks_pfam_075"]!="0","pfam_first_rel"]=as.numeric(pfam_clusters[pfam_clusters[,"counts_chunks_pfam_075"]!="0","pfam_first_count"])/
  as.numeric(pfam_clusters[pfam_clusters[,"counts_chunks_pfam_075"]!="0","counts_chunks_pfam_075"])

pfam_clusters[pfam_clusters[,"counts_chunks_pfam_075"]!="0","pfam_second_rel"]=as.numeric(pfam_clusters[pfam_clusters[,"counts_chunks_pfam_075"]!="0","pfam_second_count"])/
  as.numeric(pfam_clusters[pfam_clusters[,"counts_chunks_pfam_075"]!="0","counts_chunks_pfam_075"])

saveRDS(pfam_clusters,"cut_001eval/communities_trim350_pval5_min5struc.rds")###supplementary table


################################################################################
################################################################################
#####Selection of edges connecting modules 

pairs=readRDS("cut_001eval/recode_edges_all_FILTER.rds")
colnames(pairs)=c("Query","Target","Eval","score")

modules=readRDS("cut_001eval/result_modules_all_trim350_pval5.rds")

modules[modules[,"clust_trim"]!="","subnet_trim"]=paste(modules[modules[,"clust_trim"]!="","subnet_trim"],
                                                        modules[modules[,"clust_trim"]!="","clust_trim"],
                                                        sep=";")

###Filtering

modules=modules[modules[,"subnet_trim"]!="",]

temp=table(modules[,"subnet_trim"])
temp=temp[temp>=10]

modules=modules[modules[,"subnet_trim"]%in%names(temp),]

####Finding the edges connecting modules

recoding=modules[,"subnet_trim"]
names(recoding)=rownames(modules)

pairs=pairs[pairs[,"Query"]%in%names(recoding) & 
              pairs[,"Target"]%in%names(recoding)  ,]

pairs[,"Query"]=recoding[pairs[,"Query"]]
pairs[,"Target"]=recoding[pairs[,"Target"]]

###Keeping the ones with smaller evalues

net=graph_from_data_frame(d=pairs[,c("Query","Target")],directed=F)
E(net)$weight=as.numeric(pairs[,"Eval"])
E(net)$score=as.numeric(pairs[,"score"])

net=igraph::simplify(net,remove.loops = T,remove.multiple = T , edge.attr.comb = c(weight="min",score="max","ignore"))

saveRDS(cbind(get.edgelist(net,names = T),
              E(net)$weight,
              E(net)$score),
        "cut_001eval/recode_edges_all_comunities_trim350_pval5.rds")

################################################################################
################################################################################
################################################################################


















