###################################
###covering analysis
###################################
library(igraph)
library(stringr)
library(lpSolve)
source("~/Documents/MultiModalProject/Cover/lpcover/R/optim.R")
source("/Users/karenxu/Documents/MultiModalProject/MultiModalPC/codes/utilities.R")
pb<-txtProgressBar(style = 3)



##########################################
########set up working files     #########
##########################################
load("/Users/karenxu/Documents/MultiModalProject/MultiModalPC/data/feature_name_list.RData")
load("/Users/karenxu/Documents/MultiModalProject/MultiModalPC/data/tumor_gm.RData")
load("/Users/karenxu/Documents/MultiModalProject/MultiModalPC/data/normal_gm.RData")
#--dataset
normal_gep<-normal[feature_name_list[[1]],]#48 samples
tumor_gep<-tumor[feature_name_list[[1]],]#94 samples
normal_met<-normal[feature_name_list[[2]],]#48 samples
tumor_met<-normal[feature_name_list[[2]],]#94 samples

###################################################
###        network from pathway commons      ######
###################################################
load("/Users/karenxu/Documents/MultiModalProject/MultiModalPC/data/S0_pc_network.RData")##network_mat
network_hetero_mat = network_mat[feature_name_list[[1]],feature_name_list[[2]]]
network_gene_mat = network_mat[feature_name_list[[1]],feature_name_list[[1]]]
network_met_mat = network_mat[feature_name_list[[2]],feature_name_list[[2]]]

###################################################
###    get covering pairs on pathway commons ######
###################################################
#load univariate divergence
load("/Users/karenxu/Documents/MultiModalProject/MultiModalPC/data/div_gep_uni.RData")
load("/Users/karenxu/Documents/MultiModalProject/MultiModalPC/data/div_met_uni.RData")

#get g-m pairs
pc_gm_pairs_1step<-matrix2pairs(network_hetero_mat)
temp<-network_gene_mat%*%network_hetero_mat
temp<-ifelse(temp>=1,1,0)
pc_gm_pairs_2step<-matrix2pairs(temp)
temp<-network_gene_mat%*%temp
temp<-ifelse(temp>=1,1,0)
pc_gm_pairs_3step<-matrix2pairs(temp)
pc_gm_pairs<-rbind(pc_gm_pairs_1step,pc_gm_pairs_2step,pc_gm_pairs_3step)%>%unique()
dim(pc_gm_pairs)#1880969 pairs

#filter pairs using chisquare
p<-c()
for (i in 1:nrow(pc_gm_pairs)){
  setTxtProgressBar(pb,i/nrow(pc_gm_pairs))
  gep=abs(div_gep$Mat.div[pc_gm_pairs$Var1[i],])
  met=abs(div_met$Mat.div[pc_gm_pairs$Var2[i],])
  tab=table(gep,met)
  if(nrow(tab)==2 & ncol(tab)==2){
    p<-c(p,chisq.test(tab)$p.value)
  }else{p<-c(p,NULL)}
}
pc_gm_pairs_selected<-pc_gm_pairs[p<=0.05,]#32971 pairs
#filter pairs with rare divergence
pc_gm_pairs_selected<-pc_gm_pairs_selected%>%
  left_join(div_gep$features.div,by=c("Var1"="feature"))%>%
  dplyr::rename(prob.div.Var1=prob.div)%>%
  left_join(div_met$features.div,by=c("Var2"="feature"))%>%
  dplyr::rename(prob.div.Var2=prob.div)
prob<-c()
for (i in 1:nrow(pc_gm_pairs_selected)){
  if(i%%5000==0){print(i)}
  gep=abs(div_gep$Mat.div[pc_gm_pairs_selected$Var1[i],])
  met=abs(div_met$Mat.div[pc_gm_pairs_selected$Var2[i],])
  tab=table(gep,met)
  if(nrow(tab)==2 & ncol(tab)==2){
    prob<-c(prob,tab[2,2]/94)
  }else{prob<-c(prob,0)}
}
pc_gm_pairs_selected$prob.div.pair=prob

pc_selected_pair<-pc_gm_pairs_selected%>%filter(prob.div.pair>=0.02)
pc_selected_source<-pc_gm_pairs_selected%>%filter(prob.div.Var1>=0.02)%>%dplyr::select(Var1)%>%unique()
pc_selected_target<-pc_gm_pairs_selected%>%filter(prob.div.Var2>=0.02)%>%dplyr::select(Var2)%>%unique()
#save(pc_selected_pair,file = "pc_selected_pair.RData")#3679 #---------all final pairs after filtering
#save(pc_selected_source,file = "pc_selected_source.RData")#12245 #---------all final pairs after filtering
#save(pc_selected_target,file = "pc_selected_target.RData")#71 #---------all final pairs after filtering

###################################################
###    optimize coverings on pathway commons ######
###################################################
#construct network from 3679 pairs
Pair_name=pc_selected_pair%>%unite(col="Pair",Var1,Var2)%>%dplyr::select(Pair)
Pair_name=Pair_name$Pair
#build div matrix as pair*sample
pc_selected_pair_mat<-matrix(0,ncol = 94,nrow = nrow(pc_selected_pair),
                             dimnames = list(Pair_name,colnames(div_gep$Mat.div)))
for (i in 1:nrow(pc_selected_pair_mat)){
  if(i%%500==0){print(i)}
  div_pair_<-abs(div_gep$Mat.div[pc_selected_pair$Var1[i],])+abs(div_met$Mat.div[pc_selected_pair$Var2[i],])
  div_pair_<-ifelse(div_pair_==2,1,0)
  pc_selected_pair_mat[i,] <-div_pair_
}
set.seed(1)
#optimization step---source-target pairs  
covers_pair<-Cover_Pair(pc_selected_pair_mat)#cover pair as list object
length(covers_pair$sol)#20 pairs
cover_features<-unique(unlist(str_split(covers_pair$sol,"_")));cover_features #vector
sum(cover_features%in%rownames(normal_met))#12 metabolites
sum(cover_features%in%rownames(normal_gep))#20 genes

#covering network with source and target
CoveringNetwork_sparse<-network_mat[cover_features,cover_features]
#net_cover_pairs_df<-ggnetwork(net_cover_pairs)
#net_cover_pairs_df<-net_cover_pairs_df%>%
#  mutate(feature_type=ifelse(vertex.names%in%rownames(normal_met),"Metabolites","Gene"))

#fill up intermediate genes and metabolites
cover_features<-str_split(covers_pair$sol,"_")
graph=graph_from_adjacency_matrix(network_mat)
cover_features_extended=unique(names(unlist(lapply(cover_features, 
                                                   function(x){shortest_paths(graph,from = x[1],
                                                                              to=x[2])$vpath[[1]]}))))

CoveringNetwork<-network_mat[cover_features_extended,cover_features_extended]
save(CoveringNetwork,file = "/Users/karenxu/Documents/MultiModalProject/MultiModalPC/data/CoveringNetwork.RData")







