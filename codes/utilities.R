###functions
#1. from adjacency matrix to pairs
matrix2pairs<-function(mat){
  pairs=which(mat>0)
  r_num<-pairs%/%nrow(mat)
  c_num<-pairs%%nrow(mat)
  ind=which(c_num==0)
  c_num[ind]=nrow(mat)
  r_num[ind]=r_num[ind]-1
  pairs<-data.frame(Var1=rownames(mat)[c_num],
                    Var2=colnames(mat)[r_num+1])
  return(pairs)
}



#2.find covering pairs from matrix, mat: pair*sample
Cover_Pair<-function(mat,J=1){
  alpha = 1 - sum(colSums(mat) >= J)/ncol(mat)
  R = cover_lpSolve(mat, alpha=alpha, maxsol=10, J=J)
  mat.sol = R$sol
  union = unique(as.vector(mat.sol))
  pvals = rowMeans(mat[union, ])
  sel.sol = which.max(apply(mat.sol, 2, function(x) sum(pvals[x]) ))[1]
  sol = mat.sol[, sel.sol]
  core = union[which(rowMeans(apply(mat.sol, 2, function(x) union %in% x )) == 1)]
  alpha = sum(colSums(mat[sol, ]) >= J)/ncol(mat)
  list(mat.sol=mat.sol,
       union=union,
       core=core,#the group has max number components from union
       sol=sol,#the solution with maxmized divergence prob
       alpha=alpha)
}


#3. plot patient network
plot_patient<-function(ID,AjacencyMat=net_cover_pairs,
                       BaseNet=graph_cover,
                       BaseGraph=igraph_cover,
                       Pair_Div_Mat=pc_selected_pair_mat,
                       covers_pair=covers_pair){
  #diverged pairs--plot edges
  edge=names(which(Pair_Div_Mat[rownames(Pair_Div_Mat)%in%covers_pair$sol,ID]==1))
  edge = str_split(edge,"_")
  edge = lapply(edge, function(x){names(shortest_paths(BaseGraph,from = x[1],to=x[2])$vpath[[1]])})
  mat_temp<-matrix(0,nrow = nrow(net_cover_pairs),ncol=ncol(net_cover_pairs),
                   dimnames =  dimnames(net_cover_pairs))
  for (i in 1:length(edge)){
    mat_temp[edge[[i]],edge[[i]]]<-net_cover_pairs[edge[[i]],edge[[i]]]
  }
  set.edge.value(BaseNet,'pair_div',mat_temp)
  #diverged nodes--plot nodes
  nodes=names(which(abs(rbind(div_gep$Mat.div,div_met$Mat.div)[unique(unlist(cover_features)),ID])==1))
  BaseNet%v%'node_div'=ifelse(network.vertex.names(BaseNet)%in%nodes,1,0)
  set.seed(1)
  df=ggnetwork(BaseNet)
  p=ggplot(df,aes(x = x, y = y, xend = xend, yend = yend))+theme_blank()+
    geom_edges(color="grey",size=0.3,alpha=0.5)+
    geom_nodes(aes(color=Type),alpha=0.3,size=0.8)+
    theme(legend.position = 'bottom')+
    geom_edges(data=df%>%filter(pair_div==1),color='red',size=1,alpha=0.5)+
    geom_nodelabel(data=df%>%filter(node_div==1),aes(color = Type, label = vertex.names),
                   fontface = "bold", size=2.5,
                   alpha=0.5)
  return(p)
}



#4. label prop
diffus_label=function(network,label_mat,alpha=0.75,iter=10,difference=1e-6){
  require(SMUT)
  #prepocessing
  label_mat=label_mat[match(rownames(network),rownames(label_mat)),]
  dcol=sqrt(colSums(network))
  drow=sqrt(rowSums(network))
  network=t(t(network)/dcol)/drow
  #initialize
  snet_1=label_mat
  snet=snet_1
  #diffusion on adjacency matrix
  for(kk in 1:iter){
    snet_1<-alpha*eigenMapMatMult(network,snet)+(1-alpha)*(label_mat)
    diff=max(abs(snet_1-snet))
    print(c("iteration:",kk,"difference:",diff))
    if(diff<difference){return(snet_1)}
    snet=snet_1
  }
  return(snet_1)
}

#5. generate binary indicator by median
binary_median=function(x){ifelse(x>=median(x),1,0)}


#6. group by signature
groupbysignature<-function(mat,signature,method="ward.D",return_rownames=T){
  dend=hclust(dist(scale(t(mat[signature,]))),method = method)
  dend=cutree(dend, k = 2)
  dend=ifelse(dend==1,'Group 1','Group 2')
  if (return_rownames){df=data.frame(uuid=names(dend),group=dend)}else{
    df=data.frame(group=dend) 
  }
  return(df)
}


#7. generate candidate signature
signaturefromcluster<-function(edge_cluster,label_mat){
  candidate_signature=list()
  count=0
  for (i in 1:max(edge_cluster)){
    edge=as.numeric(str_remove(names(which(edge_cluster==i)),'edge'))
    if(length(edge)>1){candidate_signature[[count+1]]=rownames(label_mat)[which(apply(label_mat[,edge],1,sum)>0)]
    count=count+1}
  }
  return(candidate_signature)
}

#8. generate KTSP scores 
ktsp_score<-function(mat,pairs,group){
  mat_tsp = matrix(NA,ncol = ncol(mat),nrow = nrow(pairs))
  for (i in 1:nrow(pairs)){
    mat_tsp[i,] = mat[pairs[i,1],]>mat[pairs[i,2],]
  }
  colnames(mat_tsp)=colnames(mat)
  rownames(mat_tsp)=paste0(pairs[,1],'_',pairs[,2])
  mat_tsp=mat_tsp+0
  
  mat_tsp_1 = mat_tsp[,which(group==unique(group)[1])]
  mat_tsp_2 = mat_tsp[,which(group==unique(group)[2])]
  scores = abs(apply(mat_tsp_1, 1, sum)/ncol(mat_tsp_1)-apply(mat_tsp_2, 1, sum)/ncol(mat_tsp_2))
  
  return(list(mat_tsp=mat_tsp,scores=scores))
}

