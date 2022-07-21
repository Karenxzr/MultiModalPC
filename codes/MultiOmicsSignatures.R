library(survival)
library(stringr)
library(plyr)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(survminer)


source("/Users/karenxu/Documents/MultiModalProject/MultiModalPC/codes/utilities.R")
load("/Users/karenxu/Documents/MultiModalProject/MultiModalPC/data/feature_name_list.RData")
load("/Users/karenxu/Documents/MultiModalProject/MultiModalPC/data/tumor_gm.RData")
load("/Users/karenxu/Documents/MultiModalProject/MultiModalPC/data/CoveringNetwork.RData")
load("/Users/karenxu/Documents/MultiModalProject/MultiModalPC/data/survival_df.RData")



#reorder cover features
cover_features = rownames(CoveringNetwork)
cover_features_g=cover_features[cover_features%in%feature_name_list[[1]]]
cover_features_m=cover_features[cover_features%in%feature_name_list[[2]]]
cover_features<-c(cover_features_g,cover_features_m)
CoveringNetwork<-CoveringNetwork[cover_features,cover_features]###pc

#identify signatures from gm pairs by diffusion
pairs_gm=matrix2pairs(CoveringNetwork[cover_features_g,cover_features_m])
label_mat=matrix(nrow = ncol(CoveringNetwork))
for (i in 1:nrow(pairs_gm)){
  emptymat=matrix(0,nrow = ncol(CoveringNetwork),ncol = 1)
  rownames(emptymat)=c(cover_features_g,cover_features_m)
  emptymat[unlist(pairs_gm[i,]),1]=1
  label_mat=cbind(label_mat,emptymat)
}
label_mat=label_mat[,-1]
diffuse_label_df<-diffus_label(network = CoveringNetwork,
                               label_mat = label_mat,iter=100,alpha = 0.75)##diffusion

library(ggdendro)
dend=hclust(dist(t(diffuse_label_df)),method = "ward.D")
ggdendrogram(dend)+
  theme(axis.text.x = element_text(size=8))+
  geom_hline(yintercept = 0.45,color='red',linetype='dashed')+
  geom_hline(yintercept = 0.65,color='red',linetype='dashed')+
  geom_hline(yintercept = 0.85,color='red',linetype='dashed')

###get candidate signatures
colnames(diffuse_label_df)=paste0('edge',1:ncol(diffuse_label_df))
candidate_signature=unique(do.call(c,list(signaturefromcluster(edge_cluster0,label_mat),
                                          signaturefromcluster(edge_cluster1,label_mat),
                                          signaturefromcluster(edge_cluster2,label_mat))))

##########grouping
library(survival)
library(survminer)
#---------ward.D
p_sig=c()
for (i in 1:length(candidate_signature)){
  sig_div=groupbysignature(tumor,candidate_signature[[i]],'ward.D',return_rownames = T)
  sig_div$uuid=str_remove(rownames(sig_div),'-T')
  sig_div=sig_div%>%left_join(survival_df)%>%filter(!is.na(Year))
  fit0=survfit(Surv(Year,delta)~group,data = sig_div)
  p_sig=c(p_sig,surv_pvalue(fit0)$pval)
}
candidate_signature1 = candidate_signature[which(p_sig<0.05)]
#---------complete
p_sig=c()
for (i in 1:length(candidate_signature)){
  sig_div=groupbysignature(tumor,candidate_signature[[i]],'complete',return_rownames = T)
  sig_div$uuid=str_remove(rownames(sig_div),'-T')
  sig_div=sig_div%>%left_join(survival_df)%>%filter(!is.na(Year))
  fit0=survfit(Surv(Year,delta)~group,data = sig_div)
  p_sig=c(p_sig,surv_pvalue(fit0)$pval)
}
candidate_signature2 = candidate_signature[which(p_sig<0.05)]
candidate_signature = unique(c(candidate_signature1,candidate_signature2))

##############permutation
#"EGLN3"  "trans.4.hydroxyproline" "succinate"  power: 95.6%
#"IL6"       "SLC22A2"   "histamine" 94.5%
#"EGLN3""GNB1""trans.4.hydroxyproline""succinate""histamine"94.5%   ---unbalance      
#"GNB1"      "succinate" "histamine" 93.3%---unbalance

i=1;print(candidate_signature[[i]])
exp=t(tumor[candidate_signature[[i]],])
rownames(exp)=str_remove(rownames(exp),'-T')
p_sig = c()
#permutation
#signature1
for (k in 1:1000){
  set.seed(k)
  exp_perm = apply(exp, 2, sample,nrow(exp));rownames(exp_perm)=rownames(exp)
  group = groupbysignature(t(exp_perm),candidate_signature[[i]],method="ward.D")
  group = group%>%left_join(survival_df)
  fit0=survfit(Surv(Year,delta)~group,data = group)
  p_sig=c(p_sig,surv_pvalue(fit0)$pval)
}
p_sig=ifelse(is.na(p_sig),1,p_sig)
1-sum(p_sig<=0.05)/1000

#signature2
i=4;print(candidate_signature[[i]])
exp=t(tumor[candidate_signature[[i]],])
rownames(exp)=str_remove(rownames(exp),'-T')
p_sig = c()
#permutation
for (k in 1:1000){
  set.seed(k)
  exp_perm = apply(exp, 2, sample,nrow(exp));rownames(exp_perm)=rownames(exp)
  group = groupbysignature(t(exp_perm),candidate_signature[[i]],method="complete")
  group = group%>%left_join(survival_df)
  fit0=survfit(Surv(Year,delta)~group,data = group)
  p_sig=c(p_sig,surv_pvalue(fit0)$pval)
}
p_sig=ifelse(is.na(p_sig),1,p_sig)
1-sum(p_sig<=0.05)/1000


###final signatures
survival_df_signature=cbind(groupbysignature(tumor,candidate_signature[[1]],'ward.D',return_rownames = T),
              groupbysignature(tumor,candidate_signature[[4]],'complete',return_rownames = F))
colnames(survival_df_signature)[2:3]=paste0("signature",1:2)
survival_df_signature$uuid=str_remove(rownames(survival_df_signature),'-T')
survival_df_signature=survival_df_signature%>%left_join(survival_df)%>%mutate(Year = months_since_rp/12)
#signature ensemble
survival_df_signature$signature_ensemble = ifelse(survival_df_signature$signature1=='Group 2'&
                                               survival_df_signature$signature2=='Group 2','Group 2','Group 1')
save(survival_df_signature,file='/Users/karenxu/Documents/MultiModalProject/MultiModalPC/data/survival_df_signature.RData')

