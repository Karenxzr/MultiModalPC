library(doParallel)
library(adabag)
library(survival)
library(survminer)
library(biomaRt)

library(glmnet)
source("/Users/karenxu/Documents/MultiModalProject/MultiModalPC/codes/utilities.R")


#####construct surrogate gene signature
#load data
load("/Users/karenxu/Documents/MultiModalProject/MultiModalPC/data/feature_name_list.RData")
load("/Users/karenxu/Documents/MultiModalProject/MultiModalPC/data/tumor_gm.RData")
load("/Users/karenxu/Documents/MultiModalProject/MultiModalPC/data/survival_df_signature.RData")
signatures = list(signature1 = c("EGLN3","trans.4.hydroxyproline","succinate"),
                  signature2 = c("IL6","SLC22A2","histamine"))
colnames(tumor) = str_remove(colnames(tumor),'-T')
tumor_g = tumor[feature_name_list[[1]],]
tumor_m = tumor[feature_name_list[[2]],]


############################################################
#####step 1: construct dataframe for gene pool##############
############################################################
#signature1
gene_score = list()
for (i in 1:nrow(tumor_g)){
  ttest = t.test(tumor_g[i,survival_df_signature[survival_df_signature$signature1=='Group 2','uuid']],
                 tumor_g[i,survival_df_signature[survival_df_signature$signature1=='Group 1','uuid']])
  m=signatures[[1]][signatures[[1]]%in%feature_name_list[[2]]]
  cor_vec = c()
  for (k in 1:length(m)){
    cortest = cor.test(tumor_g[i,],tumor_m[m[k],])
    cor_vec = c(cor_vec,c(cortest$estimate,cortest$p.value))
  }
  gene_score[[i]]=c(t = ttest$statistic,p= ttest$p.value,cor_vec)
}
gene_score = do.call(rbind,gene_score)
rownames(gene_score)=rownames(tumor_g)
colnames(gene_score) = c('de_t','de_p','cor_trans.4.hydroxyproline',
                         'p_trans.4.hydroxyproline','cor_succinate','p_succinate')
#save(gene_score,file = 'gene_score_sig1.RData')
#signature2
gene_score = list()
for (i in 1:nrow(tumor_g)){
  ttest = t.test(tumor_g[i,survival_df_signature[survival_df_signature$signature2=='Group 2','uuid']],
                 tumor_g[i,survival_df_signature[survival_df_signature$signature2=='Group 1','uuid']])
  m=signatures[[2]][signatures[[2]]%in%feature_name_list[[2]]]
  cor_vec = c()
  for (k in 1:length(m)){
    cortest = cor.test(tumor_g[i,],tumor_m[m[k],])
    cor_vec = c(cor_vec,c(cortest$estimate,cortest$p.value))
  }
  gene_score[[i]]=c(t = ttest$statistic,p= ttest$p.value,cor_vec)
}
gene_score = do.call(rbind,gene_score)
rownames(gene_score)=rownames(tumor_g)
colnames(gene_score) = c('de_t','de_p','cor_histamine','p_histamine')
#save(gene_score,file = 'gene_score_sig2.RData')

#####load from here########
load("/Users/karenxu/Documents/MultiModalProject/MultiModalPC/data/gene_score_sig1.RData")
gene_score_sig1 = gene_score
load("/Users/karenxu/Documents/MultiModalProject/MultiModalPC/data/gene_score_sig2.RData")
gene_score_sig2 = gene_score
rm(gene_score)

##generate pairs-----sig1
gene_score_sig1 = gene_score_sig1%>%as.data.frame()%>%
  filter((de_p<=0.01 & (p_trans.4.hydroxyproline<=0.01|p_succinate<=0.01))|rownames(.)=='EGLN3')%>%
  mutate(gene = rownames(.),regulate = ifelse(de_t>=0,'up','down'))
up = gene_score_sig1%>%filter(regulate=='up')
down = gene_score_sig1%>%filter(regulate=='down')
gene_score_sig1_pairs = as.matrix(expand.grid(up = up$gene,down = down$gene));dim(gene_score_sig1_pairs)#540 pairs
#tsp pairs
ktsp_score_sig1 = ktsp_score(tumor_g,pairs=gene_score_sig1_pairs,group = survival_df_signature$signature1)
ktsp_score_sig1_pairs = ktsp_score_sig1$scores[order(ktsp_score_sig1$scores,decreasing = T)][1:50]
#data for build surrogate classifier
df_data =as.data.frame(t(ktsp_score_sig1$mat_tsp[names(ktsp_score_sig1_pairs),]))
df_data$signature1 = survival_df_signature$risk1
df_data = colwise(as.factor)(df_data)

set.seed(42)

fit1 = cv.glmnet(t(ktsp_score_sig1$mat_tsp[names(ktsp_score_sig1_pairs),]),
                 y=as.numeric(as.factor(survival_df_signature$signature1))-1,
                 family='binomial',type.measure = 'class')

##generate pairs-----sig2
gene_score_sig2 = gene_score_sig2%>%as.data.frame()%>%
  filter((de_p<=0.01 & p_histamine<=0.01)|rownames(.)%in%c('IL6','SLC22A2'))%>%
  mutate(gene = rownames(.),regulate = ifelse(de_t>=0,'up','down'))
up = gene_score_sig2%>%filter(regulate=='up')
down = gene_score_sig2%>%filter(regulate=='down')
gene_score_sig2_pairs = as.matrix(expand.grid(up = up$gene,down = down$gene));dim(gene_score_sig2_pairs)#170 pairs
#tsp pairs
ktsp_score_sig2 = ktsp_score(tumor_g,pairs=gene_score_sig2_pairs,group = sig_div$signature2)
ktsp_score_sig2_pairs = ktsp_score_sig2$scores[order(ktsp_score_sig2$scores,decreasing = T)][1:50]

set.seed(75)
fit2 = cv.glmnet(t(ktsp_score_sig2$mat_tsp[names(ktsp_score_sig2_pairs),]),
                 y=as.numeric(as.factor(survival_df_signature$signature2))-1,
                 family='binomial',type.measure = 'class')

