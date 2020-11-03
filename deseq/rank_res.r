setwd("~/Documents/infection_netowrk_meta_analysis/")
library(Biobase)

lfc_mat=read.table("./deseq/res/meta_log2FC_mat_2019.txt",header = T)
padj_mat=read.table("./deseq/res/meta_padj_mat_2019.txt",header = T)

lfc_mat[is.na(lfc_mat)]=0
padj_mat[is.na(padj_mat)]=1


keep=which((rowMin(as.matrix(padj_mat))<0.001)&(rowMax(abs(as.matrix(lfc_mat)))>1)) 


lfc_mat=lfc_mat[keep,]
padj_mat=padj_mat[keep,]

rnk_mat=matrix(NA, ncol = dim(lfc_mat)[2], nrow = dim(lfc_mat)[1])
rownames(rnk_mat)=rownames(lfc_mat)
colnames(rnk_mat)=colnames(lfc_mat)
for (i in seq(ncol(lfc_mat))){
  rnk_mat[,i]=as.numeric(rank(-lfc_mat[,i]*-log(padj_mat[,i]),ties.method = "first"))
}


select=which(rowVars(rnk_mat)>summary(rowVars(rnk_mat))[4])

rnk_mat_flt=rnk_mat[select,]

colnames(rnk_mat_flt)=c("DENV1","DENV2","HSV1","EBOLA","HRV","HIV","RSV","STAPH","IAV","HCMV")
for (i in seq(ncol(rnk_mat_flt))){
  rnk_mat_flt[,i]=as.numeric(rank(rnk_mat_flt[,i],ties.method = "first"))
}

write.table(rnk_mat_flt,"./deseq/res/rank_matrix_padj_cutoff_0001_rowVar.txt",row.names = T, col.names = T, quote = F, sep = "\t")
