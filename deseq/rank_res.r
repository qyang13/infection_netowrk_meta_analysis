setwd("./")

lfc_mat=read.table("./deseq/res/meta_log2FC_mat_2019.txt",header = T)
padj_mat=read.table("./deseq/res/meta_padj_mat_2019.txt",header = T)

lfc_mat[is.na(lfc_mat)]=0
padj_mat[is.na(padj_mat)]=1

rnk_mat=matrix(NA, ncol = dim(lfc_mat)[2], nrow = dim(lfc_mat)[1])
rownames(rnk_mat)=rownames(lfc_mat)
colnames(rnk_mat)=colnames(lfc_mat)
for (i in seq(ncol(lfc_mat))){
  rnk_mat[,i]=as.numeric(rank(-lfc_mat[,i]*-log(padj_mat[,i]),ties.method = "first"))
}

colnames(rnk_mat)=c("DENV1","DENV2","HSV1","EBOLA","HRV","HIV","RSV","STAPH","IAV","HCMV")

write.table(rnk_mat,"./deseq/res/rank_matrix.txt",row.names = T, col.names = T, quote = F, sep = "\t")
