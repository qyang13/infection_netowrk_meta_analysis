setwd("~/Documents/infection_netowrk_meta_analysis/deseq/")

library(DESeq2)

# Function definition for deseq2
runDeseq2=function(data, conds){
  # Create pesudo replicates if no replicates present
  if (ncol(data)==2){
    temp=data
    temp[,2]=round(data[,1]*runif(length(data[,1]), min=0.5, max=1.5))
    temp[,3]=data[,2]
    temp[,4]=round(data[,2]*runif(length(data[,1]), min=0.5, max=1.5))
    data=temp
    names=c('m1',"m2","m3","m4")
    colnames(data)=names
    tempC=conds
    tempC[2,]=tempC[1,]
    tempC=rbind(tempC,conds[2,])
    tempC=rbind(tempC,conds[2,])
    rownames(tempC)=names
    conds=tempC
    conds[,1]=rownames(conds)
  }
  
  cond=as.data.frame(conds)
  colnames(cond)=c('srr','condition')
  cond=cond[which(cond$srr%in%colnames(data)),]
  cond=cond[order(match(cond$srr,colnames(data))),]
  rownames(cond)=cond[,1]

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = data,
                                        colData = cond,
                                        design = ~ condition)
  
  dds$condition <- factor(dds$condition, levels = c("mock","infected"))
  dds <- DESeq2::estimateSizeFactors(dds)
  dds = DESeq2::estimateDispersions(dds)
  DESeq2::plotDispEsts(dds)
  dds = DESeq2::nbinomWaldTest (dds)
  return(data.frame(results(dds)))
}

# Meta-data
meta=read.csv("meta_data.csv", header = T)

# Read the featurecount results
cts=read.table("read_counts_unstranded.rc",header = T)
colnames(cts)=gsub("X.scratch.Users.qiya9811.meta_analysis.bam.sorted.","",colnames(cts))
colnames(cts)=gsub(".sorted.bam","",colnames(cts))
for (a in seq(2,5)) {
  cts[,a]=gsub(";.*","",cts[,a])
}
dat=(cts)
rownames(dat)=dat$Geneid
dat=dat[,-1:-6]
keep <- rowSums(dat) >= 100
dat <- dat[keep,]
dat=dat[,-(which(colSums(dat)==0))]

# Initialize the lfc and padj matrix
lfc_mat=data.frame(matrix(nrow=nrow(dat), ncol=0))
padj_mat=data.frame(matrix(nrow=nrow(dat), ncol=0))
rownames(lfc_mat)=rownames(dat)
rownames(padj_mat)=rownames(dat)
processed_SRP=''
# Iteratively analyze all datasets, populate the matrices
for (srp in unique(meta$SRP.Index)) {
  processed_SRP=c(processed_SRP, srp)
  # srp="SRP152882"
  srpIdxs=which(meta$SRP.Index%in%srp)
  
  for (cellType in unique(meta$Cell.Type[srpIdxs])) {
    cellIdxs=intersect(which(meta$Cell.Type%in%cellType),srpIdxs)
    
    for (virus in unique(meta$Pathogen[cellIdxs])) {
      print(paste("Processing:",srp,"Cell Type:",cellType,"Virus:",virus))
      
      virusIdx=intersect(which(meta$Pathogen%in%virus),cellIdxs)
      srrs=meta$SRR.Accession.Number[virusIdx]
      conds=cbind(as.character(srrs),as.character(meta$Mock.infected[virusIdx]))
      
      res=runDeseq2(dat[,which(colnames(dat)%in%as.character(srrs))],conds)
      lfc_mat=cbind(lfc_mat,res$log2FoldChange)
      padj_mat=cbind(padj_mat,res$padj)
      colnames(lfc_mat)[ncol(lfc_mat)]=paste(virus,cellType,sep="_")
      colnames(padj_mat)[ncol(padj_mat)]=paste(virus,cellType,sep="_")
    }
  }
}

# Save the matrices
write.table(lfc_mat, "meta_log2FC_mat_v2.txt",sep="\t",row.names = T, col.names = T, quote = F)
write.table(padj_mat, "meta_padj_mat_v2.txt",sep="\t",row.names = T, col.names = T, quote = F)
