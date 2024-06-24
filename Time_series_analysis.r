####do time series analysis for all genes.
#module load R/3.6.1-gccmkl
library(maSigPro)
setwd('/endosome/archive/CRI/Xu_lab/Data/RNA-seq/Erythroid_Metabolome')
norm_readcount=read.table('HTSeq_output/BM_norm_readcount_all_stage_final.txt',header=T,sep="\t",row.names=1)
sample_info=matrix(NA,dim(norm_readcount)[2],3,dimnames=list(colnames(norm_readcount),c('Time','Replicate','Group')))

sample_info[,1]=c(0,0,0,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,6,6,6,7,7,7,7)
sample_info[,2]=c(1,1,1,2,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7,7,7,7,8,8,8,8)
sample_info[,3]=1

design <- make.design.matrix(sample_info, degree = 4)
fit <- p.vector(norm_readcount, design, Q = 0.01, MT.adjust = "BH", min.obs = 3,counts=TRUE)
tstep <- T.fit(fit, step.method = "backward", alfa = 0.01)
setwd('/endosome/archive/CRI/Xu_lab/Data/RNA-seq/Erythroid_Metabolome/maSigPro')
#save.image(file = "maSigPro_tstep.RData")
# To load the data again
#load("maSigPro_tstep.RData")

sigs <- get.siggenes(tstep, rsq = 0.6, vars = "groups")
setwd('/endosome/archive/CRI/Xu_lab/Data/RNA-seq/Erythroid_Metabolome/maSigPro')
pdf('BM_maSigPro_output_rmRn45s.pdf')
sig_gene=see.genes(sigs$sig.genes$Group, show.fit = T, dis =design$dis,cluster.method="hclust" ,cluster.data = 1, k = 3,newX11=FALSE)
dev.off()
unsure=colnames(tstep$influ.info)
gene_3cluster=sig_gene$cut
output=gene_3cluster[!(names(gene_3cluster)%in%unsure)]
write.table(output,'BM_maSigPro_sig_gene_3cluster_rmRn45s.txt',row.names=T,col.names=T,sep="\t",quote=F)
####do time series analysis for all genes and identify three clusters of genes.

###########draw heatmap using three clusters of genes.
setwd('/endosome/archive/CRI/Xu_lab/Data/RNA-seq/Erythroid_Metabolome/maSigPro')
sig_gene=read.table('BM_maSigPro_sig_gene_3cluster_rmRn45s.txt',header=T,sep="\t",row.names=1)
sig_gene[sig_gene[,1]==1,]=0; sig_gene[sig_gene[,1]==2,]=1; sig_gene[sig_gene[,1]==0,]=2  #switch cluster 1 and cluster 2.
sig_gene_readcount=merge(norm_readcount,sig_gene,by.x=0,by.y=0)
sig_gene_readcount_mtx=as.matrix(sig_gene_readcount[,2:(dim(sig_gene_readcount)[2]-1)])
rownames(sig_gene_readcount_mtx)=sig_gene_readcount[,1]
#sig_gene_readcount_mtx=log2(sig_gene_readcount_mtx+1)
Replicate=c(1,1,1,2,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7,7,7,7,8,8,8,8)
stage_num=length(unique(Replicate))
median_sig_gene_readcount_mtx=matrix(NA,dim(sig_gene_readcount_mtx)[1],stage_num,dimnames=list(rownames(sig_gene_readcount_mtx),c('BFU-E','CFU-E','ProE','BasoE','PolyE','OrthoE','Retic','RBC')))
for (i in 1:stage_num){
	median_sig_gene_readcount_mtx[,i]=apply(sig_gene_readcount_mtx[,Replicate==i],1,median)
}

max_exp_per=(median_sig_gene_readcount_mtx)*100./apply(median_sig_gene_readcount_mtx,1,max)
max_exp_per=as.matrix(max_exp_per)

library("ComplexHeatmap")
library("gplots")
library("RColorBrewer")
library('pROC')
library('colorspace')
# https://jokergoo.github.io/ComplexHeatmap-reference/book/a-list-of-heatmaps.html  Heatmap instruction
source('/project/CRI/Xu_lab/Data/RNA-seq/Erythroid_Metabolome/colorRamp2.r')
col_ori = colorRamp2(c(0, max(as.matrix(log2(median_sig_gene_readcount_mtx+1)))), c("white", "red"))
col_per = colorRamp2(c(0, 100), c("white", "red"))
cluster_flag=as.vector(sig_gene_readcount[,dim(sig_gene_readcount)[2]])
log_split <- paste0("Cluster\n",cluster_flag)
log_hmap <- Heatmap(as.matrix(log2(median_sig_gene_readcount_mtx+1)), col=col_ori,cluster_row_slices=FALSE,split=log_split,row_names_gp = gpar(fontsize = 0),column_names_gp = gpar(fontsize = 10),cluster_columns = FALSE)
max_split <- log_split
percent_hmap <- Heatmap(max_exp_per, col=col_per, cluster_row_slices=FALSE,split=max_split,row_names_gp = gpar(fontsize = 0),column_names_gp = gpar(fontsize = 10),row_dend_reorder = TRUE,cluster_columns = FALSE)
pdf('BM_maSigPro_3groups_merge_replicates_white_red.pdf')
pushViewport(viewport(layout=grid.layout(nr=1, nc=2)))
  pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
	draw(log_hmap, newpage=FALSE)
  upViewport()
  pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
	draw(percent_hmap, newpage=FALSE)
  upViewport()
upViewport()
dev.off()
###########draw heatmap using three clusters of genes.

###########draw heatmap using three clusters of metabolic genes.
metab_gene=read.table('/endosome/archive/CRI/Xu_lab/Data/RNA-seq/Erythroid_Metabolome/FPKM_TPM_Table/metabolic_gene.txt',header=T,sep="\t")[,1]
median_sig_gene_readcount_mtx=median_sig_gene_readcount_mtx[rownames(median_sig_gene_readcount_mtx)%in%metab_gene,]
sig_gene_readcount=sig_gene_readcount[sig_gene_readcount[,1]%in%metab_gene,]

max_exp_per=(median_sig_gene_readcount_mtx)*100./apply(median_sig_gene_readcount_mtx,1,max)
max_exp_per=as.matrix(max_exp_per)

library("ComplexHeatmap")
library("gplots")
library("RColorBrewer")
# https://jokergoo.github.io/ComplexHeatmap-reference/book/a-list-of-heatmaps.html  Heatmap instruction
col_ori = colorRamp2(c(0, max(as.matrix(log2(median_sig_gene_readcount_mtx+1)))), c("white", "red"))
col_per = colorRamp2(c(0, 100), c("white", "red"))
cluster_flag=as.vector(sig_gene_readcount[,dim(sig_gene_readcount)[2]])
log_split <- paste0("Cluster\n",cluster_flag)
log_hmap <- Heatmap(as.matrix(log2(median_sig_gene_readcount_mtx+1)),col=col_ori, cluster_row_slices=FALSE,split=log_split,row_names_gp = gpar(fontsize = 0),column_names_gp = gpar(fontsize = 10),cluster_columns = FALSE)
max_split <- log_split
percent_hmap <- Heatmap(max_exp_per, col=col_per, cluster_row_slices=FALSE,split=max_split,row_names_gp = gpar(fontsize = 0),column_names_gp = gpar(fontsize = 10),row_dend_reorder = TRUE,cluster_columns = FALSE)
pdf('BM_maSigPro_3groups_merge_replicates_metabolic_gene_white_red.pdf')
pushViewport(viewport(layout=grid.layout(nr=1, nc=2)))
  pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
	draw(log_hmap, newpage=FALSE)
  upViewport()
  pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
	draw(percent_hmap, newpage=FALSE)
  upViewport()
upViewport()
dev.off()
###########draw heatmap using three clusters of metabolic genes.


############### do time series analysis for metabolites in BM using maSigPro and draw heatmap.
setwd('/endosome/archive/CRI/Xu_lab/Data/RNA-seq/Erythroid_Metabolome/Metabolic_Gene_List/BM')
intensity_raw=t(read.table('New_Merge_CRI70_106_BM.txt',header=F,sep="\t"))
intensity=intensity_raw[3:dim(intensity_raw)[1],2:dim(intensity_raw)[2]]
rownames(intensity)=intensity_raw[3:dim(intensity_raw)[1],1]
mode(intensity)='numeric'
col_temp=intensity_raw[2,2:dim(intensity_raw)[2]]
col_name=NULL
uni_col_name=unique(col_temp)
Time=NULL; Replicate=NULL
for (i in 1:length(uni_col_name)){
	rep_num=sum(col_temp==uni_col_name)
	col_name=c(col_name,paste0(uni_col_name[i],'_rep',1:rep_num))
	Time=c(Time,rep(i-1,rep_num))
	Replicate=c(Replicate,rep(i,rep_num))
}
colnames(intensity)=col_name
norm_readcount=intensity

library(maSigPro)
sample_info=matrix(NA,dim(norm_readcount)[2],3,dimnames=list(colnames(norm_readcount),c('Time','Replicate','Group')))
sample_info[,1]=Time
sample_info[,2]=Replicate
sample_info[,3]=1

design <- make.design.matrix(sample_info, degree = 4)
fit <- p.vector(norm_readcount, design, Q = 0.01, MT.adjust = "BH", min.obs = 3)
tstep <- T.fit(fit, step.method = "backward", alfa = 0.01)
sigs <- get.siggenes(tstep, rsq = 0, vars = "groups")
setwd('/endosome/archive/CRI/Xu_lab/Data/RNA-seq/Erythroid_Metabolome/Metabolic_Gene_List/BM/maSigPro')
pdf('BM_maSigPro_output.pdf')
sig_gene=see.genes(sigs$sig.genes$Group, show.fit = T, dis =design$dis,cluster.method="hclust" ,cluster.data = 1, k = 3,newX11=FALSE)
dev.off()
#write.table(sig_gene$cut,'BM_maSigPro_sig_gene_3cluster.txt',row.names=T,col.names=T,sep="\t",quote=F)

sig_gene=read.table('maSigPro/BM_maSigPro_sig_gene_3cluster.txt',header=T,sep="\t",row.names=1)
sig_gene[sig_gene[,1]==1,]=0; sig_gene[sig_gene[,1]==3,]=1; sig_gene[sig_gene[,1]==0,]=3  #switch cluster 1 and cluster 3.
sig_gene_readcount=merge(norm_readcount,sig_gene,by.x=0,by.y=0)
sig_gene_readcount_mtx=as.matrix(sig_gene_readcount[,2:(dim(sig_gene_readcount)[2]-1)]); 
rownames(sig_gene_readcount_mtx)=sig_gene_readcount[,1]
sig_gene_readcount_mtx[sig_gene_readcount_mtx<2^(-5)]=2^-5
sig_gene_readcount_mtx[sig_gene_readcount_mtx>2^(5)]=2^5

max_exp_per=(sig_gene_readcount_mtx)*100./apply(sig_gene_readcount_mtx,1,max)
max_exp_per=as.matrix(max_exp_per)

## draw the heatmap to show normalized readcount in merged replicate.
stage_num=length(unique(Replicate))
median_sig_gene_readcount_mtx=matrix(NA,dim(sig_gene_readcount_mtx)[1],stage_num,dimnames=list(rownames(sig_gene_readcount_mtx),c('BFU-E','CFU-E','ProE','BasoE','PolyE','OrthoE','Retic','RBC')))
for (i in 1:stage_num){
	median_sig_gene_readcount_mtx[,i]=apply(sig_gene_readcount_mtx[,Replicate==i],1,median)
}
setwd('/endosome/archive/CRI/Xu_lab/Data/RNA-seq/Erythroid_Metabolome/Metabolic_Gene_List/BM/maSigPro')
#pdf('BM_sig_metabolites_intensity_boxplot.pdf')
#boxplot(log2(median_sig_gene_readcount_mtx), main='Sig metabolites',las=2)
#dev.off()

median_sig_gene_readcount_mtx[median_sig_gene_readcount_mtx<2^(-5)]=2^-5
median_sig_gene_readcount_mtx[median_sig_gene_readcount_mtx>2^(5)]=2^5
# From the boxplot,I can tell most value is between [2^-5,2^5]. So I set the value below or over the interval to the edge value.

max_exp_per=(median_sig_gene_readcount_mtx)*100./apply(median_sig_gene_readcount_mtx,1,max)
max_exp_per=as.matrix(max_exp_per)

library("ComplexHeatmap")
library("gplots")
library("RColorBrewer")
library('pROC')
library('colorspace')
# https://jokergoo.github.io/ComplexHeatmap-reference/book/a-list-of-heatmaps.html  Heatmap instruction
source('/project/CRI/Xu_lab/Data/RNA-seq/Erythroid_Metabolome/colorRamp2.r')
col_ori = colorRamp2(c(0, 1), c("white", "blue"))
col_per = colorRamp2(c(0, 100), c("white", "blue"))
cluster_flag=as.vector(sig_gene_readcount[,dim(sig_gene_readcount)[2]])
log_split <- paste0("Cluster\n",cluster_flag)
log_hmap <- Heatmap(as.matrix(log2(median_sig_gene_readcount_mtx)), col=col_ori,cluster_row_slices=FALSE,split=log_split,row_names_gp = gpar(fontsize = 2),column_names_gp = gpar(fontsize = 10),cluster_columns = FALSE)
max_split <- log_split
percent_hmap <- Heatmap(max_exp_per, col=col_per,cluster_row_slices=FALSE,split=max_split,row_names_gp = gpar(fontsize = 2),column_names_gp = gpar(fontsize = 10),row_dend_reorder = TRUE,cluster_columns = FALSE)
pdf('BM_maSigPro_3groups_merge_replicates_blue_white.pdf')
pushViewport(viewport(layout=grid.layout(nr=1, nc=2)))
  pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
	draw(log_hmap, newpage=FALSE)
  upViewport()
  pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
	draw(percent_hmap, newpage=FALSE)
  upViewport()
upViewport()
dev.off()
############### generate a heatmap of metabolites in BM using maSigPro