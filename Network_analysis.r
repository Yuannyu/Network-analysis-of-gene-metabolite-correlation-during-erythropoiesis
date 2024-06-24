#################Integrate RNA-Seq and metabolism to do network analysis##############################
setwd('/endosome/archive/CRI/Xu_lab/Data/RNA-seq/Erythroid_Metabolome')
intensity=read.table('Metabolic_Gene_List/BM/New_Merge_CRI70_106_BM_reformat.txt',header=T,sep="\t",row.names=1)
Replicate=c(1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,5,5,5,5,5,5,6,6,6,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8)
stage_num=length(unique(Replicate))
median_intensity=matrix(NA,dim(intensity)[1],stage_num,dimnames=list(rownames(intensity),c('BFUe','CFUe','S1','S2','S3','S4','S5','S6')))
for (i in 1:stage_num){
	median_intensity[,i]=apply(intensity[,Replicate==i],1,median)
}

setwd('/endosome/archive/CRI/Xu_lab/Data/RNA-seq/Erythroid_Metabolome/HTSeq_output')
norm_readcount=read.table('BM_norm_readcount_all_stage_final.txt',header=T,sep="\t",row.names=1)
Replicate=c(1,1,1,2,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7,7,7,7,8,8,8,8)  

stage_num=length(unique(Replicate))
median_norm_readcount=matrix(NA,dim(norm_readcount)[1],stage_num,dimnames=list(rownames(norm_readcount),c('BFUe','CFUe','S1','S2','S3','S4','S5','S6')))
for (i in 1:stage_num){
	median_norm_readcount[,i]=apply(norm_readcount[,Replicate==i],1,median)
}

#exclude BFUe
median_intensity=median_intensity[,2:dim(median_intensity)[2]]
median_norm_readcount=median_norm_readcount[,2:dim(median_norm_readcount)[2]]

RNA_metab_cor_p=matrix(NA,dim(median_intensity)[1],dim(median_norm_readcount)[1],dimnames=list(rownames(median_intensity),rownames(median_norm_readcount)))
RNA_metab_cor_r=matrix(NA,dim(median_intensity)[1],dim(median_norm_readcount)[1],dimnames=list(rownames(median_intensity),rownames(median_norm_readcount)))
for (i in 1:dim(median_intensity)[1]){
	for (j in 1:dim(median_norm_readcount)[1]){
		temp=cor.test(median_intensity[i,],median_norm_readcount[j,],method='spearman',exact=FALSE)
		RNA_metab_cor_p[i,j]=temp$p.value
		RNA_metab_cor_r[i,j]=temp$estimate
	}
}
RNA_metab_cor_p[is.na(RNA_metab_cor_p)]=1

p_cutoff=0.05
for_cyto=NULL
for (i in 1:dim(RNA_metab_cor_p)[1]){
	meta=rownames(RNA_metab_cor_p)[i]
	sig_idx=RNA_metab_cor_p[i,]<p_cutoff
	cor_gene=colnames(RNA_metab_cor_p)[sig_idx]
	for_cyto_temp=cbind(rep(meta,length(cor_gene)),cor_gene,RNA_metab_cor_p[i,sig_idx],RNA_metab_cor_r[i,sig_idx])
	for_cyto=rbind(for_cyto,for_cyto_temp)
}
colnames(for_cyto)=c('Metabolites','metab_genes','P_value','R')
setwd('/endosome/archive/CRI/Xu_lab/Data/RNA-seq/Erythroid_Metabolome')
write.table(for_cyto,'Integrate_RNA_metabolism/Metabolites_cor_metab_genes_exclude_BFUe.txt',row.names=F,col.names=T,sep="\t",quote=F)

setwd('/endosome/archive/CRI/Xu_lab/Data/RNA-seq/Erythroid_Metabolome')
metab=read.table('Integrate_RNA_metabolism/CPDB_pathways_metabolites.tab',header=T,sep="\t",comment.char='')
metab_path=NULL
for (i in 1:dim(metab)[1]){
	metab_id=sapply(strsplit(unlist(strsplit(as.vector(metab[i,3]),',')),':'),function(x) x[2])
	metab_path=rbind(metab_path,cbind(metab_id,as.vector(metab[i,1])))
}	
metab_path=unique(metab_path)
metab_id=read.table('Metabolic_Gene_List/metabolites_name.txt',header=T,sep="\t")[,c(1,3)]
metab_path=merge(metab_path,metab_id,by.x='metab_id',by.y='KEGG_ID')

gene=read.table('Integrate_RNA_metabolism/CPDB_pathways_genes.tab',header=T,sep="\t",comment.char='')
gene_path=NULL
for (i in 1:dim(gene)[1]){
	gene_id=unlist(strsplit(as.vector(gene[i,4]),','))
	gene_path=rbind(gene_path,cbind(gene_id,as.vector(gene[i,1])))
}	
gene_path=unique(gene_path); colnames(gene_path)=c('gene_id','path_id')
#write.table(gene_path,'CPDB_gene2path.txt',row.names=F,col.names=T,sep="\t",quote=F)

metab_path_gene=merge(gene_path,metab_path,by.x='path_id',by.y='V2')
human2mouse=read.table('/project/CRI/Xu_lab/shared/YuannyuZhang/data/homologene/homologene_human2mouse_clean.txt',header=F,sep="\t")[,c(2,4)]
colnames(human2mouse)=c('human_geneid','mouse_geneid')
metab_path_gene=merge(metab_path_gene,human2mouse,by.x='gene_id',by.y='human_geneid')
#write.table(metab_path_gene,'CPDB_metab_path_gene.txt',row.names=F,col.names=T,sep="\t",quote=F)

metab_cor_gene=read.table('Integrate_RNA_metabolism/Metabolites_cor_metab_genes_exclude_BFUe.txt',header=T,sep="\t")
metab_attri=table(metab_cor_gene[,2]); metab_attri=metab_attri[metab_attri>0]
gene_attri=table(metab_cor_gene[,1]); gene_attri=gene_attri[gene_attri>0]
node_attri=rbind(cbind(metab_attri,1),cbind(gene_attri,2))
sorted_node_attri=node_attri[order(node_attri[,1],decreasing=T),]
head(sorted_node_attri)

metab_cor_gene_conf=unique(merge(metab_cor_gene,metab_path_gene,by.x=c('Metabolites','metab_genes'),by.y=c('Compound.Name','mouse_geneid')))
main_path=read.table('Integrate_RNA_metabolism/CPDB_human_KEGG_main_pathways.txt',header=T,sep="\t")
metab_cor_gene_conf2=unique(merge(metab_cor_gene_conf,main_path,by.x='path_id',by.y='pathway'))
write.table(metab_cor_gene_conf2,'Integrate_RNA_metabolism/Metabolites_cor_metab_genes_inKEGG_main_path_exclude_BFUe.txt',row.names=F,col.names=T,sep="\t",quote=F)

metab_attri=table(metab_cor_gene_conf[,2]); metab_attri=metab_attri[metab_attri>0]
gene_attri=table(metab_cor_gene_conf[,1]); gene_attri=gene_attri[gene_attri>0]
node_attri=rbind(cbind(metab_attri,1),cbind(gene_attri,2))
sorted_node_attri=node_attri[order(node_attri[,1],decreasing=T),]
head(sorted_node_attri)
#write.table(node_attri,'Integrate_RNA_metabolism/Metabolites_cor_metab_genes_inKEGG_attri.txt',row.names=T,col.names=T,sep="\t",quote=F)

uni_metab_cor_gene_conf2=unique(metab_cor_gene_conf2[,c(2,3)])
metab_attri=table(uni_metab_cor_gene_conf2[,1]); metab_attri=metab_attri[metab_attri>0]
gene_attri=table(uni_metab_cor_gene_conf2[,2]); gene_attri=gene_attri[gene_attri>0]
node_attri=rbind(cbind(metab_attri,1),cbind(gene_attri,2))
sorted_node_attri=node_attri[order(node_attri[,1],decreasing=T),]
head(sorted_node_attri)
#write.table(node_attri,'Integrate_RNA_metabolism/Metabolites_cor_metab_genes_inKEGG_main_path_attri.txt',row.names=T,col.names=T,sep="\t",quote=F)


metab_cor_gene_conf=read.table('Integrate_RNA_metabolism/Metabolites_cor_metab_genes_inKEGG_main_path_exclude_BFUe.txt',header=T,sep="\t")
metab_3cluster=read.table('Metabolic_Gene_List/BM/maSigPro/BM_maSigPro_sig_gene_3cluster.txt',header=T,sep="\t")
gene_3cluster=read.table('maSigPro/BM_maSigPro_sig_gene_3cluster_rmRn45s.txt',header=T,sep="\t")

metab_cor_gene_conf3=merge(metab_cor_gene_conf,gene_3cluster,by.x='metab_genes',by.y=0)
uni_metab_cor_gene_conf3=unique(metab_cor_gene_conf3[,c(1,3)])
metab_attri=table(uni_metab_cor_gene_conf3[,2]); metab_attri=metab_attri[metab_attri>0]
gene_attri=table(uni_metab_cor_gene_conf3[,1]); gene_attri=gene_attri[gene_attri>0]
node_attri=rbind(cbind(metab_attri,1),cbind(gene_attri,2))
sorted_node_attri=node_attri[order(node_attri[,1],decreasing=T),]
head(sorted_node_attri)
#write.table(node_attri,'Integrate_RNA_metabolism/Metabolites_cor_metab_genes_inKEGG_gene_3cluster_attri.txt',row.names=T,col.names=T,sep="\t",quote=F)

metab_cor_gene_conf4=merge(metab_cor_gene_conf3,metab_3cluster,by.x='Metabolites',by.y=0)
uni_metab_cor_gene_conf4=unique(metab_cor_gene_conf4[,c(2,1)])
metab_attri=table(uni_metab_cor_gene_conf4[,2]); metab_attri=metab_attri[metab_attri>0]
gene_attri=table(uni_metab_cor_gene_conf4[,1]); gene_attri=gene_attri[gene_attri>0]
node_attri=rbind(cbind(metab_attri,1),cbind(gene_attri,2))
sorted_node_attri=node_attri[order(node_attri[,1],decreasing=T),]
head(sorted_node_attri)

#colnames(node_attri)=c('Node','diffusion_output_rank','label')
#write.table(node_attri,'Integrate_RNA_metabolism/Metabolites_cor_metab_genes_inKEGG_metab_olp_3cluster_attri.txt',row.names=T,col.names=T,sep="\t",quote=F)
write.table(metab_cor_gene_conf4,'Integrate_RNA_metabolism/Metabolites_cor_metab_genes_inKEGG_main_path_olp_3cluster_exclude_BFUe.txt',row.names=F,col.names=T,sep="\t",quote=F)

#metab_attri=table(metab_cor_gene_conf4[,2]); metab_attri=metab_attri[metab_attri>0]
#gene_attri=table(metab_cor_gene_conf4[,1]); gene_attri=gene_attri[gene_attri>0]
#node_attri=rbind(cbind(metab_attri,1),cbind(gene_attri,2))
#write.table(node_attri,'Integrate_RNA_metabolism/Metabolites_cor_metab_genes_inKEGG_main_path_3cluster_attri.txt',row.names=T,col.names=T,sep="\t",quote=F)

setwd('/endosome/archive/CRI/Xu_lab/Data/RNA-seq/Erythroid_Metabolome')
network=read.table('Integrate_RNA_metabolism/Metabolites_cor_metab_genes_inKEGG_main_path_olp_3cluster_exclude_BFUe.txt',header=T,sep="\t")
metab_degree=table(as.vector(network[,1]))
metab_norm_degree=sort(round((metab_degree/length(unique(network[,2])))*100,1),decreasing=T)
gene_degree=table(as.vector(network[,2]))
gene_norm_degree=sort(round((gene_degree/length(unique(network[,1])))*100,1),decreasing=T)
metab_norm_degree_over10=metab_norm_degree[metab_norm_degree>10]
gene_norm_degree_over10=gene_norm_degree[gene_norm_degree>10]
pdf("Integrate_RNA_metabolism/norm_degree_exclude_BFUe.pdf")
plot(as.vector(metab_norm_degree),col='grey',cex=log2(as.vector(metab_norm_degree)+0.1)*0.5,pch=16,xlab='metabolites or genes',ylab="Normalized Connectivity Score")
points(as.vector(gene_norm_degree),col='grey',cex=log2(as.vector(gene_norm_degree)+0.1)*0.5,pch=16)
points(as.vector(gene_norm_degree_over10),col='red',cex=log2(as.vector(gene_norm_degree)+0.1)*0.5,pch=16)
points(as.vector(metab_norm_degree_over10),col='blue',cex=log2(as.vector(metab_norm_degree)+0.1)*0.5,pch=16)
dev.off()

pdf("Integrate_RNA_metabolism/norm_degree_exclude_BFUe_with_label.pdf")
plot(as.vector(metab_norm_degree),col='grey',cex=log2(as.vector(metab_norm_degree)+0.1)*0.5,pch=16,xlab='metabolites or genes',ylab="Normalized Connectivity Score")
points(as.vector(gene_norm_degree),col='grey',cex=log2(as.vector(gene_norm_degree)+0.1)*0.5,pch=16)
points(as.vector(gene_norm_degree_over10),col='red',cex=log2(as.vector(gene_norm_degree)+0.1)*0.5,pch=16)
points(as.vector(metab_norm_degree_over10),col='blue',cex=log2(as.vector(metab_norm_degree)+0.1)*0.5,pch=16)
text(1:length(metab_norm_degree_over10),as.vector(metab_norm_degree_over10),names(metab_norm_degree_over10),cex=0.5)
text(1:length(gene_norm_degree_over10),as.vector(gene_norm_degree_over10),names(gene_norm_degree_over10),cex=0.5)
dev.off()
norm_degree=rbind(cbind(metab_attri,1),cbind(gene_attri,2))
#################Integrate RNA-Seq and metabolism##############################