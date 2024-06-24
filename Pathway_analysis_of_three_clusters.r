###########do pathway enrichment using CPDB pathway for all genes in three clusters.
library('ggplot2')
library("stringr")
setwd('/endosome/archive/CRI/Xu_lab/Data/RNA-seq/Erythroid_Metabolome')
anno=unique(read.table('/endosome/archive/CRI/Xu_lab/Data/RNA-seq/Erythroid_Metabolome/CPDB_gene2path.txt',header=T,sep="\t"))
## manually change heme biosythesis to Heme biosythesis
main_path=read.table('Integrate_RNA_metabolism/CPDB_human_KEGG_main_pathways.txt',header=T,sep="\t")
mouse_symbol=read.table('/project/CRI/Xu_lab/shared/YuannyuZhang/data/homologene/homologene_human2mouse_clean.txt',header=F,sep="\t")[,c(2,4)]
anno=merge(anno,mouse_symbol,by.x='gene_id',by.y='V2')
uni_path=as.vector(unique(anno[,2]))
group=c('All_gene')
# use CPDB database to do pathway enrichment.
## do path enrichment using CPDB database.
norm_readcount=read.table('HTSeq_output/BM_norm_readcount_all_stage_final.txt',header=T,sep="\t",row.names=1)
norm_readcount=norm_readcount[apply(norm_readcount,1,max)>=10,]
sig_gene=read.table('maSigPro/BM_maSigPro_sig_gene_3cluster_rmRn45s.txt',header=T,sep="\t")
sig_gene[sig_gene[,1]==1,]=0; sig_gene[sig_gene[,1]==2,]=1; sig_gene[sig_gene[,1]==0,]=2  #switch cluster 1 and cluster 2.
all_gene=rownames(norm_readcount)
N=length(all_gene)
str_wrap_factor <- function(x,L) {
  levels(x) <- str_wrap(levels(x),width=L)
  x
}
##keep the order when wrap the string.
for (j in 1:3){
	enrich_path=matrix(NA,length(uni_path),7)
    enrich_path[,1]=uni_path
	k=NULL; n=NULL; M=NULL;
	for (w in 1:dim(enrich_path)[1]){
		gene_in_path=intersect(unique(as.vector(anno[anno[,2]==enrich_path[w,1],3])),all_gene)
		M[w]=length(gene_in_path)
		gene_in_cluster=rownames(sig_gene)[sig_gene[,1]==j] 
		n[w]=length(gene_in_cluster)
		k[w]=length(intersect(gene_in_path,gene_in_cluster))
	}			
	p=1-phyper(k-1,M,N-M,n)
	padj=p.adjust(p,method='BH')
	ratio=round(k/M,2)
	ER=round(k*N/(M*n),2)
	enrich_path[,2:7]=cbind(N,M,n,ER,padj,ratio)
	enrich_path=enrich_path[order(padj),]
	colnames(enrich_path)=c('Path_class','#Gene_NormRead_over10','#Gene_in_Path','#Gene_in_Cluster','Enrichment_Ratio','adjusted_P_value','Percentage')
	enrich_path=as.data.frame(enrich_path)
	enrich_path[,5]=as.numeric(as.vector(enrich_path[,5])); enrich_path[,6]=as.numeric(as.vector(enrich_path[,6])); enrich_path[,7]=as.numeric(as.vector(enrich_path[,7])); enrich_path[,3]=as.numeric(as.vector(enrich_path[,3]))
	#enrich_path[,1]=sapply(strsplit(as.vector(enrich_path[,1]),'-'),function(x) x[1])
	#enrich_path=enrich_path[enrich_path[,3]<1000,]
	#enrich_path$Path_class <- factor(enrich_path$Path_class,levels=enrich_path$Path_class)  #lock in factor level order
	enrich_path[enrich_path[,6]==0,6]=10^-14
	#enrich_path[,6]=-log10(enrich_path[,6])
	
	uni_padj=unique(enrich_path[,6])
	enrich_path_ordered=NULL
	for (n in 1:length(uni_padj)){
		enrich_path_temp=enrich_path[enrich_path[,6]==uni_padj[n],]
		enrich_path_temp_ordered=enrich_path_temp[order(enrich_path_temp[,5]),]
		enrich_path_ordered=rbind(enrich_path_ordered,enrich_path_temp_ordered)
	}
	## if p value is the same. Rank by ER score.
	
	##ranke by p value. p value as size. ER as dot color.
	for_fig=enrich_path_ordered[1:20,]; for_fig[,7]=-log10(for_fig[,6]); for_fig=for_fig[order(for_fig[,7]),]
	if (j==1){for_fig[for_fig[,5]>2,5]=2}
	if (j==2){for_fig[for_fig[,5]>3,5]=3}
	if (j==3){for_fig[for_fig[,5]>6,5]=6}
	colnames(for_fig)=c('Path_class','#Gene_NormRead_over10','#Gene_in_Path','#Gene_in_Cluster','Enrichment_Ratio','adjusted_P_value','-log10(adjusted_P_value)')
	for_fig$Path_class <- factor(for_fig$Path_class,levels=for_fig$Path_class)
	
	pdf(paste0('maSigPro/BM_',group[1],'_maSigPro_group',j,'_enrich_CPDB_path_top20_rank_by_p_value.pdf'))
	print(
		ggplot(for_fig,aes(x =-log10(adjusted_P_value) , y = str_wrap_factor(Path_class,60), colour = Enrichment_Ratio , size = -log10(adjusted_P_value) )) + 
		geom_point() +
		scale_colour_gradient(low = "mistyrose", high = "red") +
		theme_bw() +
		###changing xticks
		if(j == 1) {
			scale_x_continuous(breaks = c(4, 8, 12, 16), lim = c(2, 16))
		} else if(j == 2) {
			scale_x_continuous(breaks = c(1, 3, 5, 7), lim = c(0.5, 7))
		} else if(j == 3) {
			scale_x_continuous(breaks = c(11, 13, 15), lim = c(11, 15))
		}
	)
	dev.off()
	
	#write.table(enrich_path,paste0('maSigPro/BM_',group[1],'_maSigPro_group',j,'_enrich_CPDB_path.txt'),row.names=F,col.names=T,sep="\t",quote=F)
	
	main_enrich_path=merge(main_path,enrich_path,by.x='pathway',by.y='Path_class')
	#write.table(main_enrich_path,paste0('maSigPro/BM_',group[1],'_maSigPro_group',j,'_enrich_human_KEGG_path.txt'),row.names=F,col.names=T,sep="\t",quote=F)
}


## do pathway enrichment for metabolic genes in three clusters.
library('ggplot2')
library("stringr")
setwd('/endosome/archive/CRI/Xu_lab/Data/RNA-seq/Erythroid_Metabolome')
anno=unique(read.table('/endosome/archive/CRI/Xu_lab/Data/RNA-seq/Erythroid_Metabolome/CPDB_gene2path.txt',header=T,sep="\t"))
main_path=read.table('Integrate_RNA_metabolism/CPDB_human_KEGG_main_pathways.txt',header=T,sep="\t")
mouse_symbol=read.table('/project/CRI/Xu_lab/shared/YuannyuZhang/data/homologene/homologene_human2mouse_clean.txt',header=F,sep="\t")[,c(2,4)]
anno=merge(anno,mouse_symbol,by.x='gene_id',by.y='V2')
uni_path=as.vector(unique(anno[,2]))
group=c('metabolic_gene')
# use CPDB database to do pathway enrichment.
## do path enrichment using CPDB database.
norm_readcount=read.table('HTSeq_output/BM_norm_readcount_all_stage_final.txt',header=T,sep="\t",row.names=1)
norm_readcount=norm_readcount[apply(norm_readcount,1,max)>=10,]
metabolic=read.table('Metabolic_Gene_List/Sabatini_Metabolic_Gene_List_Human.txt',header=T,sep="\t")
human2mouse=read.table('/project/CRI/Xu_lab/shared/YuannyuZhang/data/homologene/homologene_human2mouse_clean.txt',header=F,sep="\t")
metabolic_id=unique(as.vector(human2mouse[human2mouse[,1]%in%metabolic[,2],4]))
norm_readcount=norm_readcount[rownames(norm_readcount)%in%metabolic_id,]
N=length(metabolic_id)
sig_gene=read.table('maSigPro/BM_maSigPro_sig_gene_3cluster_rmRn45s.txt',header=T,sep="\t")
sig_gene[sig_gene[,1]==1,]=0; sig_gene[sig_gene[,1]==2,]=1; sig_gene[sig_gene[,1]==0,]=2  #switch cluster 1 and cluster 2.
temp=unique(merge(sig_gene,metabolic_id,by.x=0,by.y=1))
rownames(temp)=temp[,1]; sig_gene=temp[-1]
all_gene=rownames(norm_readcount)
str_wrap_factor <- function(x,L) {
  levels(x) <- str_wrap(levels(x),width=L)
  x
}
##keep the order when wrap the string.
for (j in 1:3){
	enrich_path=matrix(NA,length(uni_path),7)
    enrich_path[,1]=uni_path
	k=NULL; n=NULL; M=NULL;
	for (w in 1:dim(enrich_path)[1]){
		gene_in_path=intersect(unique(as.vector(anno[anno[,2]==enrich_path[w,1],3])),all_gene)
		M[w]=length(gene_in_path)
		gene_in_cluster=rownames(sig_gene)[sig_gene[,1]==j] 
		n[w]=length(gene_in_cluster)
		k[w]=length(intersect(gene_in_path,gene_in_cluster))
	}			
	p=1-phyper(k-1,M,N-M,n)
	padj=p.adjust(p,method='BH')
	ratio=round(k/M,2)
	ER=round(k*N/(M*n),2)
	enrich_path[,2:7]=cbind(N,M,n,ER,padj,ratio)
	enrich_path=enrich_path[order(padj),]
	enrich_path=enrich_path[enrich_path[,3]>2,]
	colnames(enrich_path)=c('Path_class','#Gene_NormRead_over10','#Gene_in_Path','#Gene_in_Cluster','Enrichment_Ratio','adjusted_P_value','Percentage')
	enrich_path=as.data.frame(enrich_path)
	enrich_path[,5]=as.numeric(as.vector(enrich_path[,5])); enrich_path[,6]=as.numeric(as.vector(enrich_path[,6])); enrich_path[,7]=as.numeric(as.vector(enrich_path[,7])); enrich_path[,3]=as.numeric(as.vector(enrich_path[,3]))
	#enrich_path[,1]=sapply(strsplit(as.vector(enrich_path[,1]),'-'),function(x) x[1])
	enrich_path=enrich_path[enrich_path[,3]<500,]
	#enrich_path$Path_class <- factor(enrich_path$Path_class,levels=enrich_path$Path_class)  #lock in factor level order
	enrich_path[enrich_path[,6]==0,6]=10^-14
	#enrich_path[,6]=-log10(enrich_path[,6])
	if (j==3){
		enrich_path_top=enrich_path[enrich_path[,6]<1,]
		enrich_path_bottom=enrich_path[enrich_path[,6]==1,]
		enrich_path_bottom_ranked=enrich_path_bottom[order(enrich_path_bottom[,5],decreasing=T),]
		enrich_path=rbind(enrich_path_top,enrich_path_bottom_ranked)
	}
	
	uni_padj=unique(enrich_path[,6])
	enrich_path_ordered=NULL
	for (n in 1:length(uni_padj)){
		enrich_path_temp=enrich_path[enrich_path[,6]==uni_padj[n],]
		enrich_path_temp_ordered=enrich_path_temp[order(enrich_path_temp[,5],decreasing=T),]
		enrich_path_ordered=rbind(enrich_path_ordered,enrich_path_temp_ordered)
	}
	## if p value is the same. Rank by ER score.
	
    ##ranke by p value. p value as size. ER as dot color.
	for_fig=enrich_path_ordered[20:1,]; for_fig[,7]=-log10(for_fig[,6]); 
	if (j==1){for_fig[for_fig[,5]>3,5]=3}
	if (j==2){for_fig[for_fig[,5]>8,5]=8}
	if (j==3){for_fig[for_fig[,5]>20,5]=20}
	colnames(for_fig)=c('Path_class','#Gene_NormRead_over10','#Gene_in_Path','#Gene_in_Cluster','Enrichment_Ratio','adjusted_P_value','-log10(adjusted_P_value)')
	for_fig$Path_class <- factor(for_fig$Path_class,levels=for_fig$Path_class)
	pdf(paste0('maSigPro/BM_',group[1],'_maSigPro_group',j,'_enrich_CPDB_path_top20_rank_by_p_value.pdf'))
	print(
		ggplot(for_fig,aes(x =-log10(adjusted_P_value) , y = str_wrap_factor(Path_class,60), colour = Enrichment_Ratio , size = -log10(adjusted_P_value) )) + 
		geom_point() +
		scale_colour_gradient(low = "mistyrose", high = "red") +
		theme_bw() +
		###changing xticks
		if(j == 1) {
			scale_x_continuous(breaks = c(1,3,5), lim = c(0, 7))
		} else if(j == 2) {
			scale_x_continuous(breaks = c(0,1,2), lim = c(0, 3))
		} else if(j == 3) {
			scale_x_continuous(breaks = c(0, 0.5, 1), lim = c(0, 1))
		}
	)
	dev.off()	
	
	for_fig=enrich_path_ordered[20:1,]; for_fig[,7]=-log10(for_fig[,6]); 
	if (j==1){for_fig[for_fig[,5]>3,5]=3; for_fig=for_fig[-18,]; for_fig=for_fig[15:19,]; color_low="mistyrose"}
	if (j==2){for_fig[for_fig[,5]>8,5]=8; for_fig=for_fig[-17,]; for_fig=for_fig[15:19,]; color_low="mistyrose"}
	if (j==3){for_fig[for_fig[,5]>20,5]=20;for_fig=for_fig[-16,]; for_fig=for_fig[15:19,]; color_low="red"}
	colnames(for_fig)=c('Path_class','#Gene_NormRead_over10','#Gene_in_Path','#Gene_in_Cluster','Enrichment_Ratio','adjusted_P_value','-log10(adjusted_P_value)')
	for_fig$Path_class <- factor(for_fig$Path_class,levels=for_fig$Path_class)
	pdf(paste0('maSigPro/BM_',group[1],'_maSigPro_group',j,'_enrich_CPDB_path_top5_rank_by_p_value.pdf'),width=6, height=1.5,useDingbats=FALSE)
	print(
		ggplot(for_fig,aes(x =-log10(adjusted_P_value) , y = str_wrap_factor(Path_class,60), colour = Enrichment_Ratio , size = -log10(adjusted_P_value) )) + 
		geom_point() +
		scale_colour_gradient(low = color_low, high = "red") +
		theme_bw() +
		###changing xticks
		if(j == 1) {
			scale_x_continuous(breaks = c(1,3,5), lim = c(0, 7))
		} else if(j == 2) {
			scale_x_continuous(breaks = c(0,1,2), lim = c(0, 3))
		} else if(j == 3) {
			scale_x_continuous(breaks = c(0, 0.5, 1), lim = c(0, 1))
		}
	)
	dev.off()
	
	pdf(paste0('maSigPro/BM_',group[1],'_maSigPro_group',j,'_enrich_CPDB_path_top5_rank_by_p_value_color_bar.pdf'),width=6, height=3,useDingbats=FALSE)
	print(
		ggplot(for_fig,aes(x =-log10(adjusted_P_value) , y = str_wrap_factor(Path_class,60), colour = Enrichment_Ratio , size = -log10(adjusted_P_value) )) + 
		geom_point() +
		scale_colour_gradient(low = color_low, high = "red") +
		theme_bw() +
		###changing xticks
		if(j == 1) {
			scale_x_continuous(breaks = c(1,3,5), lim = c(0, 7))
		} else if(j == 2) {
			scale_x_continuous(breaks = c(0,1,2), lim = c(0, 3))
		} else if(j == 3) {
			scale_x_continuous(breaks = c(0, 0.5, 1), lim = c(0, 1))
		}
	)
	dev.off()
	
	write.table(enrich_path,paste0('maSigPro/BM_',group[1],'_maSigPro_group',j,'_enrich_CPDB_path.txt'),row.names=F,col.names=T,sep="\t",quote=F)
	
	main_enrich_path=merge(main_path,enrich_path,by.x='pathway',by.y='Path_class')
	write.table(main_enrich_path,paste0('maSigPro/BM_',group[1],'_maSigPro_group',j,'_enrich_human_KEGG_path.txt'),row.names=F,col.names=T,sep="\t",quote=F)
}

###########do pathway enrichment