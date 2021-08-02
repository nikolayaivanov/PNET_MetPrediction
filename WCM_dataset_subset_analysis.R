#
library(DESeq2)
library(tximeta)
library(pheatmap)
library(RColorBrewer)

source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')
load('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/wcm_metadata.rda')
load('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/validated_genes_DEbyMetStatus_concordantLFC_only.rda')

# import data
se = tximeta(coldata=wcm_metadata, type = "salmon")
# found matching transcriptome:
# [ Ensembl - Homo sapiens - release 97 ]

# summarize transcript-level quantifications to gene-level
gse = summarizeToGene(se)

# make DESeqDataSet object
dds = DESeqDataSet(gse, design = ~ Distant_Mets)

# make transformed count data, using variance stabilizing transformation (VST)
vsd = vst(dds, blind=FALSE)
vst_counts = as.matrix(assay(vsd))

all(colnames(vst_counts)==wcm_metadata$names) #TRUE

# Patients from whom samples CA20 and CA35 were derived had localized disease at the time of sequencing, but later developed distant mets; add that info to the metadata
wcm_metadata$Recurrence=NA
wcm_metadata$Recurrence[which(wcm_metadata$Distant_Mets=='Y')]=NA
wcm_metadata$Recurrence[which(wcm_metadata$names=='CA20')]='Yes'
wcm_metadata$Recurrence[which(wcm_metadata$names=='CA35')]='Yes'

library(AnnotationHub)
hub = AnnotationHub()
dm = query(hub, c("EnsDb", "sapiens", "97"))
edb = dm[["AH73881"]]
genes=as.data.frame(genes(edb))

mm_validated_genes_w_concordantLFCs=match(validated_genes_concordantLFC_only$Ensembl, rownames(vst_counts))

if(length(which(is.na(mm_validated_genes_w_concordantLFCs))) != 0) { 
	validated_genes_concordantLFC_only = validated_genes_concordantLFC_only[-which(is.na(mm_validated_genes_w_concordantLFCs)),]
	mm_validated_genes_w_concordantLFCs=match(validated_genes_concordantLFC_only$Ensembl, rownames(vst_counts))
 }

which(duplicated(validated_genes_concordantLFC_only$Ensembl))# none
which(duplicated(rownames(vst_counts)))# none

mm=match(rownames(vst_counts),genes$gene_id)
which(duplicated(rownames(vst_counts)))# none
which(duplicated(genes$gene_id))# none
which(is.na(mm)) # none

rownames(vst_counts) = genes$gene_name[mm]

all(colnames(vst_counts)==wcm_metadata$names) #TRUE

annotation_col=data.frame(Distant.Mets=wcm_metadata$Distant_Mets, Loc.to.Met=wcm_metadata$Recurrence)
rownames(annotation_col)=as.vector(wcm_metadata$names)

ann_colors = list(
	Loc.to.Met = c(Yes='blue'),
	Distant.Mets = c(Y='green', N='purple')
	)

pdf('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_pdfs/WCM_dataset_subset_analysis_heatmap.pdf')

pheatmap(vst_counts[mm_validated_genes_w_concordantLFCs,], color=colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(3000), main='Subset analysis; WCM dataset', 
	cluster_rows=FALSE, cluster_cols=TRUE, annotation_col=annotation_col, scale='row', fontsize_col=10, annotation_colors = ann_colors)

aa=vst_counts[mm_validated_genes_w_concordantLFCs,]

if(length(which(rowVars(aa) < 1e-20) != 0)){ 
	aa=aa[-which(rowVars(aa)< 1e-20),] 
}

pheatmap(aa, color=colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(3000), main='Subset analysis; WCM dataset', 
		cluster_rows=TRUE, cluster_cols=TRUE, annotation_col=annotation_col, scale='row', fontsize_col=10, show_rownames=TRUE, annotation_colors = ann_colors)

# pheatmap(vst_counts, color=colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(3000), main='Subset analysis; WCM dataset', 
# 	cluster_rows=FALSE, cluster_cols=TRUE, annotation_col=annotation_col, scale='row', fontsize_col=10, show_rownames=FALSE, annotation_colors = ann_colors)

dev.off()

################ PCA

# pcaData <- plotPCA(vsd, intgroup = c( "Distant_Mets"), returnData = TRUE)
# percentVar <- round(100 * attr(pcaData, "percentVar"))

library(calibrate)

pca = prcomp(t(assays(vsd)[[1]]))
#command that will return % of variance explained by each PC:
pcaVars=getPcaVars(pca)

all(rownames(pca$x) == wcm_metadata$names) #TRUE

PCs=as.matrix(pca$x)

## PCA colored by metastatic status

col=wcm_metadata$Distant_Mets
col=sub('N','purple',col)
col=sub('Y','green',col)

pdf('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_pdfs/wcm_datset_PCA_plots_metStatus.pdf')

par(mar=c(5.1,5.3,4.1,2.1))

for(i in seq(from=1, to=10, by=2)){
    x_lab=paste0('PC',i,': ',signif(pcaVars[i],2),'% variance')
    y_lab=paste0('PC',i+1,': ',signif(pcaVars[i+1],2),'% variance')
    plot(PCs[,i], PCs[,i+1], xlab=x_lab, ylab=y_lab, pch=21, cex=1.2, col='black', bg=col, cex.lab=2, cex.axis=2, xlim=c(min(PCs[,i])-10,max(PCs[,i])+10))
    textxy(PCs[,i], PCs[,i+1],wcm_metadata$names,cex =.7, offset = .7)
    legend("topleft", legend=c('Distant mets','No distant mets'), pch=15, col=c("green", "purple"),  cex=1, pt.cex=1)
}

dev.off()






################################### DRAFTS



pdf('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_pdfs/WCM_dataset_subset_analysis_heatmap_EXPERIMENTAL.pdf')

#
pheatmap(vst_counts[mm_validated_genes_w_concordantLFCs,], color=colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(3000), main='Subset analysis; WCM dataset', 
	cluster_rows=FALSE, cluster_cols=TRUE, annotation_col=annotation_col, scale='row', fontsize_col=10)

	pheatmap(vst_counts, color=colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(3000), main='Subset analysis; WCM dataset', 
	cluster_rows=FALSE, cluster_cols=TRUE, annotation_col=annotation_col, scale='row', fontsize_col=10, show_rownames=FALSE)

#
aa=vst_counts[mm_validated_genes_w_concordantLFCs,]

if(length(which(rowVars(aa)==0) != 0)){ 
	rownames(aa)[which(rowVars(aa)==0)]
	aa=aa[-which(rowVars(aa)==0),] 
}

pheatmap(aa, color=colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(3000), main='Subset analysis; WCM dataset', 
		cluster_rows=TRUE, cluster_cols=TRUE, annotation_col=annotation_col, scale='row', fontsize_col=10, show_rownames=FALSE)

#
if ( length(which(rowVars(vst_counts)==0)) != 0 ){
	pheatmap(vst_counts[-which(rowVars(vst_counts)==0),], color=colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(3000), main='Subset analysis; WCM dataset', 
	cluster_rows=TRUE, cluster_cols=TRUE, annotation_col=annotation_col, scale='row', fontsize_col=10, show_rownames=FALSE)
} else{
	pheatmap(vst_counts, color=colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(3000), main='Subset analysis; WCM dataset', 
	cluster_rows=TRUE, cluster_cols=TRUE, annotation_col=annotation_col, scale='row', fontsize_col=10, show_rownames=FALSE)
}


dev.off()





