#
library(DESeq2)
library(tximeta)
library(pheatmap)
library(RColorBrewer)

source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')
load('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/chan_metadata.rda')
load('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/validated_genes_DEbyMetStatus_concordantLFC_only.rda')

# import data
se = tximeta(coldata=chan_metadata, type = "salmon")
# found matching transcriptome:
# [ Ensembl - Homo sapiens - release 97 ]

# summarize transcript-level quantifications to gene-level
gse = summarizeToGene(se)

# make DESeqDataSet object
dds = DESeqDataSet(gse, design = ~ Distant_Mets)

# make transformed count data, using variance stabilizing transformation (VST)
vsd = vst(dds, blind=FALSE)
vst_counts = as.matrix(assay(vsd))

all(colnames(vst_counts)==chan_metadata$names) #TRUE

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

all(colnames(vst_counts)==chan_metadata$names) #TRUE

annotation_col=data.frame(Distant.Mets=chan_metadata$Distant_Mets)
rownames(annotation_col)=as.vector(chan_metadata$names)

pdf('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_pdfs/Chan_dataset_subset_analysis_heatmap.pdf')

ann_colors = list(
	Distant.Mets = c(Y='green', N='purple')
	)

pheatmap(vst_counts[mm_validated_genes_w_concordantLFCs,], color=colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(3000), main='Subset analysis; Chan et al dataset', 
	cluster_rows=FALSE, cluster_cols=TRUE, annotation_col=annotation_col, scale='row', fontsize_col=10, annotation_colors = ann_colors)

pheatmap(vst_counts, color=colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(3000), main='Subset analysis; Chan et al dataset', 
	cluster_rows=FALSE, cluster_cols=TRUE, annotation_col=annotation_col, scale='row', fontsize_col=10, show_rownames=FALSE)

dev.off()



