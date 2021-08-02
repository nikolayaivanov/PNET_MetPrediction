#

######### discovery dataset

# print TPM matrix
load('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/discovery_dataset_DESeq2_DEbyMetStatus_BUNDLE.rda') # se, gse, tpm, counts, ddsClean, rr, vsd

# add gene symbols
library(AnnotationHub)
hub = AnnotationHub()
dm = query(hub, c("EnsDb", "sapiens", "97"))
edb = dm[["AH73881"]]
genes=as.data.frame(genes(edb))

mm=match(rownames(tpm),genes$gene_id)
length(which(is.na(mm))) # 0
gene.names=as.vector(genes$gene_name[mm])
new_rownames=paste0(gene.names,' (',rownames(tpm),')')
rownames(tpm)=new_rownames

write.csv(tpm, file='/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_data_tables/discovery_dataset_tpm_values.csv')

# AGPAT2 and PLIN5
index=c(which(rr$gene=='AGPAT2'),which(rr$gene=='PLIN5'))
rr[index,]
#                  baseMean log2FoldChange     lfcSE      stat    pvalue
#                 <numeric>      <numeric> <numeric> <numeric> <numeric>
# ENSG00000169692  4191.310      -0.229166  0.372996 -0.614393  0.538956
# ENSG00000214456   132.763      -0.215759  0.583010 -0.370078  0.711325
#                      padj         chr         Ensembl        gene
#                 <numeric> <character>     <character> <character>
# ENSG00000169692  0.999998           9 ENSG00000169692      AGPAT2
# ENSG00000214456  0.999998          19 ENSG00000214456       PLIN5

######### validation dataset
rm(list=ls())

# print TPM matrix
load('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/validation_dataset_DESeq2_DEbyMetStatus_BUNDLE.rda') # se, gse, tpm, counts, ddsClean, rr, vsd

# add gene symbols
library(AnnotationHub)
hub = AnnotationHub()
dm = query(hub, c("EnsDb", "sapiens", "97"))
edb = dm[["AH73881"]]
genes=as.data.frame(genes(edb))

mm=match(rownames(tpm),genes$gene_id)
length(which(is.na(mm))) # 0
gene.names=as.vector(genes$gene_name[mm])
new_rownames=paste0(gene.names,' (',rownames(tpm),')')
rownames(tpm)=new_rownames

write.csv(tpm, file='/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_data_tables/validation_dataset_tpm_values.csv')

# AGPAT2 and PLIN5
index=c(which(rr$gene=='AGPAT2'),which(rr$gene=='PLIN5'))
rr[index,]
#                  baseMean log2FoldChange     lfcSE      stat    pvalue
#                 <numeric>      <numeric> <numeric> <numeric> <numeric>
# ENSG00000169692  1519.860       0.667119  0.345916  1.928558 0.0537858
# ENSG00000214456    57.248       0.145939  0.563546  0.258965 0.7956621
#                      padj         chr         Ensembl        gene
#                 <numeric> <character>     <character> <character>
# ENSG00000169692  0.323931           9 ENSG00000169692      AGPAT2
# ENSG00000214456  0.999621          19 ENSG00000214456       PLIN5


# NAI








