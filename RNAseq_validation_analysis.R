#
###########################################################################################
##### Validation DE analysis [localized vs metastatic]
###########################################################################################

###### Read in the data and peform differential expression analysis

load('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/validation_dataset_metadata.rda')

source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')
library(DESeq2)
library(tximeta)

# import data
se = tximeta(coldata=validation_dataset_metadata, type = "salmon")
# found matching transcriptome:
# [ Ensembl - Homo sapiens - release 97 ]

# summarize transcript-level quantifications to gene-level
gse = summarizeToGene(se)

# get TPM matrix
tpm = assays(gse)$abundance

#get count matrix
counts=assays(gse)$counts

# make DESeqDataSet object
dds = DESeqDataSet(gse, design = ~ Distant_Mets)

#perform pre-filtering to keep only rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# make transformed count data, using variance stabilizing transformation (VST)
vsd = vst(dds, blind=FALSE)

# run SVA (see chapter 8.1 of https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html)
library(sva)

dds <- estimateSizeFactors(dds) # using 'avgTxLength' from assays(dds), correcting for library size
dat <- counts(dds, normalized=TRUE)
idx = rowMeans(dat) > 1
dat = dat[idx, ]
mod = model.matrix(~ Distant_Mets, colData(dds))
mod0 <- model.matrix(~1, colData(dds))
svaobj = svaseq(dat, mod, mod0)
# Number of significant surrogate variables is:  10

colnames(svaobj$sv)=paste0('SV_',1:ncol(svaobj$sv))

colData(dds) = as(cbind(as.data.frame(colData(gse)),svaobj$sv),'DataFrame')

design(dds) = ~ SV_1 + SV_2 + SV_3 + SV_4 + SV_5 + SV_6 + SV_7 + SV_8 + SV_9 + SV_10 + Distant_Mets

# examine how well the SVA method did at recovering the batch variable (i.e. which dataset the samples originated from)

pdf('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_pdfs/validation_dataset_SVA_plots.pdf')

par(mar=c(5.1,5.3,4.1,2.1))

x=factor(dds$Dataset)

for(i in 1:10){
  boxplot(svaobj$sv[, i] ~ x, xlab='Batch (Dataset)', ylab=paste0("SV", i), main=paste0("Surrogate Variable ", i),
    cex.main=2, cex.lab=2, cex.axis=1.5, outline=FALSE, col='lightgrey', ylim=c( min(as.vector(svaobj$sv[, i])), max(as.vector(svaobj$sv[, i])) ) )
  points(as.vector(svaobj$sv[, i]) ~ jitter(as.numeric(x), amount=0.2), pch =21, col='black', bg='darkgrey', cex=1.4)
}

# for (i in 1:12) { 
#   stripchart(svaobj$sv[, i] ~ dds$Dataset, vertical = TRUE, main = paste0("Surrogate Variable ", i), ylab=paste0("SV", i),xlab='Batch (Dataset)',cex.main=2, cex.lab=2, cex.axis=1.5) 
#   abline(h = 0, lty='dashed') 
# }

dev.off()

# run DE analysis
dds=DESeq(dds)

which(is.na(mcols(dds)$betaConv)) # none

# 11 rows did not converge in beta
# omit rows that did not converge in beta (these are typically genes with very small counts and little power)
# see https://support.bioconductor.org/p/65091/
ddsClean <- dds[which(mcols(dds)$betaConv),]

# extract results
rr=results(ddsClean, alpha=0.1, contrast=c('Distant_Mets','Y','N'))
# contrast = c( the name of a factor in the design formula, name of the numerator level for the fold change, name of the denominator level for the fold change)
summary(rr)
# out of 25837 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 873, 3.4%
# LFC < 0 (down)     : 785, 3%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

# add gene symbols, chr, & Ensembl gene IDs
library(AnnotationHub)
hub = AnnotationHub()
dm = query(hub, c("EnsDb", "sapiens", "97"))
edb = dm[["AH73881"]]
genes=as.data.frame(genes(edb))

mm=match(rownames(rr), genes$gene_id)
length(which(is.na(mm))) # 0
rr$chr=as.vector(genes$seqnames[mm])
rr$Ensembl=as.vector(rownames(rr))
rr$gene=as.vector(genes$gene_name[mm])

# save releveant data
save(se, gse, tpm, counts, ddsClean, rr, vsd, file='/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/validation_dataset_DESeq2_DEbyMetStatus_BUNDLE.rda')
# load('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/validation_dataset_DESeq2_DEbyMetStatus_BUNDLE.rda') # se, gse, tpm, counts, ddsClean, rr, vsd

###########################################################################################
##### Downstream analysis
###########################################################################################

library(DESeq2)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')
load('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/validation_dataset_metadata.rda')
load('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/validation_dataset_DESeq2_DEbyMetStatus_BUNDLE.rda') # se, gse, tpm, counts, ddsClean, rr, vsd

################ PCA

# pcaData <- plotPCA(vsd, intgroup = c( "Distant_Mets"), returnData = TRUE)
# percentVar <- round(100 * attr(pcaData, "percentVar"))

library(calibrate)

pca = prcomp(t(assays(vsd)[[1]]))
#command that will return % of variance explained by each PC:
pcaVars=getPcaVars(pca)

all(rownames(pca$x) == validation_dataset_metadata$names) #TRUE

PCs=as.matrix(pca$x)

## PCA colored by metastatic status

col=validation_dataset_metadata$Distant_Mets
col=sub('N','forestgreen',col)
col=sub('Y','red',col)

pdf('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_pdfs/validation_datset_PCA_plots_metStatus.pdf')

par(mar=c(5.1,5.3,4.1,2.1))

for(i in seq(from=1, to=10, by=2)){
    x_lab=paste0('PC',i,': ',signif(pcaVars[i],2),'% variance')
    y_lab=paste0('PC',i+1,': ',signif(pcaVars[i+1],2),'% variance')
    plot(PCs[,i], PCs[,i+1], xlab=x_lab, ylab=y_lab, pch=21, cex=1.2, col='black', bg=col, cex.lab=2, cex.axis=2, xlim=c(min(PCs[,i])-10,max(PCs[,i])+10))
    textxy(PCs[,i], PCs[,i+1],validation_dataset_metadata$names,cex =.7, offset = .7)
    legend("topright", legend=c('Distant mets','No distant mets'), pch=15, col=c("red", "forestgreen"),  cex=1, pt.cex=1)
}

dev.off()

## PCA colored by study

col=c('blue','darkred')
palette(col)
study=as.factor(validation_dataset_metadata$Dataset)

pdf('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_pdfs/validation_datset_PCA_plots_studyID.pdf')

par(mar=c(5.1,5.3,4.1,2.1))

for(i in seq(from=1, to=10, by=2)){
    x_lab=paste0('PC',i,': ',signif(pcaVars[i],2),'% variance')
    y_lab=paste0('PC',i+1,': ',signif(pcaVars[i+1],2),'% variance')
    plot(PCs[,i], PCs[,i+1], xlab=x_lab, ylab=y_lab, pch=21, cex=1.2, col='black', bg=study, cex.lab=2, cex.axis=2, xlim=c(min(PCs[,i])-10,max(PCs[,i])+10))
    #textxy(PCs[,i], PCs[,i+1],validation_dataset_metadata$names,cex =.7, offset = .7)
    legend("topright", legend=levels(study), pch=15, col=col,  cex=1, pt.cex=1, title='Study')
}

dev.off()

################ Make table of DE genes

# print results
sig_results=as.data.frame(rr[which(rr$padj<=0.1),])

# order results by LFC
oo=order(sig_results$log2FoldChange)
sig_results=sig_results[oo,]

# make output table
out=data.frame(gene=sig_results$gene, chr=sig_results$chr, Ensembl=sig_results$Ensembl, log2FoldChange=sig_results$log2FoldChange, FDR=sig_results$padj) 

nrow(out) # 1658
length(which(out$log2FoldChange > 0)) # 873 genes overexpressed in samples with distant mets (relative to localized samples)
length(which(out$log2FoldChange < 0)) # 785 genes are underexpressed in samples with distant mets (relative to localized samples)

# How many DE genes are TFs?

out$TF=FALSE

TFs=read.csv("/athena/masonlab/scratch/users/nai2008/Human_TFs_DatabaseExtract_v_1.01.csv")
TFs$Ensembl_ID=as.vector(TFs$Ensembl_ID)
nrow(TFs) # 2765
length(unique(TFs$Ensembl_ID)) # 2765
which(is.na(TFs$Ensembl_ID)) # 0

mm=match(out$Ensembl,TFs$Ensembl_ID)
which(duplicated(out$Ensembl)) # none
which(duplicated(TFs$Ensembl_ID)) # none

out$TF[which(!is.na(mm))]=TRUE
length(which(out$TF==TRUE)) # 178

save(out, file='/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/validation_dataset_DE_genes_byMetStatus.rda')
# load('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/validation_dataset_DE_genes_byMetStatus.rda')

# print table of DE genes
write.csv(out,file="/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_data_tables/validation_dataset_DE_genes_byMetStatus.csv", row.names=FALSE)

################ Perform gene set over-representation analysis (ORA)

library(goseq)
load('/athena/masonlab/scratch/users/nai2008/items_for_goseq_analysis.rda') # gene2cat_GOandKEGG, KEGG_term_names, median_tx_lengths, cat2gene_GO, cat2gene_KEGG

if(length(which(is.na(rr$padj))) != 0) { 
  rr_mod=rr[-which(is.na(rr$padj)),]
} else { rr_mod = rr }

indicator=rep(0, times=nrow(rr_mod))
indicator[which(rr_mod$padj<=0.1)]=1
aa=indicator
names(aa)=rr_mod$Ensembl

mm = match(names(aa), median_tx_lengths$gene_EnsemblID)
bias.data = median_tx_lengths$median_length[mm]
pwf = nullp(aa, 'hg38', 'ensGene', bias.data = bias.data)

pdf('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_pdfs/validation_dataset_goseq_pwf_plot.pdf')
plotPWF(pwf)
dev.off()

GO.KEGG.wall=goseq(pwf,"hg38","ensGene", gene2cat = gene2cat_GOandKEGG, test.cats=c("GO:CC", "GO:BP", "GO:MF", "KEGG"))
GO.KEGG.wall$over_represented_FDR=p.adjust(GO.KEGG.wall$over_represented_pvalue, method="BH")

GO.KEGG.wall$ontology[grep('path:hsa', GO.KEGG.wall$category)]='KEGG'

index = grep('path:hsa', GO.KEGG.wall$category)

for (i in 1:length(index)){
    mm=match(GO.KEGG.wall$category[index[i]], KEGG_term_names$KEGG_ID)
    GO.KEGG.wall$term[index[i]] = KEGG_term_names$KEGG_term[mm]
}

length(which(GO.KEGG.wall$over_represented_FDR<=0.1)) # 13
GO.KEGG.wall_sig = GO.KEGG.wall[which(GO.KEGG.wall$over_represented_FDR<=0.1),]

# Add DE genes in each GO/KEGG category

GO.KEGG.wall_sig_withoutGenes = GO.KEGG.wall_sig

library(AnnotationHub)
hub = AnnotationHub()
dm = query(hub, c("EnsDb", "sapiens", "97"))
edb = dm[["AH73881"]]
genes=as.data.frame(genes(edb))

ens.gene.map = data.frame(gene_id=genes$gene_id, gene_name=genes$gene_name)

length(names(cat2gene_GO)) == length(unique(names(cat2gene_GO))) #TRUE
length(names(cat2gene_KEGG)) == length(unique(names(cat2gene_KEGG))) #TRUE

GO.KEGG.wall_sig$genes_Ensembl=NA
GO.KEGG.wall_sig$genes=NA

for (i in 1:nrow(GO.KEGG.wall_sig)){

  cat=GO.KEGG.wall_sig$category[i]

  if (length(grep('GO',cat)) == 1){

      m.cat=match(cat, names(cat2gene_GO))

      if(is.na(m.cat)){print('error: m.cat does not match (GO)')} else {

          possible_genes=cat2gene_GO[[m.cat]]

          m.genes=match(possible_genes,names(aa))

          if( length(which(!is.na(m.genes)))==0 ){print('error: m.genes are all <NA> (GO)')} else {

              if (length(which(is.na(m.genes)))>0){ possible_genes= possible_genes[-which(is.na(m.genes))] }

              m.genes=match(possible_genes,names(aa))
              subset=aa[m.genes]
              DE_genes=subset[which(subset==1)]
              GO.KEGG.wall_sig$genes_Ensembl[i]=paste(names(DE_genes),collapse=';')

              m.ens=match(names(DE_genes),ens.gene.map$gene_id)
              GO.KEGG.wall_sig$genes[i]=paste(ens.gene.map$gene_name[m.ens],collapse=';')
              }
      }
  } else if (length(grep('path:hsa',cat)) == 1){

      m.cat=match(cat, names(cat2gene_KEGG))

      if(is.na(m.cat)){print('error: m.cat does not match (KEGG)')} else {

          possible_genes=cat2gene_KEGG[[m.cat]]

          m.genes=match(possible_genes,names(aa))

          if(length(which(!is.na(m.genes) == 0))){print('error: m.genes are all <NA> (KEGG)')} else{

              if (length(which(is.na(m.genes)))>0){ possible_genes= possible_genes[-which(is.na(m.genes))] }

              m.genes=match(possible_genes,names(aa))
              subset=aa[m.genes]
              DE_genes=subset[which(subset==1)]
              GO.KEGG.wall_sig$genes_Ensembl[i]=paste(names(DE_genes),collapse=';')

              m.ens=match(names(DE_genes),ens.gene.map$gene_id)
              GO.KEGG.wall_sig$genes[i]=paste(ens.gene.map$gene_name[m.ens],collapse=';')
          }
      }
  }
}

write.csv(GO.KEGG.wall_sig, file = '/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_data_tables/validation_dataset_DEbyMetStatus_ORA.csv')

################ Make TPM plots of validated genes with concordant LFCs b/w the discovery and validation datasets

load('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/validated_genes_DEbyMetStatus_concordantLFC_only.rda')
log2_tpm_plus1=log2(tpm+1)

x=as.vector(ddsClean$Distant_Mets)
x=sub('Y','Distant mets',x)
x=sub('N','No distant mets',x)
x=factor(x, levels=c('No distant mets','Distant mets'))

mycol=as.vector(ddsClean$Distant_Mets)
mycol[which(mycol=="N")]='purple'
mycol[which(mycol=="Y")]='green'

pdf('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_pdfs/validation_dataset_TPM_plots_of_validated_genes_with_concordantLFCs.pdf')

for (i in 1: nrow(validated_genes_concordantLFC_only)){

  par(mar=c(5.1,5.3,4.1,2.1))

  index=which(rr$Ensembl==validated_genes_concordantLFC_only$Ensembl[i])

  zag1=paste0(rr$gene[index], ' (',rr$Ensembl[index],')')
  zag2=as.expression(bquote(log[2]~"FC" == .(signif(rr$log2FoldChange[index],2))))
  zag3=paste0("FDR = ",signif(rr$padj[index],2))

  log2_tpm_plus1_subset=as.vector(log2_tpm_plus1[which(rownames(log2_tpm_plus1)==validated_genes_concordantLFC_only$Ensembl[i]),])

  boxplot(as.vector(log2_tpm_plus1_subset)~x, , xlab='Tumor status', ylab= as.expression(bquote(log[2]~"(TPM+1)")), main=zag1, cex.main=2, cex.lab=2, cex.axis=1.5, outline=FALSE, col='lightgrey', ylim=c( min(as.vector(log2_tpm_plus1_subset)), max(as.vector(log2_tpm_plus1_subset)) ) )
  points(as.vector(log2_tpm_plus1_subset) ~ jitter(as.numeric(x), amount=0.2), pch =21, col='black', bg=mycol, cex=1.4)
  legend(x='topright',legend=c(zag2,zag3), bty='n')

}

dev.off()

################ UCHL1 expression

#UCHL1 is ENSG00000154277

rr[which(rr$Ensembl=='ENSG00000154277'),]
#                  baseMean log2FoldChange     lfcSE      stat    pvalue
#                 <numeric>      <numeric> <numeric> <numeric> <numeric>
# ENSG00000154277   4506.25      -0.127661    0.4153 -0.307394  0.758544
#                      padj         chr         Ensembl        gene
#                 <numeric> <character>     <character> <character>
# ENSG00000154277  0.999621           4 ENSG00000154277       UCHL1


################ Clustering analysis

load('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/validated_genes_DEbyMetStatus_concordantLFC_only.rda')
load('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/validation_dataset_metadata.rda')

# Patient 429 (gave rise to primary tumor sample 29_S111) and patient 492 (gave rise to primary tumor sample 37_S133) had localized disease at the time of sequencing, but developed distant mets to the liver after resection
# add that info to the metadata
validation_dataset_metadata$Loc.to.Met=NA
validation_dataset_metadata$Loc.to.Met[which(validation_dataset_metadata$Distant_Mets=='Y')]=NA
validation_dataset_metadata$Loc.to.Met[which(validation_dataset_metadata$names=='29_S111')]='Yes'
validation_dataset_metadata$Loc.to.Met[which(validation_dataset_metadata$names=='37_S133')]='Yes'

###### Read in the data

source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')
library(DESeq2)
library(tximeta)

# import data
se = tximeta(coldata=validation_dataset_metadata, type = "salmon")
# found matching transcriptome:
# [ Ensembl - Homo sapiens - release 97 ]

# summarize transcript-level quantifications to gene-level
gse = summarizeToGene(se)

# get TPM matrix
tpm = assays(gse)$abundance

#get count matrix
counts=assays(gse)$counts

# make DESeqDataSet object
dds = DESeqDataSet(gse, design = ~ Distant_Mets)

#perform pre-filtering to keep only rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# make a transformed count matrix, using variance stabilizing transformation (VST)
vsd = vst(dds, blind=FALSE)
vst_counts = as.matrix(assay(vsd))

# regress out the batch variable
library(jaffelab)

mod = model.matrix(~Distant_Mets + factor(Dataset), data=as.data.frame(colData(dds))) 
clean_vst_counts=cleaningY(vst_counts, mod, P=2)

library(AnnotationHub)
hub = AnnotationHub()
dm = query(hub, c("EnsDb", "sapiens", "97"))
edb = dm[["AH73881"]]
genes=as.data.frame(genes(edb))

mm_validated_genes_w_concordantLFCs=match(validated_genes_concordantLFC_only$Ensembl, rownames(clean_vst_counts))

if(length(which(is.na(mm_validated_genes_w_concordantLFCs))) != 0) { 
  validated_genes_concordantLFC_only = validated_genes_concordantLFC_only[-which(is.na(mm_validated_genes_w_concordantLFCs)),]
  mm_validated_genes_w_concordantLFCs=match(validated_genes_concordantLFC_only$Ensembl, rownames(clean_vst_counts))
 }

which(duplicated(validated_genes_concordantLFC_only$Ensembl))# none
which(duplicated(rownames(clean_vst_counts)))# none

mm=match(rownames(clean_vst_counts),genes$gene_id)
which(duplicated(rownames(clean_vst_counts)))# none
which(duplicated(genes$gene_id))# none
which(is.na(mm)) # none

rownames(clean_vst_counts) = genes$gene_name[mm]

all(colnames(clean_vst_counts)==validation_dataset_metadata$names) #TRUE

library(pheatmap)
library(RColorBrewer)

annotation_col=data.frame(Distant.Mets=validation_dataset_metadata$Distant_Mets, Loc.to.Met=validation_dataset_metadata$Loc.to.Met)
rownames(annotation_col)=as.vector(validation_dataset_metadata$names)

ann_colors = list(
  Loc.to.Met = c(Yes='blue'),
  Distant.Mets = c(Y='green', N='purple')
  )

pdf('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_pdfs/validation_analysis_heatmaps_pearsonCorrelation.pdf')

pheatmap(clean_vst_counts[mm_validated_genes_w_concordantLFCs,], color=colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(3000), 
  main='Validation Dataset; Pearson correlation',
  clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",
  cluster_rows=TRUE, cluster_cols=TRUE, annotation_col=annotation_col, scale='row', fontsize_col=5, annotation_colors = ann_colors)

dev.off()

tt=clean_vst_counts[mm_validated_genes_w_concordantLFCs,]
shapiro.test(as.vector(tt[5,]))














