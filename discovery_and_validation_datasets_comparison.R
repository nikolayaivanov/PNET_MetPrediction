#
######################## Load data

# Discovery dataset
load('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/discovery_dataset_DE_genes_byMetStatus.rda')
DE_genes_discovery=out
dd=nrow(DE_genes_discovery)
dd # 403
dd_tf=length(which(DE_genes_discovery$TF==TRUE))
dd_tf # 26
DE_genes_discovery_TFsOnly=DE_genes_discovery[which(DE_genes_discovery$TF==TRUE),]
which(duplicated(DE_genes_discovery$Ensembl)) # none

# Validation dataset
load('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/validation_dataset_DE_genes_byMetStatus.rda')
DE_genes_validation=out
vv=nrow(DE_genes_validation)
vv # 1658
vv_tf=length(which(DE_genes_validation$TF==TRUE))
vv_tf # 178
DE_genes_validation_TFsOnly=DE_genes_validation[which(DE_genes_validation$TF==TRUE),]
which(duplicated(DE_genes_validation$Ensembl)) # none

## Overlaps b/w DE genes
mm=match(DE_genes_discovery$Ensembl, DE_genes_validation$Ensembl)
oo=length(which(!is.na(mm)))
oo # 34

## Overlaps b/w DE TFs
mm=match(DE_genes_discovery_TFsOnly$Ensembl, DE_genes_validation_TFsOnly$Ensembl)
oo_tf=length(which(!is.na(mm)))
oo_tf # 3

DE_genes_discovery_TFsOnly[which(!is.na(mm)),]
# MYT1L, IRX3, STAT4

######################## Make Venn diagram comparing discovery and validation analyses

library(VennDiagram)

pdf('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_pdfs/DiscoveryAndValidationDatasetsComparison_VennDiagram_all_genes.pdf')

overridePairwise=TRUE

draw.pairwise.venn(area1=dd , area2=vv, cross.area=oo, category = c('Discovery\nAnalysis','Validation\nAnalysis'),
euler.d = TRUE,
scaled = TRUE,
inverted = FALSE,
ext.text = TRUE,
ext.percent = c(0.05, 0.05, 0.005),
lwd = rep(2, 2),
lty = rep("solid", 2),
col = rep("black", 2),
fill = c('blue','green'),
alpha = rep(0.5, 2),
label.col = rep("black", 3),
cex = c(1.5, 1.3, 1.5),
fontface = rep("plain", 3),
fontfamily = rep("sans", 3),
cat.pos = c(-180, 180),
cat.dist = rep(0.05, 2),
cat.cex = rep(2, 2),
cat.col = rep("black", 2),
cat.fontface = rep("bold", 2),
cat.fontfamily = rep("sans", 2),
cat.just = rep(list(c(0.5, 0.5)), 2),
cat.default.pos = "outer",
cat.prompts = TRUE,
ext.pos = rep(0, 2),
ext.dist = rep(0, 2),
ext.line.lty = "solid",
ext.length = rep(0.95, 2),
ext.line.lwd = 1,
rotation.degree = 0,
rotation.centre = c(0.5, 0.5),
ind = TRUE,
sep.dist = 0.05,
offset = 0, cex.prop = NULL,
print.mode = "raw",
sigdigs = 3)

dev.off()

pdf('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_pdfs/DiscoveryAndValidationDatasetsComparison_VennDiagram_TFs_only.pdf')

overridePairwise=TRUE

draw.pairwise.venn(area1=dd_tf , area2=vv_tf, cross.area=oo_tf, category = c('Discovery\nAnalysis','Validation\nAnalysis'),
euler.d = TRUE,
scaled = TRUE,
inverted = FALSE,
ext.text = TRUE,
ext.percent = c(0.05, 0.05, 0.005),
lwd = rep(2, 2),
lty = rep("solid", 2),
col = rep("black", 2),
fill = c('blue','green'),
alpha = rep(0.5, 2),
label.col = rep("black", 3),
cex = c(1.5, 1.3, 1.5),
fontface = rep("plain", 3),
fontfamily = rep("sans", 3),
cat.pos = c(-180, 180),
cat.dist = rep(0.05, 2),
cat.cex = rep(2, 2),
cat.col = rep("black", 2),
cat.fontface = rep("bold", 2),
cat.fontfamily = rep("sans", 2),
cat.just = rep(list(c(0.5, 0.5)), 2),
cat.default.pos = "outer",
cat.prompts = TRUE,
ext.pos = rep(0, 2),
ext.dist = rep(0, 2),
ext.line.lty = "solid",
ext.length = rep(0.95, 2),
ext.line.lwd = 1,
rotation.degree = 0,
rotation.centre = c(0.5, 0.5),
ind = TRUE,
sep.dist = 0.05,
offset = 0, cex.prop = NULL,
print.mode = "raw",
sigdigs = 3)

dev.off()

######################## Make output tables

## Make a table of genes DE by sex in both the validation and discovery PNET cohorts
mm=match(DE_genes_discovery$Ensembl, DE_genes_validation$Ensembl)
validated_genes=DE_genes_discovery[which(!is.na(mm)),c(1:3,6)]
nrow(validated_genes) # 34

m=match(validated_genes$Ensembl, DE_genes_discovery$Ensemb)
validated_genes$log2FoldChange_discovery = DE_genes_discovery$log2FoldChange[m]
validated_genes$FDR_discovery = DE_genes_discovery$FDR[m]

m=match(validated_genes$Ensembl, DE_genes_validation$Ensembl)
validated_genes$log2FoldChange_validation = DE_genes_validation$log2FoldChange[m]
validated_genes$FDR_validation = DE_genes_validation$FDR[m]

validated_genes=validated_genes[order(validated_genes$log2FoldChange_discovery, decreasing=TRUE),]

validated_genes_concordantLFC_only=validated_genes[which(sign(validated_genes$log2FoldChange_discovery)==sign(validated_genes$log2FoldChange_validation)),]
nrow(validated_genes_concordantLFC_only) # 20

write.csv(validated_genes_concordantLFC_only, file='/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_data_tables/validated_genes_DEbyMetStatus_concordantLFC_only.csv', row.names=FALSE)
save(validated_genes_concordantLFC_only, file='/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/validated_genes_DEbyMetStatus_concordantLFC_only.rda')

######################## Volcano plots

library(AnnotationHub)
hub = AnnotationHub()
dm = query(hub, c("EnsDb", "sapiens", "97"))
edb = dm[["AH73881"]]
genes=as.data.frame(genes(edb))

source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')
library(DESeq2)

load('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/discovery_dataset_DESeq2_DEbyMetStatus_BUNDLE.rda') # se, gse, tpm, counts, ddsClean, rr, vsd
rr_discovery=rr
if(length(which(is.na(rr_discovery$padj))) != 0 ) { rr_discovery=rr_discovery[-which(is.na(rr_discovery$padj)),] }

load('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/validation_dataset_DESeq2_DEbyMetStatus_BUNDLE.rda') # se, gse, tpm, counts, ddsClean, rr, vsd
rr_validation=rr
if(length(which(is.na(rr_validation$padj))) != 0 ){ rr_validation=rr_validation[-which(is.na(rr_validation$padj)),] }

load('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/validated_genes_DEbyMetStatus_concordantLFC_only.rda')

library(EnhancedVolcano)

pdf('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_pdfs/DE_analysis_byMetStatus_DiscoveryAndValidationDatasetsComparison_VolcanoPlots.pdf')

a=as.expression(bquote("NS"~"(FDR">='0.1)'))
b=as.expression(bquote("|"~log[2]~"FC"~"|"~">"~"1"))
c=as.expression(bquote('FDR'<'0.1'))
d=as.expression(bquote("|"~log[2]~"FC"~"|"~">"~"1"~"&"~'FDR'<'0.1'))

# discovery analysis

EnhancedVolcano(toptable = as.data.frame(rr_discovery[,c(2,6)]),
	lab = rr_discovery$gene,
	x = 'log2FoldChange',
	y = 'padj',
	labSize = 3,
	xlab = bquote(~Log[2] ~ "Fold Change"),
	ylab = bquote(~-Log[10]~italic(FDR)),
	pCutoff = 0.1,
	FCcutoff=1,
	title='Metastatic vs. localized PNETs',
	subtitle='Discovery analysis',
	caption=NULL,
	legendLabels=c(a,b,c,d),
	legendPosition="none",
	border='full',
	col=c('grey30', 'forestgreen', 'royalblue', 'red2'),
	colAlpha=1/2,
	#xlim=c(XX,XX),
	ylim=c(min(-log10(rr_discovery$padj)),max(-log10(rr_discovery$padj))),
	shape=19,
	legendLabSize=8,
	legendIconSize=15,
	#selectLab=gene_symbol,
	#labCol=color,
	labFace='bold',
	# labhjust = 0,
	# labvjust = -1,
	# drawConnectors=TRUE
	)

EnhancedVolcano(toptable = as.data.frame(rr_discovery[,c(2,6)]),
	lab = rr_discovery$gene,
	x = 'log2FoldChange',
	y = 'padj',
	labSize = 3,
	xlab = bquote(~Log[2] ~ "Fold Change"),
	ylab = bquote(~-Log[10]~italic(FDR)),
	pCutoff = 0.1,
	FCcutoff=1,
	title='Metastatic vs. localized PNETs',
	subtitle='Discovery analysis',
	caption=NULL,
	legendLabels=c(a,b,c,d),
	legendPosition="none",
	border='full',
	col=c('grey30', 'forestgreen', 'royalblue', 'red2'),
	colAlpha=1/2,
	#xlim=c(XX,XX),
	ylim=c(min(-log10(rr_discovery$padj)),max(-log10(rr_discovery$padj))),
	shape=19,
	legendLabSize=8,
	legendIconSize=15,
	selectLab=validated_genes_concordantLFC_only$gene,
	#labCol=validated_genes_concordantLFC_only$label_color,
	labFace='bold',
	# labhjust = 0,
	# labvjust = -1,
	# drawConnectors=TRUE
	)

# validation analysis

EnhancedVolcano(toptable = as.data.frame(rr_validation[,c(2,6)]),
	lab = rr_validation$gene,
	x = 'log2FoldChange',
	y = 'padj',
	labSize = 3,
	xlab = bquote(~Log[2] ~ "Fold Change"),
	ylab = bquote(~-Log[10]~italic(FDR)),
	pCutoff = 0.1,
	FCcutoff=1,
	title='Metastatic vs. localized PNETs',
	subtitle='Validation analysis',
	caption=NULL,
	legendLabels=c(a,b,c,d),
	legendPosition="none",
	border='full',
	col=c('grey30', 'forestgreen', 'royalblue', 'red2'),
	colAlpha=1/2,
	#xlim=c(XX,XX),
	ylim=c(min(-log10(rr_validation$padj)),max(-log10(rr_validation$padj))),
	shape=19,
	legendLabSize=8,
	legendIconSize=15,
	#selectLab=gene_symbol,
	#labCol=color,
	labFace='bold',
	# labhjust = 0,
	# labvjust = -1,
	# drawConnectors=TRUE
	)

EnhancedVolcano(toptable = as.data.frame(rr_validation[,c(2,6)]),
	lab = rr_validation$gene,
	x = 'log2FoldChange',
	y = 'padj',
	labSize = 3,
	xlab = bquote(~Log[2] ~ "Fold Change"),
	ylab = bquote(~-Log[10]~italic(FDR)),
	pCutoff = 0.1,
	FCcutoff=1,
	title='Metastatic vs. localized PNETs',
	subtitle='Validation analysis',
	caption=NULL,
	legendLabels=c(a,b,c,d),
	legendPosition="none",
	border='full',
	col=c('grey30', 'forestgreen', 'royalblue', 'red2'),
	colAlpha=1/2,
	#xlim=c(XX,XX),
	ylim=c(min(-log10(rr_validation$padj)),max(-log10(rr_validation$padj))),
	shape=19,
	legendLabSize=8,
	legendIconSize=15,
	selectLab=validated_genes_concordantLFC_only$gene,
	#labCol=validated_genes_concordantLFC_only$label_color,
	labFace='bold',
	# labhjust = 0,
	# labvjust = -1,
	# drawConnectors=TRUE
	)

#plot the legend

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
a=as.expression(bquote("NS"~"(FDR">='0.1)'))
b=as.expression(bquote("|"~log[2]~"Fold"~"Change"~"|"~">"~"1"))
c=as.expression(bquote('FDR'<'0.1'))
d=as.expression(bquote("|"~log[2]~"Fold"~"Change"~"|"~">"~"1"~"&"~'FDR'<'0.1'))
title=expression(bold('Legend'))
legend("center", legend=c(a,b,c,d), pch=19, col=c("grey30", "forestgreen", "royalblue", "red2"),  cex=1.2, pt.cex=2, title=title)


dev.off()













