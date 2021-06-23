#
###########################################################################################
##### Make metadata tables for the datasets
###########################################################################################

source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

## Chan et al PNET dataset
chan_metadata=read.csv('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/Chan_etal_PNETs_RNAseq/Metadata/master_clinicopathological_metadata_file_SUBSET.csv')

chan_metadata$Sex=sub('F','Female',chan_metadata$Sex)
chan_metadata$Sex=sub('M','Male',chan_metadata$Sex)

colnames(chan_metadata)[1]='names'

table(chan_metadata$Distant_Mets)
#  N  Y
# 14 10

lf=list.files('/athena/masonlab/scratch/users/nai2008/PNET/Chan_etal_NatComm_2018/Salmon_quants',full.names=TRUE, recursive=TRUE)
lf_quant=lf[grep('/quant.sf',lf)]

id=ss(lf_quant,'/',10)

mm=match(chan_metadata$ID_for_Salmon_files,id)
which(duplicated(chan_metadata$ID_for_Salmon_files)) # none
which(duplicated(id)) # none
which(is.na(mm)) # none
which(duplicated(mm)) # none

chan_metadata$files=lf_quant[mm]
chan_metadata$Distant_Mets=factor(chan_metadata$Distant_Mets,levels=c('N','Y'))

save(chan_metadata,file='/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/chan_metadata.rda')

## Scarpa et al PNET dataset
scarpa_metadata=read.csv('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/Scarpa_etal_PNETs_RNAseq/Metadata/master_clinicopathological_metadata_file_SUBSET.csv')

scarpa_metadata$LN_Mets=sub('yes','Y',scarpa_metadata$LN_Mets)
scarpa_metadata$LN_Mets=sub('no','N',scarpa_metadata$LN_Mets)

scarpa_metadata$Distant_Mets=sub('yes','Y',scarpa_metadata$Distant_Mets)
scarpa_metadata$Distant_Mets=sub('no','N',scarpa_metadata$Distant_Mets)

colnames(scarpa_metadata)[1]='names'

table(scarpa_metadata$Distant_Mets)
#  N  Y
# 27  2

lf=list.files('/athena/masonlab/scratch/users/nai2008/PNET/Scarpa_etal_EGAD00001003336_PNETs_RNAseq/Salmon_quants',full.names=TRUE, recursive=TRUE)
lf_quant=lf[grep('/quant.sf',lf)]

id=ss(lf_quant,'/',10)

mm=match(scarpa_metadata$ID_for_Salmon_files,id)
which(duplicated(scarpa_metadata$ID_for_Salmon_files)) # none
which(duplicated(id)) # none
which(is.na(mm)) # none
which(duplicated(mm)) # none

scarpa_metadata$files=lf_quant[mm]
scarpa_metadata$Distant_Mets=factor(scarpa_metadata$Distant_Mets,levels=c('N','Y'))

save(scarpa_metadata,file='/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/scarpa_metadata.rda')

## Scott at al PNET dataset (James Howe is the study PI)
scott_metadata=read.csv('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/Scott_etal_PNETs_RNAseq_JamesHowe_dbGaP/Metadata/master_clinicopathological_metadata_file_SUBSET.csv')

scott_metadata$Sex=sub('F','Female',scott_metadata$Sex)
scott_metadata$Sex=sub('M','Male',scott_metadata$Sex)

scott_metadata$LN_Mets=sub('Yes','Y',scott_metadata$LN_Mets)
scott_metadata$LN_Mets=sub('No','N',scott_metadata$LN_Mets)

scott_metadata$Distant_Mets=sub('yes','Y',scott_metadata$Distant_Mets)
scott_metadata$Distant_Mets=sub('no','N',scott_metadata$Distant_Mets)

colnames(scott_metadata)[1]='names'

table(scott_metadata$Distant_Mets)
 # N  Y
 # 9 11

lf=list.files('/athena/masonlab/scratch/users/nai2008/PNET/Scott_etal_PNETs_RNAseq_JamesHowe_dbGaP/Salmon_quants',full.names=TRUE, recursive=TRUE)
lf_quant=lf[grep('/quant.sf',lf)]

id=ss(lf_quant,'/',10)

mm=match(scott_metadata$ID_for_Salmon_files,id)
which(duplicated(scott_metadata$ID_for_Salmon_files)) # none
which(duplicated(id)) # none
which(is.na(mm)) # none
which(duplicated(mm)) # none

scott_metadata$files=lf_quant[mm]
scott_metadata$Distant_Mets=factor(scott_metadata$Distant_Mets,levels=c('N','Y'))

save(scott_metadata,file='/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/scott_metadata.rda')

## Our WCM PNET dataset
wcm_metadata=read.csv('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/WCM_PNETs_RNAseq/Metadata/master_clinicopathological_metadata_file_SUBSET.csv')

wcm_metadata$Sex=sub('F','Female',wcm_metadata$Sex)
wcm_metadata$Sex=sub('M','Male',wcm_metadata$Sex)

wcm_metadata$LN_Mets=sub('yes','Y',wcm_metadata$LN_Mets)
wcm_metadata$LN_Mets=sub('no','N',wcm_metadata$LN_Mets)

wcm_metadata$Distant_Mets=sub('yes','Y',wcm_metadata$Distant_Mets)
wcm_metadata$Distant_Mets=sub('no','N',wcm_metadata$Distant_Mets)

colnames(wcm_metadata)[1]='names'

table(wcm_metadata$Distant_Mets)
#  N  Y
# 17  5

lf=list.files('/athena/masonlab/scratch/users/nai2008/PNET/WCM_PNET_RNAseq/Salmon_quants',full.names=TRUE, recursive=TRUE)
lf_quant=lf[grep('/quant.sf',lf)]

id=ss(lf_quant,'/',10)

mm=match(wcm_metadata$ID_for_Salmon_files,id)
which(duplicated(wcm_metadata$ID_for_Salmon_files)) # none
which(duplicated(id)) # none
which(is.na(mm)) # none
which(duplicated(mm)) # none

wcm_metadata$files=lf_quant[mm]
wcm_metadata$Distant_Mets=factor(wcm_metadata$Distant_Mets,levels=c('N','Y'))

save(wcm_metadata,file='/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/wcm_metadata.rda')

# Combine WCM w/ Chan to make the discovery dataset; combine Scott w/ Scarpa to make the validation dataset
# (this way we balance out the sample size and power b/w the discovery and validation analyses)

discovery_dataset_metadata=rbind(wcm_metadata, chan_metadata)
discovery_dataset_metadata$Distant_Mets=factor(discovery_dataset_metadata$Distant_Mets,levels=c('N','Y'))
nrow(discovery_dataset_metadata) # 46
table(discovery_dataset_metadata$Distant_Mets)
#  N  Y
# 31 15


validation_dataset_metadata=rbind(scott_metadata, scarpa_metadata)
validation_dataset_metadata$Distant_Mets=factor(validation_dataset_metadata$Distant_Mets,levels=c('N','Y'))
nrow(validation_dataset_metadata) #49
table(validation_dataset_metadata$Distant_Mets)
#  N  Y
# 36 13

save(discovery_dataset_metadata, file='/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/discovery_dataset_metadata.rda')
#load('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/discovery_dataset_metadata.rda')

save(validation_dataset_metadata, file='/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/validation_dataset_metadata.rda')
#load('/athena/masonlab/scratch/users/nai2008/PNET_FinnertyProject/_rdas/validation_dataset_metadata.rda')








