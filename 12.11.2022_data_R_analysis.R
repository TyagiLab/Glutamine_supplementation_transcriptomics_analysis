
library(DESeq2)
library(dplyr)
library(biomaRt)

importSTARgenecounts <- function(data_folder, strandedness) {
  
  # Pulls out files the Gene counts file names 
  files <- list.files(data_folder, pattern = ".*ReadsPerGene.out.tab",full.names = T, recursive = FALSE)
  
  # Makes a sample_names variable using the file names by removing the beginning and end
  sample_names<-sub(data_folder,"",sub(".ReadsPerGene.out.tab", "",files))
  
  # Name the files by their sample names as metadata
  names(files)<-sample_names
  
  # This is a bit messy but does the job
  # It loads in the first values to make a tibble of the appropriate size
  temp_STAR_table<-read.table(files[[1]],col.names = c('ensembl','total.mapped','for.mapped','rev.mapped'))
  
  # Then it loads in all the rest
  for (name in names(files)) {
    if (strandedness == 'reverse')
    {
      temp_STAR_table[name] <- read.table(files[name],col.names = c('ensembl','total.mapped','for.mapped','rev.mapped'))$rev.mapped
      }
    else if (strandedness == 'forward')
    {
      temp_STAR_table[name] <- read.table(files[name],col.names = c('ensembl','total.mapped','for.mapped','rev.mapped'))$for.mapped
    }
    else
    {
      stop('Please put \'reverse\' or \'forward\' for strandedness')
    }
    }
  
  # Then it removes the initialization stuff
  STAR_table <- subset(temp_STAR_table,select=-c(total.mapped,for.mapped,rev.mapped))
  
  # Then remove the N_data
  STAR_table <- STAR_table[-1:-4, ]
  
  # Make the index the gene names and remove it
  rownames(STAR_table) <- STAR_table$ensembl
  STAR_table$ensembl <- NULL
  
  return(STAR_table)
}

tpm <- function (counts_data_frame, gene_lengths_file) {
  # Gene lengths file needs to be from GTFtools
  # Installed GTFtools with `pip install gtftools` 
  # then with this command to calculate effective length for 
  # TPM calculation 
  # `gtftools -l gencode.vM31.primary_assembly.annotation.gene_lengths.txt gencode.vM31.primary_assembly.annotation.gtf`
  
  # Reads in gene lengths
  gene_lengths<-read.table(gene_lengths_file, header = TRUE)
  
  # Adds gene names from rownames to a column
  counts_data_frame$ensembl<-rownames(counts_data_frame)
  # Combines the genes from the gene lengths file and the counts dataframe
  # Note that this will remove genes missing from either
  counts_table_temp_w_lengths <- counts_data_frame %>%
    dplyr::inner_join(gene_lengths[c('gene','median')], by = c('ensembl' = 'gene'))
  
  # Then remove the gene column and divide each row by it's gene length in kb
  # Choosing to use median, shouldn't matter too much since cross gene comparison isn't a thing and shouldn't be done really
  rownames(counts_table_temp_w_lengths)<-counts_table_temp_w_lengths$ensembl
  counts_table_temp_w_lengths$ensembl<-NULL
  FPK_counts<-(counts_table_temp_w_lengths/counts_table_temp_w_lengths$median)*1000
  
  # Remove the gene length merged column
  FPK_counts$median<-NULL
  
  TPM_counts<-as.data.frame((prop.table(as.matrix(FPK_counts),2))*10e6)
  
  return(TPM_counts)
}

# Formatting count data from alignments folder into a dataframe
STAR_table <- importSTARgenecounts("D:/12.11.2022_data/alignments/", 'reverse')

# Using biomart for annotations since it's up to date compared to annotables

ensembl <- useEnsembl(biomart = "genes")
ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl)
gene_info<-getBM(c('ensembl_gene_id','external_gene_name','description'),mart=ensembl)

# Make a temporary variable (which is the same as the output and will be overwritten)
annot_STAR_table<-STAR_table
# Make a column with the ensembl Gene IDs (no version for compatibility)
annot_STAR_table$ensembl<-gsub("\\..*", "", rownames(annot_STAR_table))
# Join the annotables 
annot_STAR_table<- annot_STAR_table %>% 
  dplyr::inner_join(gene_info, by = c("ensembl" = "ensembl_gene_id")) %>% 
  dplyr::select("external_gene_name",'ensembl', everything())

# Saving the gene counts to a .csv. Note that these have not been normalized for reads per sample or for gene length,
write.csv(STAR_table, file = "D:/LanboData/12.11.2022_data/R_analysis_outputs/raw_counts_STAR.csv", row.names = F)
write.csv(annot_STAR_table, file = "D:/LanboData/12.11.2022_data/R_analysis_outputs/annotated_raw_counts_STAR.csv", row.names = F)

# Calculating TPM
STAR_TPM<-tpm(STAR_table,"D:/12.11.2022_data/alignments/gencode.vM31.primary_assembly.annotation.gene_lengths.txt")

# Make a temporary variable (which is the same as the output and will be overwritten)
annot_STAR_TPM<-STAR_TPM
# Make a column with the ensembl Gene IDs (no version for compatibility)
annot_STAR_TPM$ensembl<-gsub("\\..*", "", rownames(annot_STAR_TPM))
# Join the annotables 
annot_STAR_TPM<- annot_STAR_TPM %>% 
  dplyr::inner_join(gene_info, by = c("ensembl" = "ensembl_gene_id")) %>% 
  dplyr::select("external_gene_name",'ensembl', everything())

# Saving TPM to file
write.csv(STAR_TPM, file = "D:/LanboData/12.11.2022_data/R_analysis_outputs/TPM_counts_STAR.csv", row.names = F)
write.csv(annot_STAR_TPM, file = "D:/LanboData/12.11.2022_data/R_analysis_outputs/annotated_TPM_counts_STAR.csv", row.names = F)




# Full conditions table
# This was never used, but was useful to look at haha.
# condition_table <- data.frame(name = colnames(STAR_table),
#            infected=c(rep("uninf", 4),rep("inf", 16),rep("uninf", 4),rep("inf", 8),rep("uninf", 4),rep("inf", 12)),
#            days=c(rep(0, 4),25,rep(55, 4),rep(25, 3),rep(1, 4),rep(0.167, 4),rep(0, 4),rep(1, 4),rep(0.167, 4),rep(0, 4),rep(55, 4),rep(12, 4),rep(25, 4)),
#            glutamine=c(rep('no', 4),rep('yes', 8),rep('no', 36)),
#            inhibitor=c(rep('no', 20),rep('yes', 12),rep('no', 16)))
# 
# rownames(condition_table) <- condition_table$name
# condition_table$name<-NULL


# Make a conditions table for DESeq2 and label each main condition 
# Much more useful just to make categories of samples you'd like to compare

condition_table <- data.frame(name = colnames(STAR_table),
                              condition=c(rep("C0", 4),"C1_25",rep("C1_55", 4),rep("C1_25", 3),rep("C24", 4),rep("C4", 4),rep("R0", 4),rep("R24", 4),rep("R4", 4),rep("U", 4),rep("V1_55", 4),rep("V1_12", 4),rep("V1_25", 4)))
# Formatting the data frame for input
rownames(condition_table) <- condition_table$name
condition_table$name<-NULL

# Make the DESeq DataSet object for DESeq2
dds <- DESeqDataSetFromMatrix(countData = STAR_table,
                              colData = condition_table,
                              design = ~ condition)       # Can get rid of batch by putting `design = ~ condition + batch` if you have batch info.
# Run it
dds <- DESeq(dds)

# Parse the results for the comparisons we are interested in
res_C1_25vsV1_25 <- results(dds, contrast=c("condition","C1_25","V1_25"))
res_C1_55vsV1_55 <- results(dds, contrast=c("condition","C1_55","V1_55"))
res_UvsV1_12 <- results(dds, contrast=c("condition","U","V1_12"))
res_UvsV1_25 <- results(dds, contrast=c("condition","U","V1_25"))
res_UvsV1_55 <- results(dds, contrast=c("condition","U","V1_55"))

res_V1_55vsV1_25 <- results(dds, contrast=c("condition","V1_55","V1_25"))


# Add ensembl gene IDs without version as a column
res_C1_25vsV1_25$ensembl<-gsub("\\..*", "", rownames(res_C1_25vsV1_25))
res_C1_55vsV1_55$ensembl<-gsub("\\..*", "", rownames(res_C1_55vsV1_55))
res_UvsV1_12$ensembl<-gsub("\\..*", "", rownames(res_UvsV1_12))
res_UvsV1_25$ensembl<-gsub("\\..*", "", rownames(res_UvsV1_25))
res_UvsV1_55$ensembl<-gsub("\\..*", "", rownames(res_UvsV1_55))

res_V1_55vsV1_25$ensembl<-gsub("\\..*", "", rownames(res_V1_55vsV1_25))


# Now change the ensembl gene IDs to symbols for readability 
annotated_res_C1_25vsV1_25 <- as.data.frame(res_C1_25vsV1_25) %>% 
  dplyr::arrange(padj) %>% 
  dplyr::inner_join(gene_info, by = c("ensembl" = "ensembl_gene_id"))
annotated_res_C1_55vsV1_55 <- as.data.frame(res_C1_55vsV1_55) %>% 
  dplyr::arrange(padj) %>% 
  dplyr::inner_join(gene_info, by = c("ensembl" = "ensembl_gene_id"))
annotated_res_UvsV1_12 <- as.data.frame(res_UvsV1_12) %>% 
  dplyr::arrange(padj) %>% 
  dplyr::inner_join(gene_info, by = c("ensembl" = "ensembl_gene_id"))
annotated_res_UvsV1_25 <- as.data.frame(res_UvsV1_25) %>% 
  dplyr::arrange(padj) %>% 
  dplyr::inner_join(gene_info, by = c("ensembl" = "ensembl_gene_id"))
annotated_res_UvsV1_55 <- as.data.frame(res_UvsV1_55) %>% 
  dplyr::arrange(padj) %>% 
  dplyr::inner_join(gene_info, by = c("ensembl" = "ensembl_gene_id"))

annotated_res_V1_55vsV1_25 <- as.data.frame(res_V1_55vsV1_25) %>% 
  dplyr::arrange(padj) %>% 
  dplyr::inner_join(gene_info, by = c("ensembl" = "ensembl_gene_id"))

# Save DESeq2 outputs to .csv files
write.csv(annotated_res_C1_25vsV1_25, file = "D:/12.11.2022_data/R_analysis_outputs/C1_25vsV1_25_STAR_DESeq2.csv", row.names = F)
write.csv(annotated_res_C1_55vsV1_55, file = "D:/12.11.2022_data/R_analysis_outputs/C1_55vsV1_55_STAR_DESeq2.csv", row.names = F)
write.csv(annotated_res_UvsV1_12, file = "D:/12.11.2022_data/R_analysis_outputs/UvsV1_12_STAR_DESeq2.csv", row.names = F)
write.csv(annotated_res_UvsV1_25, file = "D:/12.11.2022_data/R_analysis_outputs/UvsV1_25_STAR_DESeq2.csv", row.names = F)
write.csv(annotated_res_UvsV1_55, file = "D:/12.11.2022_data/R_analysis_outputs/UvsV1_55_STAR_DESeq2.csv", row.names = F)

write.csv(annotated_res_V1_55vsV1_25, file = "D:/12.11.2022_data/R_analysis_outputs/V1_55vsV1_25_STAR_DESeq2.csv", row.names = F)

#############################################
# fGSEA analysis using mouse MSigDB files
#############################################

library(fgsea)

fGSEA_on_DESeq2 <- function (gmt_folder, annotated_DESeq2_result) {
  gmt_files <- list.files(gmt_folder, pattern = ".*gmt",full.names = T, recursive = FALSE)
  sample_names<-sub(gmt_folder,"",sub(".gmt", "",files))
  names(gmt_files)<-sample_names
  
  for (file in gmt_files) {
    gmt_pathways <- gmtPathways(file)
    fgseaRes <- fgsea(pathways = gmt_pathways, 
                           stats = annotated_DESeq2_result$external_gene_name,
                           minSize=10,
                           maxSize=500)
  }
}

########################
# For new PCA

STAR_table2<-cbind(STAR_table[37:40],STAR_table[45:48], STAR_table[5:12])
condition_table2 <- data.frame(name = colnames(STAR_table2),
                              condition=c(rep("V1_55", 4),rep("V1_25", 4),"C_25",rep("C_55", 4),rep("C_25", 3)))
# Formatting the data frame for input
rownames(condition_table2) <- condition_table2$name
condition_table2$name<-NULL

# Make the DESeq DataSet object for DESeq2
dds <- DESeqDataSetFromMatrix(countData = STAR_table2,
                              colData = condition_table2,
                              design = ~ condition)       # Can get rid of batch by putting `design = ~ condition + batch` if you have batch info.
# Run it
dds <- DESeq(dds)

# VSD transform and PCA
vsd <- vst(dds, blind=FALSE)

plotPCA(vsd, intgroup=c("condition", "condition"))

## Getting the data for plotting with other programs 
#pca_data<-plotPCA(vsd, intgroup=c("condition", "condition"))

write.csv(STAR_table2, file = "D:/12.11.2022_data/R_analysis_outputs/25_55_w_and_wo_glut_for_pca_data.csv", row.names = F)

##################################
# Another PCA

STAR_table2<-STAR_table[33:48]
condition_table2 <- data.frame(name = colnames(STAR_table2),
                               condition=c(rep("U", 4),rep("V1_55", 4),rep("V1_12", 4),rep("V1_25", 4)))
# Formatting the data frame for input
rownames(condition_table2) <- condition_table2$name
condition_table2$name<-NULL

# Make the DESeq DataSet object for DESeq2
dds <- DESeqDataSetFromMatrix(countData = STAR_table2,
                              colData = condition_table2,
                              design = ~ condition)       # Can get rid of batch by putting `design = ~ condition + batch` if you have batch info.
# Run it
dds <- DESeq(dds)

# VSD transform and PCA
vsd <- vst(dds, blind=FALSE)

plotPCA(vsd, intgroup=c("condition", "condition"))

pca_data<-plotPCA(vsd, intgroup=c("condition", "condition"))

pca_data

write.csv(STAR_table2, file = "D:/12.11.2022_data/R_analysis_outputs/25_55_w_and_wo_glut_for_pca_data.csv", row.names = F)


