library(fgsea)
library(tidyverse)
library(stringr)

run_fgsea_on_DESeq2_output_file <- function(filename,pathway_filename) {
  library(fgsea)
  library(tidyverse)

  res <- read_csv(filename)
  
  res2 <- res %>% # This will select only name and stat, remove NA's, make sure they're distinct by combining and averaging them if multiple
    dplyr::select(external_gene_name, stat) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(external_gene_name) %>% 
    summarize(stat=mean(stat))
  
  ranks <- deframe(res2)
  
  pathways <- gmtPathways(pathway_filename)
  
  fgseaRes <- fgsea(pathways=pathways, stats=ranks)
  
  return(fgseaRes)
}

comparison_files <- list.files("D:/12.11.2022_data/R_analysis_outputs/", pattern = ".*DESeq2.csv",full.names = T, recursive = FALSE)

pathway_gmt_files <- list.files("D:/12.11.2022_data/R_analysis_outputs/GSEA_analysis/MSigDB_gmt_files/", pattern = ".*gmt",full.names = T, recursive = FALSE)

for (comparison_file in comparison_files) {
  
  for (gmt_file in pathway_gmt_files) {

    fgseaRes <- run_fgsea_on_DESeq2_output_file(comparison_file,gmt_file)
    
    comparison_prefix <- sub('\\.csv$', '', basename(comparison_file))
    gmt_prefix <- sub('\\..*', '', basename(gmt_file))
    
    gsea_output_filename <- paste0(dirname(comparison_file),'/GSEA_analysis/GSEA_results/',gmt_prefix,'.',comparison_prefix,"_GSEA_result.csv")
    gsea_image_filename <- paste0(dirname(comparison_file),'/GSEA_analysis/GSEA_results/',gmt_prefix,'.',comparison_prefix,"_GSEA_top.png")
    
    ## This makes a new matrix which converts the leading edge vectors to text
    formated_matrix<-fgseaRes[order(padj), ]
    formated_matrix[, 8] <- apply(formated_matrix[, 8, drop = FALSE], 1, function(x) paste(x, collapse = ","))
    
    ## Then saves the formatted matrix 
    write_csv(formated_matrix, gsea_output_filename)
    
    most_sig_up <- fgseaRes[order(padj), ] %>%
      filter(NES > 0) %>%
      head(n = 5)
    
    most_sig_down <- fgseaRes[order(padj), ] %>%
      filter(NES < 0) %>%
      head(n = 5)
    
    sig_up_and_down <- rbind(most_sig_up, most_sig_down)
    
    
    ggplot(sig_up_and_down, aes(reorder(pathway, NES), NES)) +
      geom_col(aes(fill=padj<0.005)) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
           title="Hallmark pathways NES from GSEA") + 
      theme_minimal() +
      scale_x_discrete(labels = function(x) str_wrap(gsub("_", " ", x, fixed=TRUE), width = 30)) 
    
    ggsave(gsea_image_filename, dpi = 'retina', bg = "white")
  }
  
}

# Rerun just for the one comparison Lanbo asked me
comparison_file<- "D:/12.11.2022_data/R_analysis_outputs/V1_55vsV1_25_STAR_DESeq2.csv"

pathway_gmt_files <- list.files("D:/12.11.2022_data/R_analysis_outputs/GSEA_analysis/MSigDB_gmt_files/", pattern = ".*gmt",full.names = T, recursive = FALSE)


for (gmt_file in pathway_gmt_files) {
  
  fgseaRes <- run_fgsea_on_DESeq2_output_file(comparison_file,gmt_file)
  
  comparison_prefix <- sub('\\.csv$', '', basename(comparison_file))
  gmt_prefix <- sub('\\..*', '', basename(gmt_file))
  
  gsea_output_filename <- paste0(dirname(comparison_file),'/GSEA_analysis/GSEA_results/',gmt_prefix,'.',comparison_prefix,"_GSEA_result.csv")
  gsea_image_filename <- paste0(dirname(comparison_file),'/GSEA_analysis/GSEA_results/',gmt_prefix,'.',comparison_prefix,"_GSEA_top.png")
  
  ## This makes a new matrix which converts the leading edge vectors to text
  formated_matrix<-fgseaRes[order(padj), ]
  formated_matrix[, 8] <- apply(formated_matrix[, 8, drop = FALSE], 1, function(x) paste(x, collapse = ","))
  
  ## Then saves the formatted matrix 
  write_csv(formated_matrix, gsea_output_filename)
  
  most_sig_up <- fgseaRes[order(padj), ] %>%
    filter(NES > 0) %>%
    head(n = 5)
  
  most_sig_down <- fgseaRes[order(padj), ] %>%
    filter(NES < 0) %>%
    head(n = 5)
  
  sig_up_and_down <- rbind(most_sig_up, most_sig_down)
  
  
  ggplot(sig_up_and_down, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=padj<0.005)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="Hallmark pathways NES from GSEA") + 
    theme_minimal() +
    scale_x_discrete(labels = function(x) str_wrap(gsub("_", " ", x, fixed=TRUE), width = 30)) 
  
  ggsave(gsea_image_filename, dpi = 'retina', bg = "white")
}
