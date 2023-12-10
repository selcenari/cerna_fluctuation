#Required Packages
library(tidyverse)
library(tidygraph)
library(igraph)
library(Matrix)
library(SPONGE)

### This file just an example of sparse partial correlation workflow for a particular tissue (BRCA). 

#miRNA expression data of All normal tissues:
readRDS("TCGA_data_normal/all_tcga_normal_mirna.RDS")-> all_tcga_normal_mirna

#Gene expression data of All normal tissues:
readRDS("TCGA_data_normal/all_tcga_normal.RDS")-> all_tcga_normal

#miRNA:ceRNA pairs from ENCORI
readRDS("TCGA_data_normal/starbase_pairs.RDS")-> starbase_pairs

#Gene information data including descriptors and several annotation:
readRDS("TCGA_data_normal/data_info_gen.RDS")-> data_info_gen

#miRNA information data including identifier and transcript info:
readRDS("TCGA_data_normal/data_info_mirna.RDS")-> data_info_mirna



# selection of sample to be analyzed:
projects <-"TCGA-BRCA"

# Detection of common samples in miRNA and ceRNA expression data for specified tissue:
intersect(
  (data_info_gen %>% 
     filter(project==projects) %>% 
     dplyr::select(cases, project, cases.submitter_id, sample.submitter_id))$sample.submitter_id,
  (data_info_mirna %>% 
     filter(project==projects) %>% 
     dplyr::select(cases, project, cases.submitter_id, sample.submitter_id))$sample.submitter_id
)->project_samples

# Selection of miRNA expression data for specified tissue:
all_tcga_normal_mirna %>% 
  select(mirbase_id, project_samples) -> project_mirna

as.data.frame(project_mirna)-> project_mirna

rownames(project_mirna) <- project_mirna$mirbase_id

project_mirna <- dplyr::select(project_mirna, -mirbase_id)


# Selection of gene expression data for specified tissue:
all_tcga_normal %>% 
  select(project_samples) -> project_gene

as.matrix(t(project_gene))->project_gene
as.matrix(t(project_mirna))->project_mirna



# Selection of miRNA:ceRNA pairs and conversion to binary interaction data.
starbase_pairs%>%
  filter(geneID %in% colnames(project_gene), miRNAid %in% colnames(project_mirna))%>%
  mutate(pair = 1)%>%
  pivot_wider(names_from = miRNAid, values_from=pair, values_fill=0)->starbase_mat


column_to_rownames(starbase_mat, var = "geneID")-> starbase_mat

#Transformation to matrice
as.matrix(starbase_mat)->starbase_mat

#IMPORTANT: Making dimensions equal in all input sets
project_mirna[,colnames(starbase_mat)]->project_mirna

project_gene[,rownames(starbase_mat)]->project_gene


#Run sponge analysis
genes_miRNA_log <- sponge_gene_miRNA_interaction_filter(
  gene_expr = project_gene,
  mir_expr = project_mirna,
  mir_predicted_targets = starbase_mat)

ceRNA_interactions_log <- sponge(gene_expr = project_gene,
                                 mir_expr = project_mirna,
                                 mir_interactions = genes_miRNA_log)

mscor_null_model_log <- sponge_build_null_model(number_of_samples = nrow(project_gene))


ceRNA_interactions_sign_log <- sponge_compute_p_values(sponge_result = ceRNA_interactions_log,
                                                       null_model = mscor_null_model_log)

saveRDS(ceRNA_interactions_sign_log, "TCGA_BRCA_ceRNA_interactions_starbase.RDS")
