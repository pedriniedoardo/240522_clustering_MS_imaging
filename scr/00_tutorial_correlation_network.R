# libraries ---------------------------------------------------------------
library(tidyverse)
library(corrplot)
library(ComplexHeatmap)
library(circlize)
library(cluster)
library(factoextra)
library(magick)

# library(future)
# 
# # plan before planning
# plan()
# # change the planning
# plan("multisession", workers = 16)
# # confirm the planning
# plan()

# read in the data --------------------------------------------------------
# read in the list of metabolic genes
df <- read_csv("../../data/Metscape Metabolic enzyme list.csv")

# rename the columns
df_test <- df %>% 
  dplyr::rename(entrez = `Entrez ID Metabolic Genes/Enzymes Metscape`,
                kegg_enzyme_id = `KEGG Enzyme  code`)

# conver the gene names to gene symbols to pull the data from the MOFA object
ah <- AnnotationHub()
#query the database to identify the version of interest
query(ah, pattern = c("Homo Sapiens", "EnsDb"))

# copy the specific code ofr the database of interest
edb <- ah[["AH109606"]]
columns(edb)

# use the mapids function to pull the annotation of interest
df_tot <- 
  df_test %>%
  mutate(entrez = as.character(entrez)) %>% 
  mutate(biotype = mapIds(edb,
                          keys = entrez,
                          column = c("GENEBIOTYPE"),
                          keytype = "ENTREZID",
                          multiVals = "first")
  ) %>%
  mutate(symbol = mapIds(edb,
                         keys = entrez,
                         column = c("SYMBOL"),
                         keytype = "ENTREZID",
                         multiVals = "first")) %>%
  mutate(ensembl = mapIds(edb,
                          keys = entrez,
                          column = c("GENEID"),
                          keytype = "ENTREZID",
                          multiVals = "first")) 

# not all the entrez ids have been converted to gene symbols check the missing ones
df_tot %>% 
  dplyr::filter(is.na(symbol)) %>% 
  print(n=40)
# I have checked some of them and the are retired gene symbols

# correlation analysis metabolites ----------------------------------------
# pull the scaled expression values
MOFAobject <- readRDS(file = "../../out/object/new/MOFA2_Sample_adj_top_featureScale_final.rds")

df_exp <- get_data(MOFAobject)

# pull the metabolic genes from the rna dataset
df_rna <- df_exp$mRNA$group1 %>% 
  data.frame() %>% 
  rownames_to_column() %>% 
  dplyr::filter(rowname %in% df_tot$symbol)

# pull all the metabolites
df_met <- df_exp$metabolic$group1 %>% 
  data.frame() %>% 
  rownames_to_column()

# run samsple wise correlation
test01 <- t(df_exp$metabolic$group1)
test02 <- t(df_exp$mRNA$group1)

test <- cor(test01,test02)
saveRDS(test,file = "../../out/object/new/correlation_metabolites_genes.rds")
write_tsv(test %>% 
            data.frame() %>% 
            rownames_to_column("rowname"),file = "../../out/table/new/correlation_metabolites_genes.tsv")

# check the resulting matrix
dim(test)
test[1:5,1:5]

# Perform hierarchical clustering on the correlation matrix
# dist_matrix <- as.dist(1 - test)  # Convert correlation matrix to distance matrix
# hc_result <- hclust(dist_matrix)  # Hierarchical clustering
dist_matrix <- dist(test)
hc_result <- hclust(dist_matrix)

# -------------------------------------------------------------------------
# Determine the number of modules (clusters)
# reference: http://www.sthda.com/english/articles/29-cluster-validation-essentials/96-determiningthe-optimal-number-of-clusters-3-must-know-methods/
# Elbow method 
fviz_nbclust(test, kmeans, method = "wss",k.max = 20) +
  # geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")

# Silhouette method
# fviz_nbclust(test, kmeans, method = "silhouette")+   labs(subtitle = "Silhouette method")

# Gap statistic # nboot = 50 to keep the function speedy.  # recommended value: nboot= 500 for your analysis. # Use verbose = FALSE to hide computing progression. 
set.seed(123)
fviz_nbclust(test, kmeans, nstart = 25,  method = "gap_stat", nboot = 20,k.max = 30)+   labs(subtitle = "Gap statistic method")
# -------------------------------------------------------------------------
num_modules <- 7

# Cut the dendrogram to obtain modules
module_labels <- cutree(hc_result, k = num_modules)
module_labels_fix <- paste0("m_",module_labels)
# Convert module labels to colors
# module_colors <- colorRamp2(seq(min(module_labels), max(module_labels)),
#                             colors = c("red", "green", "blue"))(module_labels)

# build the annotation object 
# column_ha <- HeatmapAnnotation(treat = sample_ordered, 
#                                col = list(treat = c("VEH" = "green", "BMP9" = "gray")))  
column_ha <- HeatmapAnnotation(module = module_labels_fix)

# Create a heatmap with module assignments
hm <- Heatmap(t(test),
              name = "Correlation",
              col = colorRamp2(c(-1, 0, 1), colors = c("blue", "white", "red")),
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 8),
              column_names_side = "top",
              column_names_gp = gpar(fontsize = 8),
              row_title = "mRNA",
              column_title = "Met",
              # clustering_distance_rows = dist_matrix2,
              clustering_distance_columns = dist_matrix,
              # clustering_method_rows = "none",
              # clustering_method_columns = "none",
              row_dend_reorder = FALSE,
              column_dend_reorder = FALSE,
              row_title_gp = gpar(fontsize = 10, fontface = "bold"),
              column_title_gp = gpar(fontsize = 10, fontface = "bold"),
              # row_annotation = anno_text(labels = module_labels,
              #                            col = module_colors,
              #                            fontface = "bold",
              #                            size = 8),
              # row_annotation_side = "left",
              # row_annotation_width = unit(1, "cm"),
              show_column_names = FALSE,
              show_row_names = FALSE,
              top_annotation = column_ha)

pdf("../../out/image/new/05_correlation_metabolites_genes.pdf",height = 5,width = 7)
draw(hm)
dev.off()

# read in the manual annotation
LUT_metabol <- read_tsv("../../out/table/new/LUT_metabolites_tot_manual.tsv") %>% 
  mutate(feature = paste0(rowname))

# pull the features belonging to the same module
df_modules <- data.frame(feature = names(module_labels),module = module_labels_fix) %>% 
  left_join(LUT_metabol,by = "feature") %>% 
  arrange(module)

# save the table with the module annotation
df_modules %>% 
  write_tsv("../../out/table/new/05_correlation_metabolites_modules_genes.tsv")

# save the correlation matrix
saveRDS(test,"../../out/object/new/05_correlation_metabolites_modules_genes.rds")

# -------------------------------------------------------------------------
# try to run enrichment analysis per module of metabolites
# read in the annotation
kegg_metabolites <- readRDS("../../out/object/kegg_t2m_multiGSEA_HMDB.rds")

# build the reference annotation
TERM2GENE_HMDB <- kegg_metabolites %>% 
  data.frame() %>% 
  rownames_to_column("term") %>% 
  pivot_longer(names_to = "feature",values_to = "presence",-term) %>% 
  dplyr::filter(presence == 1) 

# split the list of metabolites into modules
list_modules <- df_modules %>% 
  split(f=.$module)

df_ORA <- lapply(list_modules, function(df_module){
  res <- enricher(df_module$ID, TERM2GENE=TERM2GENE_HMDB,minGSSize = 2)
  return(summary(res))
}) %>% 
  bind_rows(.id = "module_id")

df_ORA %>%
  write_csv("../../out/table/new/05_correlation_metabolites_modules_genes_ORA.tsv")