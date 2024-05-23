# AIM ---------------------------------------------------------------------
# explore the dataset


# libraries ---------------------------------------------------------------
library(tidyverse)
library(mixOmics)
library(GGally)
library(patchwork)

# read in the data --------------------------------------------------------
df <- read_csv("../../data/Brussels_LesionAge_Database.csv")

# move the id of the samples to the rownames
df_tot <- df %>%
  mutate(rowname = paste0(id,"_",subject)) %>%
  dplyr::select(-c(id,subject)) %>%
  column_to_rownames("rowname") %>%
  # select only the continous variables
  dplyr::select(-c(sex,diagnosis,diagnosis_binary,PRLStatus,LongT1Cutoff,ConfLevel)) %>%
  # lesion age days and years are redoundant keep only one
  dplyr::select(-c(LesionAge_Years))

# filter out the non complete cases
df_tot_noNA <- df_tot %>%
  filter(complete.cases(df_tot))

dim(df_tot_noNA)

# define the LUT of the samples
LUT_sample <- df %>%
  mutate(rowname = paste0(id,"_",subject)) %>%
  dplyr::select(-c(id,subject)) %>%
  column_to_rownames("rowname") %>%
  # select only the continous variables
  dplyr::select(sex,diagnosis,diagnosis_binary,PRLStatus,LongT1Cutoff,ConfLevel) %>%
  mutate(sex = case_when(sex == 0 ~ "Male",
                         sex == 1 ~ "Female"),
         diagnosis = case_when(diagnosis == 0 ~ "RRMS",
                               diagnosis == 1 ~ "SPMS",
                               diagnosis == 2 ~ "PPMS"),
         diagnosis_binary = case_when(diagnosis_binary == 0 ~ "RRMS",
                                      diagnosis_binary == 1 ~ "PMS"),
         PRLStatus = case_when(PRLStatus == 0 ~ "non-PRL",
                               PRLStatus == 1 ~ "PRL"),
         LongT1Cutoff = case_when(LongT1Cutoff == 0 ~ "short-T1",
                                  LongT1Cutoff == 1 ~ "Long_T1"),
         ConfLevel = case_when(ConfLevel == 0 ~ "No PRL",
                               ConfLevel == 2 ~ "Doubtful",
                               ConfLevel == 3 ~ "Neutral",
                               ConfLevel == 4 ~ "Certain",
                               ConfLevel == 5 ~ "Very Certain PRL")) %>%
  rownames_to_column("sample")

write_tsv(LUT_sample,"../../out/table/LUT_sample.tsv")

# EDA ---------------------------------------------------------------------
lapply(LUT_sample[,-1], as.factor) %>%
  data.frame() %>%
  summary()

# check the distribution of the varibales after center and scaling
df_tot_noNA_scale <- scale(df_tot_noNA, center=T, scale=T)

# show the distribution of the variables
df_tot_noNA_scale %>%
  data.frame() %>%
  rownames_to_column() %>%
  pivot_longer(names_to = "variable",values_to = "scaled",-rowname) %>%
  ggplot(aes(x=scaled))+geom_density()+facet_wrap(~variable)+
  theme_bw() +
  theme(strip.background = element_blank())

df_tot_noNA_scale %>%
  data.frame() %>%
  rownames_to_column() %>%
  pivot_longer(names_to = "variable",values_to = "scaled",-rowname) %>%
  ggplot(aes(x=scaled))+geom_density(aes(col=variable))+
  theme_bw() +
  theme(strip.background = element_blank())

# see correlation between all the variables
df_tot_noNA_scale %>%
  data.frame() %>%
  ggpairs(upper = "blank")+
  theme_bw() +
  theme(strip.background = element_blank())

# scree plot
tune_pca <- tune.pca(df_tot_noNA, scale = TRUE)
plot(tune_pca)

# run PCA pull all the components
pca_all <- pca(df_tot_noNA, center = TRUE, scale = TRUE,ncomp = ncol(df_tot_noNA))

# plot sample PCA for component 1 and 2
plotIndiv(pca_all,
          comp = c(1:2),   # Specify components to plot
          ind.names = TRUE, # Show row names of samples
          title = 'bruxelles')

# see all the combinations for all the components
pca_all$variates$X %>%
  data.frame() %>%
  # filter only the top 5 PC
  # dplyr::select(PC1,PC2,PC3,PC4,PC5) %>%
  ggpairs(upper = "blank")+
  theme_bw() +
  theme(strip.background = element_blank())

# add some metadata from the original dataset to scatter plot component 1 and 2
id <- colnames(LUT_sample)[-1]
list_plot <- lapply(id,function(x){
  
  p1 <- pca_all$variates$X %>%
    data.frame() %>%
    rownames_to_column("sample") %>%
    left_join(LUT_sample,by = "sample") %>%
    ggplot(aes(x=PC1,y=PC2))+geom_point(aes(col=.data[[x]]))+
    theme_bw() +
    theme(strip.background = element_blank())
  
  return(p1)
})

wrap_plots(list_plot)

# pca_all$variates$X %>%
#   data.frame() %>%
#   rownames_to_column("sample") %>%
#   left_join(LUT_sample,by = "sample") %>%
#   ggplot(aes(x=PC1,y=PC2))+geom_point(aes(col=.data[["diagnosis"]]))+
#   theme_bw() +
#   theme(strip.background = element_blank())
# 
# pca_all$variates$X %>%
#   data.frame() %>%
#   rownames_to_column("sample") %>%
#   left_join(LUT_sample,by = "sample") %>%
#   ggplot(aes(x=PC1,y=PC2))+geom_point(aes(col=diagnosis))+
#   theme_bw() +
#   theme(strip.background = element_blank())
# 
# pca_all$variates$X %>%
#   data.frame() %>%
#   rownames_to_column("sample") %>%
#   left_join(LUT_sample,by = "sample") %>%
#   ggplot(aes(x=PC1,y=PC2))+geom_point(aes(col="diagnosis"))+
#   theme_bw() +
#   theme(strip.background = element_blank())
# 
# pca_all$variates$X %>%
#   data.frame() %>%
#   rownames_to_column("sample") %>%
#   left_join(LUT_sample,by = "sample") %>%
#   ggplot(aes(x=PC1,y=PC2))+geom_point(aes_string(col="diagnosis"))+
#   theme_bw() +
#   theme(strip.background = element_blank())


# plot circle plot for component 1 and 2
plotVar(pca_all, comp = c(1, 2),
        var.names = TRUE,
        cex = 3,         # To change the font size
        # cutoff = 0.5,  # For further cutoff
        title = 'PCA comp 1 - 2')

# biplot
test <- biplot(pca_all)

test$
test$layers[[4]]

list_data <- ggplot_build(test)

ggplot() +
  geom_point(data = list_data$data[[1]],aes(x,y),alpha=0.5) +
  geom_segment(data = list_data$data[[3]],aes(x = x,y = y,xend = xend,yend = yend),linetype = 1,col="red",alpha=0.5,arrow = arrow(length=unit(.2, 'cm'))) +
  ggrepel::geom_text_repel(data = list_data$data[[4]],aes(x,y,label = label),alpha=1,col="red",force = 0.5)+theme_bw()+
  xlab("PC1")+
  ylab("PC2")

# clustering --------------------------------------------------------------
# library("factoextra")
res <- get_clust_tendency(df_tot_noNA_scale, 40, graph = TRUE)

# Let’s see the Hopskin statistic
res$hopkins_stat

# Now we can visualize the dissimilarity matrix
print(res$plot)

# kmean
set.seed(1234)
# Ok now we can compute the gap statistic
gap_stat <- clusGap(df_tot_noNA_scale, FUN = kmeans, nstart = 25,K.max = 10, B = 100)

# Plot the result
# library(factoextra)
fviz_gap_stat(gap_stat)

# A 4-cluster solution is suggested by the gap statistic.
# The method NbClust() [in the NbClust] package can also be used. Calculate the clustering using the k-means algorithm. Clustering with K-means with k = 4:

# kmean -------------------------------------------------------------------
set.seed(1234)
km.res <- kmeans(df_tot_noNA_scale, 4, nstart = 25)
head(km.res$cluster, 20)

# Now visualize clusters using factoextra
fviz_cluster(km.res, df_tot_noNA_scale,ellipse = T)+theme_bw()

sil <- silhouette(km.res$cluster, dist(df_tot_noNA_scale))
rownames(sil) <- rownames(df_tot_noNA_scale)
head(sil[, 1:3])

fviz_silhouette(sil)

# plot as heatmap
table(km.res$cluster)
column_ha1 <- HeatmapAnnotation(module = km.res$cluster,
                                col = list(module = c("1" = "#F8766D",
                                                      "2" = "#7CAE00",
                                                      "3" = "#00BFC4",
                                                      "4" = "#C77CFF")))

ht1 <- Heatmap(t(df_tot_noNA_scale), cluster_rows = F,
               cluster_columns = F,
               show_heatmap_legend = TRUE,top_annotation = column_ha1,
               column_split = km.res$cluster)

draw(ht1)

# plot the individula variables per cluster
df_tot_noNA_scale %>%
  data.frame() %>%
  rownames_to_column("sample") %>%
  mutate(cluster = as.factor(km.res$cluster)) %>%
  pivot_longer(names_to = "variable",values_to = "value",-c(sample,cluster)) %>%
  ggplot(aes(x=cluster,y=value))+geom_boxplot(outlier.shape = NA)+geom_point(position = position_jitter(width = 0.1))+facet_wrap(~variable,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank())

# ht1 <- Heatmap(t(df_tot_noNA_scale), cluster_rows = F,
#                cluster_columns = T,
#                show_heatmap_legend = TRUE,top_annotation = column_ha1,column_km = 4)
# 
# draw(ht1)

# eclust kmean ------------------------------------------------------------
set.seed(1234)
km.res2 <- eclust(df_tot_noNA_scale, "kmean")
head(km.res2$cluster, 20)

# Now visualize clusters using factoextra
fviz_cluster(km.res2, df_tot_noNA_scale,ellipse = T)+theme_bw()

sil2 <- silhouette(km.res2$cluster, dist(df_tot_noNA_scale))
rownames(sil2) <- rownames(df_tot_noNA_scale)
head(sil2[, 1:3])

fviz_silhouette(sil2)

# plot as heatmap
table(km.res2$cluster)
show_col(hue_pal()(2))
column_ha12 <- HeatmapAnnotation(module = km.res2$cluster,
                                 col = list(module = c("1" = "#F8766D",
                                                       "2" = "#00BFC4")))

ht12 <- Heatmap(t(df_tot_noNA_scale), cluster_rows = F,
                cluster_columns = F,
                show_heatmap_legend = TRUE,top_annotation = column_ha12,
                column_split = km.res2$cluster)

draw(ht12)

# plot the individula variables per cluster
df_tot_noNA_scale %>%
  data.frame() %>%
  rownames_to_column("sample") %>%
  mutate(cluster = as.factor(km.res2$cluster)) %>%
  pivot_longer(names_to = "variable",values_to = "value",-c(sample,cluster)) %>%
  ggplot(aes(x=cluster,y=value))+geom_boxplot(outlier.shape = NA)+geom_point(position = position_jitter(width = 0.1))+facet_wrap(~variable,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank())

# Hierarchical clustering eclust ------------------------------------------
# Hierarchical clustering using eclust() Enhanced hierarchical clustering
res.hc <- eclust(df_tot_noNA_scale, "hclust")
res.hc

# fviz_dend(res.hc, rect = TRUE)
modules <- cutree(res.hc,k = 4)
column_dend <- as.dendrogram(res.hc)
# plot the dendrogram
fviz_dend(column_dend, rect = TRUE,k = 4,k_colors = c("#C77CFF","#7CAE00","#00BFC4","#F8766D"))

show_col(hue_pal()(4))
column_ha2 <- HeatmapAnnotation(module = modules,
                               col = list(module = c("1" = "#F8766D",
                                                     "2" = "#7CAE00",
                                                     "3" = "#00BFC4",
                                                     "4" = "#C77CFF")))

ht2 <- Heatmap(t(df_tot_noNA_scale), cluster_rows = F,
               cluster_columns = column_dend,
               show_heatmap_legend = TRUE,top_annotation = column_ha2)

draw(ht2)

ht22 <- Heatmap(t(df_tot_noNA_scale), cluster_rows = F,
                cluster_columns = F,
                show_heatmap_legend = TRUE,top_annotation = column_ha2,column_split = modules)

draw(ht22)

# Now visualize clusters using factoextra on the pca
fviz_cluster(object = list(data = df_tot_noNA_scale,cluster = modules),ellipse = T)+theme_bw()

# explore variables per cluster
df_tot_noNA_scale %>%
  data.frame() %>%
  rownames_to_column("sample") %>%
  mutate(cluster = as.factor(modules)) %>%
  pivot_longer(names_to = "variable",values_to = "value",-c(sample,cluster)) %>%
  ggplot(aes(x=cluster,y=value))+geom_boxplot(outlier.shape = NA)+geom_point(position = position_jitter(width = 0.1))+facet_wrap(~variable,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank())
