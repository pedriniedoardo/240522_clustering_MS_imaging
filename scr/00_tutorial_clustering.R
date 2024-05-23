# # the tutorial comes form the following reference
# # https://finnstats.com/index.php/2022/05/31/clustering-example-step-by-step-methods-in-r/?utm_source=ReviveOldPost&utm_medium=social&utm_campaign=ReviveOldPost
# library(cluster)
# library(factoextra)
# library(tidyverse)
# 
# # A thorough cluster analysis can be carried out in three steps, as listed below:
# 
# # Preparation of data -----------------------------------------------------
# # Identifying the likelihood of grouping (i.e., the clusterability of the data)
# # Choosing the best number of clusters
# # Cluster analysis (e.g., k-means, pam) or hierarchical clustering can be computed using partitioning or hierarchical clustering.
# # Analyses of clustering need to be validated. plan with silhouettes
# # We give R scripts to conduct all of these processes here.
# 
# # Preparation of data -----------------------------------------------------
# # We’ll utilize the USArrests demo data set. We begin by using the scale() method to standardize the data:
# # Let’s load the data set
# data(USArrests)
# 
# # Standardize the data set
# df <- scale(USArrests)
# head(df)
# 
# # plot the effect of the standardization of the data
# df_raw <- USArrests %>% 
#   # mutate(test_murder=(Murder - mean(Murder))/sd(Murder)) %>% 
#   rownames_to_column("state") %>% 
#   pivot_longer(names_to = "var",values_to = "val",-state) %>% 
#   group_by(var) %>% 
#   mutate(test_norm =(val - mean(val))/sd(val))
# 
# df_raw
# 
# df_scale <- df %>% 
#   data.frame() %>% 
#   rownames_to_column("state") %>% 
#   pivot_longer(names_to = "var",values_to = "val_scale",-state)
# 
# df_scale
# 
# # join the two objects
# test <- left_join(df_raw,df_scale,by=c("state","var"))
# test
# 
# # notice that the scale function is z scoring the variables
# test %>% 
#   ggplot(aes(x=val_scale,y=test_norm))+geom_point()
# 
# # notice that the scale function is simply scaling the variable, is not really balancing the different entries
# test %>% 
#   ggplot(aes(x=var,y=val))+geom_violin()
# 
# test %>% 
#   ggplot(aes(x=var,y=val_scale))+geom_violin()
# 
# # Assessing the clusterability --------------------------------------------
# # The [factoextra package] function get_clust_tendency() can be used. It calculates the Hopkins statistic and offers a visual representation.
# 
# # library("factoextra")
# res <- get_clust_tendency(df, 40, graph = TRUE)
# 
# # Let’s see the Hopskin statistic
# res$hopkins_stat
# 
# # Now we can visualize the dissimilarity matrix
# print(res$plot)
# 
# # The Hopkins statistic has a value of much less than 0.5, indicating that the data is highly clusterable. The ordered dissimilarity image also contains patterns, as can be observed (i.e., clusters).
# 
# # Calculate how many clusters there are in the data -----------------------
# # We’ll use the function clusGap() [cluster package] to compute gap statistics for estimating the optimal number of clusters because k-means clustering requires specifying the number of clusters to generate.
# # Customer Segmentation K Means Cluster »
# # The gap statistic plot is visualised using the function fviz_gap_stat()  [factoextra].
# # library("cluster")
# set.seed(1234)
# # Ok now we can compute the gap statistic
# gap_stat <- clusGap(df, FUN = kmeans, nstart = 25,K.max = 10, B = 100)
# 
# # Plot the result
# # library(factoextra)
# fviz_gap_stat(gap_stat)
# 
# # A 3-cluster solution is suggested by the gap statistic.
# # The method NbClust() [in the NbClust] package can also be used. Calculate the clustering using the k-means algorithm. Clustering with K-means with k = 3:
# # Compute k-means
# set.seed(1234)
# km.res <- kmeans(df, 3, nstart = 25)
# head(km.res$cluster, 20)
# 
# # Now visualize clusters using factoextra
# fviz_cluster(km.res, USArrests)
# 
# # Statistics on cluster validation: Examine the cluster silhouette graph. Remember that the silhouette compares an object’s similarity to other items in its own cluster to those in a neighboring cluster (Si). Si values range from 1 to -1. If the value of Si is near one, the object is well clustered. A Si value near -1 suggests that the object is poorly grouped and that assigning it to a different cluster would most likely enhance the overall results. 
# 
# # Cluster Analysis in R » Unsupervised Approach » -------------------------
# sil <- silhouette(km.res$cluster, dist(df))
# rownames(sil) <- rownames(USArrests)
# head(sil[, 1:3])
# 
# fviz_silhouette(sil)
# 
# # There are also samples with negative silhouette values, as may be observed. The following are some typical inquiries:
# # What kinds of samples are these? What cluster do they belong to?
# # This can be deduced from the output of the silhouette() function as follows:
# neg_sil_index <- which(sil[, "sil_width"] < 0)
# sil[neg_sil_index, , drop = FALSE]
# 
# # eclust() is a function that improves clustering analysis ----------------
# # When compared to traditional clustering analysis packages, the function eclust()[factoextra package] offers various advantages:
# # It streamlines the clustering analysis workflow. In a single line function call, it may perform hierarchical clustering and partitioning clustering. The function eclust() computes the gap statistic for predicting the correct number of clusters automatically. It delivers silhouette information automatically. It creates stunning graphs with ggplot2 and eclust K-means clustering () Compute k-means
# 
# res.km <- eclust(df, "kmeans", nstart = 25)
# res.km
# 
# # Gap statistic plot
# fviz_gap_stat(res.km$gap_stat)
# 
# # Silhouette plot
# fviz_silhouette(res.km)
# 
# # Hierarchical clustering using eclust() Enhanced hierarchical clustering
# res.hc <- eclust(df, "hclust")
# res.hc
# 
# fviz_dend(res.hc, rect = TRUE)
# 
# modules <- cutree(res.hc,k = 3)
# column_dend <- as.dendrogram(res.hc)
# show_col(hue_pal()(3))
# column_ha <- HeatmapAnnotation(module = modules,
#                                col = list(module = c("1" = "#F8766D",
#                                                      "3" = "#00BA38",
#                                                      "2" = "#619CFF")))
# 
# ht1 <- Heatmap(t(df), cluster_rows = F,
#                cluster_columns = column_dend,
#                show_heatmap_legend = TRUE,top_annotation = column_ha)
# 
# draw(ht1)
# 
# 
# ht2 <- Heatmap(t(df), cluster_rows = F,
#                cluster_columns = F,
#                show_heatmap_legend = TRUE,top_annotation = column_ha,column_split = modules)
# 
# draw(ht2)
