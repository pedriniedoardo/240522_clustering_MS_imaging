# # PCA ---------------------------------------------------------------------
# # 1. Load the nutrimouse data from the mixOmics R package and investigate its structure.
# library(mixOmics)
# # A data object provided by an R package can be loaded with data. Its structure can be obainted with str, length, dim, etc.
# data("nutrimouse")
# ## display the structure of the nutrimouse object
# str(nutrimouse)
# 
# ## check dimensions
# lapply(nutrimouse, dim) # apply function dim to each element in list nutrimouse
# lapply(nutrimouse, length) # apply function length to each element in list nutrimouse
# 
# # 2. Take the gene expression dataset in samples x variables matrix format. Investigate variables’ distribution.
# ## get gene expression data structure
# class(nutrimouse$gene)
# dim(nutrimouse$gene)
# rownames(nutrimouse$gene)[1:10]
# colnames(nutrimouse$gene)[1:10]
# ## check if there are missing values
# any(is.na(nutrimouse$gene))
# 
# ## investigate each variable
# summary(nutrimouse$gene[, 1:4])
# 
# colors <- rainbow(20, alpha=1)
# plot(density(scale(nutrimouse$gene[, 1], center=T, scale=F)), 
#      col=colors[1], xlim=c(-0.5,0.5), ylim=c(0,8),
#      main='Density plot of the first 20 genes')
# invisible(sapply(2:20, function(i) {
#   lines(density(scale(nutrimouse$gene[, i], center=T, scale=F)), col=colors[i])
# }))
# 
# 
# # 3. Perform PCA and investigate variances, sample distribution and variable relationship with plots.
# # A number of methods in different R packages can perform PCA, e.g. stats::prcomp, stats::princomp, mixOmics::pca, multiblock::pca, psych::principal, FactoMineR::PCA, etc.
# 
# pca.res <- prcomp(nutrimouse$gene, center=TRUE, scale.=F)
# names(pca.res)
# summary(pca.res)
# 
# # Variances = eigenvalues of the covariance matrix = (standard deviation)^2.
# 
# variances <- pca.res$sdev^2
# variances
# 
# # Scree plot: plot of variances.
# 
# screeplot(pca.res, npcs=length(variances), type='lines')
# screeplot(pca.res, npcs=length(variances), type='barplot')
# barplot(variances, xlab='PC', ylab='Variance', names.arg=1:length(variances))
# 
# 
# # Scree plot on variance percentage.
# 
# varPercent <- variances/sum(variances) * 100
# barplot(varPercent, xlab='PC', ylab='Percent Variance', names.arg=1:length(varPercent))
# 
# 
# # Scores: sample coordinates in the new reference (rotated axes or principal components).
# 
# scores <- pca.res$x
# str(scores)
# 
# # Score plot: plot of sample distribution.
# 
# PCx <- "PC1"
# PCy <- "PC2"
# plot(scores[, PCx], scores[, PCy], xlab=PCx, ylab=PCy, pch=16)
# text(scores[, PCx], scores[, PCy]-0.05, rownames(scores), col='blue', cex=0.7)
# 
# # Loadings: contributions of variables to principal components (eigenvectors of covariance matrix).
# loadings <- pca.res$rotation
# str(loadings)
# 
# # Loading plot: plot of variables’ contribution, revealing their relationship.
# plot(loadings[, PCx], loadings[, PCy], type='n', main="Loadings")
# arrows(0, 0, loadings[, PCx], loadings[, PCy], xlab=PCx, ylab=PCy,
#        length=0.1, angle=20, col=rgb(0,0,1,alpha=apply(loadings[, c(PCx, PCy)], 1, norm, "2")))
# text(loadings[, PCx], loadings[, PCy], rownames(loadings), col='grey', cex=0.7)
# 
# 
# Both score and loading plots can be plotted altogether with the biplot function.
# 
# ## biplot
# biplot(pca.res, expand=1, cex=c(0.5, 0.7), col=c("gray50", "red"))
# 
# library(factoextra)
# fviz_pca_biplot(pca.res, repel = TRUE,
#                 col.var = "blue", # Variables color
#                 habillage = nutrimouse$genotype,
#                 addEllipses = T,
#                 legend="none"
# )
# 
# 
# # 4. Visually investigate the sample distribution with coloring by metadata or expression of certain genes.
# # The samples can be colored with some metadata, e.g genotype or diet,
# plot(scores[, PCx], scores[, PCy], main="Scores",
#      col=c(1:nlevels(nutrimouse$diet))[nutrimouse$diet],
#      pch=c(17,19)[nutrimouse$genotype],
#      xlab=paste0(PCx,
#                  " (",
#                  round((summary(pca.res)$importance)[2, PCx], 2), ")"),
#      ylab=paste0(PCy,
#                  " (",
#                  round((summary(pca.res)$importance)[2, PCy], 2), ")")
# )
# legend("topright", title="genotype",
#        legend=levels(nutrimouse$genotype),
#        pch=c(17,19), cex=0.7)
# legend("bottomright", title="diet",
#        legend=levels(nutrimouse$diet),
#        col=c(1:5), cex=0.7, pch=16)
# 
# # or by some gene expression.
# nbreaks <- 5
# plot(scores[, PCx], scores[, PCy], xlab=PCx, ylab=PCy, 
#      pch=c(17,19)[nutrimouse$genotype], 
#      col=colorRampPalette(c('red','blue'))(nbreaks)[as.numeric(cut(nutrimouse$gene$ALDH3,breaks = nbreaks))])
# 
# 
