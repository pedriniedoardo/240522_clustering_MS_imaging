# # load the data -----------------------------------------------------------
# library(mixOmics)
# data(multidrug)
# X <- multidrug$ABC.trans
# dim(X) # Check dimensions of data
# 
# # PCA ---------------------------------------------------------------------
# # 3.2.1 Choose the number of components
# tune.pca.multi <- tune.pca(X, ncomp = 10, scale = TRUE)
# plot(tune.pca.multi)
# 
# # 3.2.2 PCA with fewer components
# final.pca.multi <- pca(X, ncomp = 3, center = TRUE, scale = TRUE)
# final.pca.multi$var.tot
# final.pca.multi$prop_expl_var$X
# 
# # 3.2.3 Identify the informative variables
# head(selectVar(final.pca.multi, comp = 1)$value)
# 
# # 3.2.4 Sample plots
# plotIndiv(final.pca.multi,
#           comp = c(1, 2),   # Specify components to plot
#           ind.names = TRUE, # Show row names of samples
#           group = multidrug$cell.line$Class,
#           title = 'ABC transporters, PCA comp 1 - 2',
#           legend = TRUE, legend.title = 'Cell line')
# 
# plotVar(final.pca.multi, comp = c(1, 2),
#         var.names = TRUE,
#         cex = 3,         # To change the font size
#         # cutoff = 0.5,  # For further cutoff
#         title = 'Multidrug transporter, PCA comp 1 - 2')
# 
# biplot(final.pca.multi, group = multidrug$cell.line$Class, 
#        legend.title = 'Cell line')
