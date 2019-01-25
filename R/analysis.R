##################################
# R source code file for analyzing TCGAOV data
# Created by Sahir, July 21, 2016
# Updated: September 28th, 2016
# now being run on mammouth cluster
# NOTE: this is on the TCGAOV repo, master branch
# 
# 
##################################


# rm(list = ls())


# Get data ----------------------------------------------------------------

# source("/home/bhatnaga/coexpression/rda/tcgaov/data.R")
source(paste(Sys.getenv("PBS_O_WORKDIR"),"data.R", sep="/"))


# Create training and test set indices ------------------------------------

Y <- DT_final$OS %>% as.matrix() %>% drop

any(is.na(Y))

trainIndex <- caret::createDataPartition(Y, p = 1, list = FALSE, times = 1) %>% drop

# not really used... need to do LOOCV 
# testIndex <- which(seq_len(length(Y)) %ni% trainIndex)
testIndex <- trainIndex
# all(c(trainIndex, testIndex) %in% seq_len(length(Y)))

# we also need the exposure variable as a vector
# needs to be 0 and 1
E <- DT_final$E



# Cluster -----------------------------------------------------------------

pp <- cluster_kmeans(data = as.matrix(DT_final[,-c("rn","subtype","E","status"),with=F]), 
                     exposure = E,
                     response = Y,
                     min_cluster_size = 50,
                     train_index = trainIndex, 
                     test_index = testIndex,
                     cluster_distance = "tom",
                     cluster_method = "hclust",
                     cut_method = "dynamic",
                     agglomeration_method = "average",
                     distance_method = "euclidean",
                     eclust_add = TRUE,
                     eclust_distance = "difftom",
                     nPC = 1)

# save(pp, file = "/home/bhatnaga/coexpression/rda/tcgaov/PC_TCGAOV_881_signature_genes_TOM_DIFFTOM.RData")
save(pp, file = paste(Sys.getenv("PBS_O_WORKDIR"),"PC_TCGAOV_881_signature_genes_TOM_DIFFTOM.RData", sep="/"))

# just calling this res5k, even though its actually only 881 genes
res5k <- pp


# Different input data (all based on TOM, not using correlations, no test set) ----

# 5k most variable 
pc_train_5k_eclust <- res5k$clustersAddon$PC
pc_train_5k_clust <- res5k$clustersAll$PC
avg_train_5k_eclust <- res5k$clustersAddon$averageExpr
avg_train_5k_clust <- res5k$clustersAll$averageExpr
original_train_5k <- as.matrix(DT_final[,-c("rn","subtype","E","status"),with=F])


# eclust and clust interaction data 
pc_train_5k_eclust_interaction <- prepare_data(data = cbind(pc_train_5k_eclust, 
                                                            pheno = Y[trainIndex],
                                                            income = E[trainIndex]),
                                               response = "pheno", exposure = "income")

pc_train_5k_clust_interaction <- prepare_data(data = cbind(pc_train_5k_clust, 
                                                           pheno = Y[trainIndex],
                                                           income = E[trainIndex]),
                                              response = "pheno", exposure = "income")

avg_train_5k_eclust_interaction <- prepare_data(data = cbind(avg_train_5k_eclust, 
                                                             pheno = Y[trainIndex],
                                                             income = E[trainIndex]),
                                                response = "pheno", exposure = "income")

avg_train_5k_clust_interaction <- prepare_data(data = cbind(avg_train_5k_clust, 
                                                            pheno = Y[trainIndex],
                                                            income = E[trainIndex]),
                                               response = "pheno", exposure = "income")

# original interaction data
original_train_5k_interaction <- prepare_data(data = cbind(as.data.frame(original_train_5k), 
                                                           pheno = Y[trainIndex],
                                                           income = E[trainIndex]),
                                              response = "pheno", exposure = "income")


# Fit models --------------------------------------------------------------

# fitControl <-  trainControl(method = "repeatedcv",
#                             number = 5,
#                             repeats = 3,
#                             verboseIter = TRUE)

fitControl <-  trainControl(method = "boot632",
                            number = 100,
                            # repeats = 3,
                            verboseIter = TRUE)

# Define the candidate models to test
marsGrid <- expand.grid(.degree = 1:2, .nprune = 1000)

#############################################################################
#               MARS Original                                            ----
#############################################################################
set.seed(1056)
mars_original_5k <- train(original_train_5k_interaction$X[,original_train_5k_interaction$main_effect_names], 
                          original_train_5k_interaction$Y,
                          method = "earth",
                          trace = 1, nk = 1000, keepxy = TRUE, pmethod = "backward",
                          tuneGrid = marsGrid,
                          trControl = fitControl)

#############################################################################
#               MARS Clust                                              ----
#############################################################################
set.seed(1056)
# marsGridClust <- expand.grid(.degree = 1, .nprune = 1000)
mars_pc_5k_clust <- train(pc_train_5k_clust_interaction$X[, pc_train_5k_clust_interaction$main_effect_names], 
                          pc_train_5k_clust_interaction$Y,
                          method = "earth",
                          trace = 1, nk = 1000, keepxy = TRUE, pmethod = "backward",
                          tuneGrid = marsGrid,
                          trControl = fitControl)

set.seed(1056)
mars_avg_5k_clust <- train(avg_train_5k_clust_interaction$X[, avg_train_5k_clust_interaction$main_effect_names], 
                           avg_train_5k_clust_interaction$Y,
                           method = "earth",
                           trace = 1, nk = 1000, keepxy = TRUE, pmethod = "backward",
                           tuneGrid = marsGrid,
                           trControl = fitControl)

#############################################################################
#                     MARS EClust                                        ----
#############################################################################
set.seed(1056)
mars_pc_5k_eclust <- train(pc_train_5k_eclust_interaction$X[, pc_train_5k_eclust_interaction$main_effect_names], 
                           pc_train_5k_eclust_interaction$Y,
                           method = "earth",
                           trace = 1, nk = 1000, keepxy = TRUE, pmethod = "backward",
                           # Explicitly declare the candidate models to test
                           tuneGrid = marsGrid,
                           trControl = fitControl)

set.seed(1056)
mars_avg_5k_eclust <- train(avg_train_5k_eclust_interaction$X[, avg_train_5k_eclust_interaction$main_effect_names], 
                            avg_train_5k_eclust_interaction$Y,
                            method = "earth",
                            trace = 1, nk = 1000, keepxy = TRUE, pmethod = "backward",
                            # Explicitly declare the candidate models to test
                            tuneGrid = marsGrid,
                            trControl = fitControl)


#############################################################################
#               Bagged MARS Original                                    ----
#############################################################################
set.seed(1056)
bag_mars_original_5k <- train(original_train_5k_interaction$X[,original_train_5k_interaction$main_effect_names], 
                              original_train_5k_interaction$Y,
                              method = "bagEarth", B = 50,
                              trace = 1, nk = 1000, keepxy = TRUE, pmethod = "backward",
                              tuneGrid = marsGrid,
                              trControl = fitControl)


#############################################################################
#               Bagged MARS Clust                                        ----
#############################################################################
set.seed(1056)
bag_mars_pc_5k_clust <- train(pc_train_5k_clust_interaction$X[, pc_train_5k_clust_interaction$main_effect_names], 
                              pc_train_5k_clust_interaction$Y,
                              method = "bagEarth", B = 50,
                              trace = 1, nk = 1000, keepxy = TRUE, pmethod = "backward",
                              # Explicitly declare the candidate models to test
                              tuneGrid = marsGrid,
                              trControl = fitControl)

set.seed(1056)
bag_mars_avg_5k_clust <- train(avg_train_5k_clust_interaction$X[, avg_train_5k_clust_interaction$main_effect_names], 
                               avg_train_5k_clust_interaction$Y,
                               method = "bagEarth", B = 50,
                               trace = 1, nk = 1000, keepxy = TRUE, pmethod = "backward",
                               # Explicitly declare the candidate models to test
                               tuneGrid = marsGrid,
                               trControl = fitControl)

#############################################################################
#               Bagged MARS EClust                                        ----
#############################################################################
set.seed(1056)
bag_mars_pc_5k_eclust <- train(pc_train_5k_eclust_interaction$X[, pc_train_5k_eclust_interaction$main_effect_names], 
                               pc_train_5k_eclust_interaction$Y,
                               method = "bagEarth", B = 50,
                               trace = 1, nk = 1000, keepxy = TRUE, pmethod = "backward",
                               # Explicitly declare the candidate models to test
                               tuneGrid = marsGrid,
                               trControl = fitControl)

set.seed(1056)
bag_mars_avg_5k_eclust <- train(avg_train_5k_eclust_interaction$X[, avg_train_5k_eclust_interaction$main_effect_names], 
                                avg_train_5k_eclust_interaction$Y,
                                method = "bagEarth", B = 50,
                                trace = 1, nk = 1000, keepxy = TRUE, pmethod = "backward",
                                # Explicitly declare the candidate models to test
                                tuneGrid = marsGrid,
                                trControl = fitControl)


# #############################################################################
# #               Lasso Original                                            ----
# #############################################################################
# 
# # Lasso + ECLUST ----
# # apparently this does tune over both alpha and lambda, but just need to provide lambda
# # http://stats.stackexchange.com/questions/69638/does-caret-train-function-for-glmnet-cross-validate-for-both-alpha-and-lambda?rq=1
# # lassoGrid <- expand.grid(.alpha = seq(0,1,length.out = 5), .lambda = seq(min(glmnet(x = kl_5k$X, y = kl_5k$Y)$lambda),
# #                                                               max(glmnet(x = kl_5k$X, y = kl_5k$Y)$lambda),length.out = 5))
# # glmnet(x = kl_5k$X, y = kl_5k$Y)$lambda
# 
# set.seed(1056)
# lasso_original_5k <- train(x = original_train_5k_interaction$X, 
#                            y = original_train_5k_interaction$Y,
#                            method = "glmnet",
#                            tuneLength = 20,
#                            trControl = fitControl)
# 
# #############################################################################
# #               Lasso Clust                                            ----
# #############################################################################
# 
# set.seed(1056)
# lasso_pc_5k_clust <- train(x = pc_train_5k_clust_interaction$X, 
#                            y = pc_train_5k_clust_interaction$Y,
#                            method = "glmnet",
#                            tuneLength = 20,
#                            trControl = fitControl)
# 
# set.seed(1056)
# lasso_avg_5k_clust <- train(x = avg_train_5k_clust_interaction$X, 
#                             y = avg_train_5k_clust_interaction$Y,
#                             method = "glmnet",
#                             tuneLength = 20,
#                             trControl = fitControl)
# 
# #############################################################################
# #               Lasso EClust                                            ----
# #############################################################################
# 
# set.seed(1056)
# lasso_pc_5k_eclust <- train(x = pc_train_5k_eclust_interaction$X, 
#                             y = pc_train_5k_eclust_interaction$Y,
#                             method = "glmnet",
#                             tuneLength = 20,
#                             trControl = fitControl)
# 
# set.seed(1056)
# lasso_avg_5k_eclust <- train(x = avg_train_5k_eclust_interaction$X, 
#                              y = avg_train_5k_eclust_interaction$Y,
#                              method = "glmnet",
#                              tuneLength = 20,
#                              trControl = fitControl)




resamps <- resamples(list(mars_original_5k=mars_original_5k,
                          mars_pc_5k_clust=mars_pc_5k_clust,
                          mars_avg_5k_clust=mars_avg_5k_clust,
                          mars_pc_5k_eclust=mars_pc_5k_eclust,
                          mars_avg_5k_eclust=mars_avg_5k_eclust,
                          bag_mars_pc_5k_clust=bag_mars_pc_5k_clust,
                          bag_mars_avg_5k_clust=bag_mars_avg_5k_clust,
                          bag_mars_pc_5k_eclust=bag_mars_pc_5k_eclust,
                          bag_mars_avg_5k_eclust=bag_mars_avg_5k_eclust))
# lasso_original_5k=lasso_original_5k,
# lasso_pc_5k_clust=lasso_pc_5k_clust,
# lasso_avg_5k_clust=lasso_avg_5k_clust,
# lasso_pc_5k_eclust=lasso_pc_5k_eclust,
# lasso_avg_5k_eclust=lasso_avg_5k_eclust))

# save(resamps, file = "/home/bhatnaga/coexpression/rda/tcgaov/TCGAOV_resamps_881_signature_genes.RData")
save(resamps, file = paste(Sys.getenv("PBS_O_WORKDIR"),"TCGAOV_resamps_881_signature_genes.RData", sep="/"))

# summary(resamps)
# trellis.par.set(caretTheme())
# dotplot(resamps)
# xyplot(resamps, what = "BlandAltman")
# difValues <- diff(resamps)
# difValues
# summary(difValues)
# bwplot(difValues,scales = list(x = list(relation = "free")))
# dotplot(difValues,scales = list(x = list(relation = "free")))
# dev.off()
# trellis.par.set("theme1")

# png("/home/bhatnaga/coexpression/rda/tcgaov/boxplot_resamps_tcgaov_881_signature_genes.png")
png(paste(Sys.getenv("PBS_O_WORKDIR"),"boxplot_resamps_tcgaov_881_signature_genes.png", sep="/"))
bwplot(resamps, horizontal=T, scales = list(x = list(relation = "free")))
dev.off()

# png("/home/bhatnaga/coexpression/rda/tcgaov/dotplot_resamps_tcgaov_881_signature_genes.png")
png(paste(Sys.getenv("PBS_O_WORKDIR"),"dotplot_resamps_tcgaov_881_signature_genes.png", sep="/"))
dotplot(resamps, scales = list(x = list(relation = "free")))
dev.off()


# 
# 
# caret:::bwplot.resamples(resamps, scales =list(x = list(relation = "free")))
# splom(resamps)
# densityplot(resamps)
# xyplot(resamps)
# caret:::densityplot.resamples(resamps, auto.key = list(columns = 3),pch = "|")
# # to combine all the principal components
# pcTrain_TOM <- pp$clustersAddon$PC
# dim(pcTrain_TOM)
# head(pcTrain_TOM)
# avgTrain_TOM <- pp$clustersAddon$averageExpr
# 
# # to combine all the principal components
# pcTrain_corr <- pp2$clustersAddon$PC
# dim(pcTrain_corr)
# head(pcTrain_corr)
# avgTrain_corr <- pp2$clustersAddon$averageExpr
# 
# varexp_PC1_TOM <- pp$clustersAddon$varExplained[seq(1, length(pp$clustersAddon$varExplained), by = 2)]
# varexp_PC2_TOM <- pp$clustersAddon$varExplained[seq(2, length(pp$clustersAddon$varExplained), by = 2)]
# 
# varexp_PC1_corr <- pp2$clustersAddon$varExplained[seq(1, length(pp2$clustersAddon$varExplained), by = 2)]
# varexp_PC2_corr <- pp2$clustersAddon$varExplained[seq(2, length(pp2$clustersAddon$varExplained), by = 2)]
# 
# 
# dTOM <- data.frame(index = seq_len(length(varexp_PC1_TOM)), varexp_PC1_TOM, varexp_PC2_TOM) %>%
#   gather(type, value, -index) %>%
#   separate(type, c("measure", "PC", "type"))
# dcorr <- data.frame(index = seq_len(length(varexp_PC1_corr)), varexp_PC1_corr, varexp_PC2_corr) %>%
#   gather(type, value, -index) %>%
#   separate(type, c("measure", "PC", "type"))
# 
# var_expl_data <- rbind(dTOM, dcorr)
# 
# p <- ggplot(var_expl_data, aes(x = index, y = value, color = PC))
# p + geom_point(size=2) + facet_wrap(~type) + ylab("variance explained") + theme_bw()
# 
# # plot(seq_len(length(var_exp_pc1)), var_exp_pc1, pch = 19, col = "red",
# #      ylim = c(0, max(var_exp_pc1, var_exp_pc2)),
# #      xlab = "cluster index", ylab = "variance explained")
# # points(seq_len(length(var_exp_pc1)), var_exp_pc2, pch = 19, col = "blue")
# # legend("topright", c("PC1", "PC2"), col = c("red","blue"),
# #        pch = c(19,19), bg = "gray90")
# 
# colnames(pcTrain_TOM) <- gsub("\\.","_",colnames(pcTrain_TOM))
# # colnames(pcTrain_TOM) <- paste0("PC",colnames(pcTrain_TOM))
# datt <- cbind(pcTrain_TOM, Y = Y[trainIndex], age = DT.pheno.placenta$Age_gestationnel, sex = DT.pheno.placenta$Sexe)
# colnames(datt)
# str(datt)
# 
# 
# 
# # sent to Celia July 13
# bouchard_2PC_10K_probes_TOM_DIFFTOM <- cbind(pcTrain_TOM,
#                                             bmizscore = Y[trainIndex],
#                                             gdstatus = E[trainIndex])
# 
# write.csv(bouchard_2PC_10K_probes_TOM_DIFFTOM, file = "bouchard_2PC_10K_TOM_DIFFTOM.csv",
#           quote = FALSE,row.names = FALSE)
# head(bouchard_2PC_10K_probes_TOM_DIFFTOM)
# 
# 
# max_heat <- max(c(max(pcTrain_TOM[which(pp$etrain==0),]),max(pcTrain_TOM[which(pp$etrain==1),])))
# min_heat <- min(c(min(pcTrain_TOM[which(pp$etrain==0),]),min(pcTrain_TOM[which(pp$etrain==1),])))
# 
# pheatmap::pheatmap(t(pcTrain_TOM[which(pp$etrain==0),]),
#                    clustering_method = "average",
#                    color = viridis(100),
#                    breaks = seq(min_heat, max_heat, length.out = 101))
# pheatmap::pheatmap(t(pcTrain_TOM[which(pp$etrain==1),]),
#                    clustering_method = "average",
#                    color = viridis(100),
#                    breaks = seq(min_heat, max_heat, length.out = 101))
# 
# 
# library(ComplexHeatmap)
# require(circlize)
# 
# cm <- colorRamp2(seq(min_heat, max_heat, length.out = 100), viridis(100))
# ht1 = Heatmap(t(pcTrain_TOM[which(pp$etrain==0),]),
#               name = "E=0",
#               # col = viridis(10),
#               col = cm,
#               # column_title = "E = 0 : Age [4.8, 11.3]",
#               # column_title = "Income_Level: 1-7",
#               column_title = "NGD",
#               show_row_names = FALSE)
# ht2 = Heatmap(t(pcTrain_TOM[which(pp$etrain==1),]),
#               name = "E=1",
#               # col = viridis(10),
#               col = cm,
#               # column_title = "E = 1 : Age [11.3, 18]",
#               # column_title = "Income_Level: 8-10",
#               column_title = "GD",
#               show_row_names = FALSE)
# ht1 + ht2




