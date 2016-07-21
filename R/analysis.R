##################################
# R source code file for analyzing TCGAOV data
# Created by Sahir, July 21, 2016
# Updated:
# NOTE: this is on the TCGAOV repo, master branch
# 
# 
##################################


rm(list = ls())


# Get data ----------------------------------------------------------------

source("R/data.R")


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
                     nPC = 2)


# to combine all the principal components
pcTrain_TOM <- pp$clustersAddon$PC
dim(pcTrain_TOM)
head(pcTrain_TOM)
avgTrain_TOM <- pp$clustersAddon$averageExpr

# to combine all the principal components
pcTrain_corr <- pp2$clustersAddon$PC
dim(pcTrain_corr)
head(pcTrain_corr)
avgTrain_corr <- pp2$clustersAddon$averageExpr

varexp_PC1_TOM <- pp$clustersAddon$varExplained[seq(1, length(pp$clustersAddon$varExplained), by = 2)]
varexp_PC2_TOM <- pp$clustersAddon$varExplained[seq(2, length(pp$clustersAddon$varExplained), by = 2)]

varexp_PC1_corr <- pp2$clustersAddon$varExplained[seq(1, length(pp2$clustersAddon$varExplained), by = 2)]
varexp_PC2_corr <- pp2$clustersAddon$varExplained[seq(2, length(pp2$clustersAddon$varExplained), by = 2)]


dTOM <- data.frame(index = seq_len(length(varexp_PC1_TOM)), varexp_PC1_TOM, varexp_PC2_TOM) %>%
  gather(type, value, -index) %>%
  separate(type, c("measure", "PC", "type"))
dcorr <- data.frame(index = seq_len(length(varexp_PC1_corr)), varexp_PC1_corr, varexp_PC2_corr) %>%
  gather(type, value, -index) %>%
  separate(type, c("measure", "PC", "type"))

var_expl_data <- rbind(dTOM, dcorr)

p <- ggplot(var_expl_data, aes(x = index, y = value, color = PC))
p + geom_point(size=2) + facet_wrap(~type) + ylab("variance explained") + theme_bw()

# plot(seq_len(length(var_exp_pc1)), var_exp_pc1, pch = 19, col = "red",
#      ylim = c(0, max(var_exp_pc1, var_exp_pc2)),
#      xlab = "cluster index", ylab = "variance explained")
# points(seq_len(length(var_exp_pc1)), var_exp_pc2, pch = 19, col = "blue")
# legend("topright", c("PC1", "PC2"), col = c("red","blue"),
#        pch = c(19,19), bg = "gray90")

colnames(pcTrain_TOM) <- gsub("\\.","_",colnames(pcTrain_TOM))
# colnames(pcTrain_TOM) <- paste0("PC",colnames(pcTrain_TOM))
datt <- cbind(pcTrain_TOM, Y = Y[trainIndex], age = DT.pheno.placenta$Age_gestationnel, sex = DT.pheno.placenta$Sexe)
colnames(datt)
str(datt)



# sent to Celia July 13
bouchard_2PC_10K_probes_TOM_DIFFTOM <- cbind(pcTrain_TOM,
                                            bmizscore = Y[trainIndex],
                                            gdstatus = E[trainIndex])

write.csv(bouchard_2PC_10K_probes_TOM_DIFFTOM, file = "bouchard_2PC_10K_TOM_DIFFTOM.csv",
          quote = FALSE,row.names = FALSE)
head(bouchard_2PC_10K_probes_TOM_DIFFTOM)


max_heat <- max(c(max(pcTrain_TOM[which(pp$etrain==0),]),max(pcTrain_TOM[which(pp$etrain==1),])))
min_heat <- min(c(min(pcTrain_TOM[which(pp$etrain==0),]),min(pcTrain_TOM[which(pp$etrain==1),])))

pheatmap::pheatmap(t(pcTrain_TOM[which(pp$etrain==0),]),
                   clustering_method = "average",
                   color = viridis(100),
                   breaks = seq(min_heat, max_heat, length.out = 101))
pheatmap::pheatmap(t(pcTrain_TOM[which(pp$etrain==1),]),
                   clustering_method = "average",
                   color = viridis(100),
                   breaks = seq(min_heat, max_heat, length.out = 101))


library(ComplexHeatmap)
require(circlize)

cm <- colorRamp2(seq(min_heat, max_heat, length.out = 100), viridis(100))
ht1 = Heatmap(t(pcTrain_TOM[which(pp$etrain==0),]),
              name = "E=0",
              # col = viridis(10),
              col = cm,
              # column_title = "E = 0 : Age [4.8, 11.3]",
              # column_title = "Income_Level: 1-7",
              column_title = "NGD",
              show_row_names = FALSE)
ht2 = Heatmap(t(pcTrain_TOM[which(pp$etrain==1),]),
              name = "E=1",
              # col = viridis(10),
              col = cm,
              # column_title = "E = 1 : Age [11.3, 18]",
              # column_title = "Income_Level: 8-10",
              column_title = "GD",
              show_row_names = FALSE)
ht1 + ht2




