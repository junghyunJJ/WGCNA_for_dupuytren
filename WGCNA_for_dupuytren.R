rm(list = ls())


library(WGCNA)
library(stringr)

load("WGCNA_GSE75152.rdata")

################################################################################
### select Soft Threshold (power) for WGCNA ####################################
################################################################################

powers <- c(c(1:11), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(t(wgcna_GSE75152), powerVector = powers,networkType = "signed")
# Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
# 1      1 0.075800 19.200          0.934  2070.0   2080.00   2130
# 2      2 0.000217 -0.204          0.954  1170.0   1170.00   1320
# 3      3 0.056400 -1.280          0.916   717.0    711.00    930
# 4      4 0.106000 -1.120          0.896   469.0    460.00    708
# 5      5 0.208000 -1.150          0.915   322.0    312.00    566
# 6      6 0.338000 -1.220          0.934   231.0    219.00    468
# 7      7 0.479000 -1.320          0.950   171.0    158.00    396
# 8      8 0.601000 -1.380          0.966   129.0    117.00    341
# 9      9 0.688000 -1.440          0.974   101.0     88.00    298
# 10    10 0.740000 -1.510          0.971    79.6     67.40    264
# 11    11 0.782000 -1.540          0.975    64.1     52.10    235
# 12    12 0.815000 -1.570          0.979    52.4     40.60    211 ***
# 13    14 0.861000 -1.590          0.983    36.3     25.50    174
# 14    16 0.879000 -1.630          0.977    26.2     16.50    146
# 15    18 0.897000 -1.630          0.980    19.5     11.10    124
# 16    20 0.910000 -1.620          0.983    15.0      7.71    107

par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"),cex.axis=0.8)
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=0.7,col="red")
abline(h=0.8,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"),cex.axis=0.8)
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.7,col="red")
par(mfrow = c(1,1))

################################################################################
#### WGCNA analysis ############################################################
################################################################################

adjacencyA1 <- adjacency(t(wgcna_GSE75152),power=12,type="signed")
diag(adjacencyA1) <- 0
dissTOMA1 <- 1-TOMsimilarity(adjacencyA1, TOMType="signed")
geneTree <- hclust(as.dist(dissTOMA1), method="average")

minModSize <- 50 # Modules are at least 50 genes large
dthresh <- 0.15 # MEs are no more than 0.85 correlated, if they are then the modules are merged and the ME is re-calculated
ds <- 2 # deep split parameter to determine how finely to cut the tree

tree1 <- cutreeHybrid(dendro = geneTree, pamStage=FALSE,
                      minClusterSize = minModSize,cutHeight = 0.99999,
                      deepSplit = ds, distM = dissTOMA1)
merged1 <- mergeCloseModules(exprData = t(wgcna_GSE75152),
                             colors = tree1$labels,
                             cutHeight = dthresh)
modules <- labels2colors(merged1$colors)


modules <- str_replace(modules,"^red$","imsi")
modules <- str_replace(modules,"^turquoise$","red")
modules <- str_replace(modules,"imsi","turquoise")

modules <- str_replace(modules,"^blue$","imsi")
modules <- str_replace(modules,"^greenyellow$","blue")
modules <- str_replace(modules,"imsi","greenyellow")

table(modules)

plotDendroAndColors(geneTree, modules,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
