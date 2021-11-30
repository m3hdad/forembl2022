library(WGCNA)
library(Biobase)
library(RColorBrewer)
library(tidyverse)
library(ggplot2)
library(gplots)
library(reshape2)
library(plotly)

options(StringsAsFactors = FALSE)
load("obj/eset.RData")

tissuetype <- factor(pData(eset)$tissue)
ntissuetype <- nlevels(tissuetype)
tissuecols <- colorRampPalette(brewer.pal(ntissuetype, "Set1"))(ntissuetype)
coltissuetype <- tissuecols[tissuetype]

#raw density
esetlog2 <- eset
exprs(esetlog2) <- log2(exprs(eset))
x <- melt(as.matrix(exprs(esetlog2)))
colnames(x) <- c("geneid","sample","value")
x$experiment <- substr(x$sample,1,nchar(as.character(x$sample)) - 4)
x <- left_join(x,pData(esetlog2),by = "experiment")
g <- ggplot(x, aes(x = value, color = tissue)) +
  geom_density()
#  facet_wrap(~tissue)

ggplotly(g)
#Removed 3161671 rows containing non-finite values (stat_density). 
# diff when preprocess

# Raw sample clustering
pdf(file = "plots/heatmapsampleClusteringraw.pdf", width = 11.69, height = 8.27)
#sizeGrWindow(12,9)
#par(cex = 1)
par(mar = c(3,1,9,1))
heatmap.2(cor(exprs(eset)),
          RowSideColors = coltissuetype,
          cexRow = 0.5,
          cexCol = 0.5,
          trace = "none", 
          keysize = 1,
          main = "Sample Correlations (raw)")

dev.off()

#Manual cleanups
load("obj/eset.RData")
table(pData(eset)$tissue)

low_count_mask <- apply(exprs(eset), 1, function(x){sum(x < 0.005) > length(x) * 0.95})
eset95 <- eset[!low_count_mask,]
dim(eset95)

sprintf("Removing %d low-count genes (%d remaining).",
        sum(low_count_mask), 
        sum(!low_count_mask))

# Log trans
eset95log2 <- eset95
exprs(eset95log2) <- log2(exprs(eset95) + 1)

# 95%log2 sample clustering
tissuetype <- factor(pData(eset95log2)$tissue)
ntissuetype <- nlevels(tissuetype)
tissuecols <- colorRampPalette(brewer.pal(ntissuetype, "Set1"))(ntissuetype)
coltissuetype <- tissuecols[tissuetype]

pdf(file = "plots/heatmapsampleClustering95log2.pdf", width = 11.69, height = 8.27)
#sizeGrWindow(12,9)
#par(cex = 1)
par(mar = c(3,1,9,1))
heatmap.2(cor(exprs(eset95log2)),
          RowSideColors = coltissuetype,
          cexRow = 0.5,
          cexCol = 0.5,
          trace = "none", 
          keysize = 1,
          main = "Sample Correlations 95% log2")

dev.off()

# Density plot
x <- melt(as.matrix(exprs(eset95log2)))
colnames(x) <- c("geneid","sample","value")
x$experiment <- substr(x$sample,1,nchar(as.character(x$sample)) - 4)
x <- left_join(x,pData(eset95log2),by = "experiment")
g <- ggplot(x, aes(x = value, color = tissue)) +
  geom_density()
#  facet_wrap(~tissue)

ggplotly(g)

dim(eset95log2)

# Remove non diff genes

low_var_mask <- apply(exprs(eset95log2), 1, var) > 0
table(low_var_mask)

dim(fData(eset95log2))
fnames <- colnames(fData(eset95log2))
length(fnames)
deval <- fData(eset95log2)[,8:length(fnames)]
nDEG <- rowSums(is.na(deval))
hist(nDEG)
head(nDEG)
x <- nDEG == length(fnames)
x <- nDEG == length(fnames) - 6
table(x)



# WGCNA cleanups
load("obj/eset.RData")
tdat <- t(exprs(eset))
gsg <- goodSamplesGenes(tdat,verbose = 5)
gsg$allOK

if (!gsg$allOK) {
  # Optionally print the genes and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(colnames(tdat)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(tdat)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  geset <- eset[gsg$goodGenes, gsg$goodSamples]
}
#clean
rm(tdat)

# WGCNA cleanups of 95log2
tdat <- t(exprs(eset95log2))
gsg <- goodSamplesGenes(tdat,verbose = 5)
gsg$allOK

if (!gsg$allOK) {
  # Optionally print the genes and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(colnames(tdat)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(tdat)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  geset95log2 <- eset95log2[gsg$goodGenes, gsg$goodSamples]
}
#clean
rm(tdat)


# dim of
dim(eset)
dim(geset)
dim(eset95log2)

# Sample clustering WGCNA cleanups
tissuetype <- factor(pData(eset95log2)$tissue)
ntissuetype <- nlevels(tissuetype)
tissuecols <- brewer.pal(n = ntissuetype, name = "Set1")
coltissuetype <- tissuecols[tissuetype]

pdf(file = "plots/heatmapsampleClusteringgeset.pdf", width = 11.69, height = 8.27)
#sizeGrWindow(12,9)
#par(cex = 1)
par(mar = c(3,1,9,1))
heatmap.2(cor(exprs(geset)),
          RowSideColors = coltissuetype,
          cexRow = 0.5,
          cexCol = 0.5,
          trace = "none", 
          keysize = 1,
          main = "Sample Correlations (geset)")

dev.off()

#clustering samples
sampleTree <- hclust(dist(t(exprs(eset95log2))), method = "average")

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
require(dendextend)
dend <- as.dendrogram(sampleTree)

#Get tissue types and numbers
tissuetype <- factor(pData(eset95log2)$tissue)
ntissuetype <- nlevels(tissuetype)
#Set colors for tissues and assign
tissuecols <- brewer.pal(n = ntissuetype, name = "Set1")
coltissuetype <- tissuecols[tissuetype]

#Update attributes of dendrogram (leaves color)
#labels_colors(dend) <- coltissuetype[order.dendrogram(dend)]

#Output
#sizeGrWindow(12,9)
pdf(file = "plots/sampleClustering.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(14,7,7,7))
plot(dend)
colored_bars(colors = coltissuetype, dend = dend, rowLabels = "Tissue")
legend("topright", legend = levels(tissuetype), fill = tissuecols)
dev.off()

#clean
rm(coltissuetype,ntissuetype,tissuetype,tissuecols)

# Choose a set of soft-thresholding powers
#powers <- c(c(1:10), seq(from = 12, to=20, by=2))
powers <- c(1:12,seq(13,23,2))
# Call the network topology analysis function
sft <- pickSoftThreshold(t(exprs(eset95log2)),
                         powerVector = powers,
                         verbose = 5,
                         nBreaks = 10,
                         blockSize = length(rownames(eset95log2)))
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 <- 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2^",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


net <- blockwiseModules(# Input data
                        datExpr = t(exprs(eset95log2)),
                        weights = NULL,
                        # Data checking options
                        checkMissingData = T,
                        # Options for splitting data in blocks
                        blocks = NULL,
                        maxBlockSize = length(rownames(eset95log2)),
                        # Adjacency function options
                        power = 9,
                        networkType = "unsigned",
                        # Topological overlap options
                        TOMType = "unsigned",
                        TOMDenom = "min",
                        suppressTOMForZeroAdjacencies = F,
                        suppressNegativeTOM = F,
                        # Saving or returning TOM
                        saveTOMs = TRUE,
                        saveTOMFileBase = "obj/physcoTOM", 
                        # Basic Tree cut options
                        minModuleSize = min(20, ncol(t(exprs(eset95log2))) / 2),
                        # Advance Tree cut options
                        pamRespectsDendro = FALSE,
                        # Module merging options
                        mergeCutHeight = 0.25,
                        impute = T,
                        trapErrors = F,
                        # Gene reassignment, module trimming, and module "significance" criteria
                        reassignThreshold = 0, 
                        # Output options
                        numericLabels = TRUE,
                        # Options controlling behaviour 
                        verbose = 3)
save(net, file = "./obj/netauto.RData")
load("./obj/netauto.RData")
gc()
table(net$colors)
class(net)
mod12 <- exprs(eset95log2)[net$colors == 12,]
mod12
dim(mod12)




# post ana
load("obj/netauto.RData")
load("obj/physco-networkConstruction-auto.RData")
load("obj/physcoTOM-block.1.RData")



# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors <- labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "obj/physco-networkConstruction-auto.RData")

# Define numbers of genes and samples
nGenes <- nrow(exprs(eset95log2))
nSamples <- ncol(exprs(eset95log2))
# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(t(exprs(eset95log2)), moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, pData(eset95log2), use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(pData(geset)),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))













#################################################################### 
####################################################################
####################################################################
####################################################################
##########               STEP BY STEP                  #############
####################################################################
####################################################################
####################################################################
####################################################################

enableWGCNAThreads()


powers <- c(1:12,seq(12,25,2))
sft <- pickSoftThreshold(t(exprs(eset95log2)), powerVector = powers, verbose = 5,blockSize = length(rownames(eset95log2)))

# Plot the results:


powers <- c(1:12,seq(13,23,2))
# Call the network topology analysis function
sft <- pickSoftThreshold(t(exprs(eset95log2)),
                         powerVector = powers,
                         verbose = 5,
                         nBreaks = 10,
                         blockSize = length(rownames(eset95log2)))
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 <- 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2^",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower <- 9
adjacency <- adjacency(#col correspond to genes, rows to samples
                       datExpr = t(exprs(eset95log2)),
                       selectCols = NULL,
                       # network type: "unsigned", "signed", "signed hybrid", "distance"
                       type = "unsigned",
                       power = softPower,
                       # function to be used to calculate co-expression similarity for correlation networks,
                       # default: pearson, any func returning values -1 to 1 can be used!
                       corFnc = "cor", corOptions = list(use = "p"),
                       weights = NULL,
                       distFnc = "dist", distOptions = "method = 'euclidean'",
                       weightArgNames = c("weights.x", "weights.y")
                       )


# Turn adjacency into topological overlap
TOM <- TOMsimilarity(adjacency)


dissTOM <- 1-TOM


# Call the hierarchical clustering function
geneTree <- hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram)
#sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high:
minModuleSize <- 30
# Module identification using dynamic tree cut:
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)

table(dynamicMods)



# Convert numeric lables into colors
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList <- moduleEigengenes(t(exprs(eset95log2)), colors = dynamicColors)
MEs <- MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss <- 1-cor(MEs)
# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")
# Plot the result
#sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres <- 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")

# Call an automatic merging function
merge <- mergeCloseModules(t(exprs(eset95log2)), dynamicColors, cutHeight = MEDissThres, verbose = 3)

      
# The merged module colors
mergedColors <- merge$colors
# Eigengenes of the new merged modules:
mergedMEs <- merge$newMEs

sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

dev.off()

# Rename to moduleColors
moduleColors <- mergedColors
# Construct numerical labels corresponding to the colors
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder)-1
MEs <- mergedMEs
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "../objects/FemaleLiver-02-networkConstruction-stepByStep.RData")