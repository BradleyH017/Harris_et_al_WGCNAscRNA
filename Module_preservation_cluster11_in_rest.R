# Bradley Feb 2022
# Alternative assessment of trans-eQTL relatedness across clusters for Smillie et al scRNAseq
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/ModulePreservation/Tutorials/ TUTORIAL 2

# Set up
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/")
source('analysis.r')
library('Seurat')
library(ggsci)
library(rstatix)
library(optparse)
library(patchwork)
library(fastcluster)
library(dynamicTreeCut)
library(WGCNA)

# Set up

# generate the cluster list
clusters <- seq(1:12)
clusters <- clusters-1
clusters <- paste("cluster_", clusters, sep = "")

# Remove cluster 11 - This is the basis for all of the analysis
#clusters <- clusters[-which(clusters == "cluster_11")]

#~~~~~ 1.  Load the expression - This needs to be raw and then pseudo-bulked
# THINK: Using the pseudo-bulked 5000 HVG in analysis previously did not work, only 32 common genes - TRY THIS AGAIN
expr0 <- vector("list", length = length(clusters))
names(expr0) <- clusters


# Define pseudo bulk function 
tmm_PB <- function(object, assay, z_score, by){
                splitObj <- SplitObject(object, split.by = "Sample");
                sampexpr <- vector("list", length = length(splitObj));
                indx <- seq(1:length(splitObj));
                for(x in indx){sampexpr[[x]] = rowSums(splitObj[[x]][[assay]]@counts)};
                bulk <- data.frame(matrix(unlist(sampexpr), nrow=length(sampexpr), byrow=T));
                colnames(bulk) = names(sampexpr[[1]]);
                rownames(bulk) = names(splitObj);
                #apply scaling factor normalisation to each sample
                bulk <- as.data.frame(t(bulk))
                dge <- DGEList(counts = bulk);
                dge <- calcNormFactors(dge, method = "TMM");
                tmm <- edgeR::cpm(dge, normalized.lib.sizes = T,log = T);
                if(z_score == T){
                        genes <- rownames(tmm);
                        samples <- colnames(tmm)
                        if( by == "samples"){
                                z_sc <- apply(tmm, 2, FUN=scale, center = T, scale = T);
                                rownames(z_sc) <- genes;
                                return(z_sc)
                        }
                        if( by == "genes"){
                                z_sc <- apply(tmm, 1, FUN=scale, center = T, scale = T);
                                rownames(z_sc) <- samples;
                                z_sc <- data.frame(t(z_sc));
                                return(z_sc)
                        }
                } else {
                        return(tmm)
                }
        };
 
# Load object
pathOut = c("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/BH_analysis/Seurat/August_2021/Epithelial")
load(paste(pathOut, "objects/Processed_seur.integrated.clusteredRds", sep = "/"))


# Define seurlist and use this to extract expression for genes in each cluster. Pseudo-bulk as have done before
seurlist <- vector("list", length = length(clusters))
Idents(seur.integrated) <- "seurat_clusters"
cluster_index <- gsub("cluster_", "", clusters)
for(cluster in seq_along(clusters)){
	seurlist[[cluster]] <- subset(x = seur.integrated, idents=as.character(cluster_index[cluster]))
	expr0[[cluster]] <- vector("list", length = 1)
	names(expr0[[cluster]]) <- "data"
	expr0[[cluster]]$data <- tmm_PB(object = seurlist[[cluster]], assay = "RNA", z_score = T, by = "genes")
	# Exclude C11orf53, COLCA1, COLCA2
#	expr0[[cluster]]$data <- expr0[[cluster]]$data[-which(rownames(expr0[[cluster]]$data) == "C11orf53"),]
#	expr0[[cluster]]$data <- expr0[[cluster]]$data[-which(rownames(expr0[[cluster]]$data) == "COLCA1"),]
#	expr0[[cluster]]$data <- expr0[[cluster]]$data[-which(rownames(expr0[[cluster]]$data) == "COLCA2"),]
}

# ~~~~~~ 2. Restrict to common genes
# Have used all genes so can skip this


# ~~~~~~ 3. Loading genes that belog to the selected interseting module from cluster 11
gI <- read.csv("BH_analysis/Seurat/August_2021/Epithelial/WGCNA_within_clusters/cluster_11/tables/geneInfo.csv", row.names=1)
cluster11_black <- rownames(gI[gI$moduleColor == "black",])
# Could subset these here for trans-eQTL targets in this cluster?
trans <- c("LRMP", "SH2D6", "PSTPIP2", "HTR3E", "TRPM5", "HTR3C", "ALOX5", "OGDHL", "BMX", "MATK", "SH2D7", "PIK3CG", "PLCG2", "PTGS1", "IL17RB", "AZGP1", "GNG13", "CAMP", "ANKHD1", "GIN1", "SPAG6")
cluster11_black <- cluster11_black[cluster11_black %in% trans] 
# Add C11orf53
cluster11_black <- c(cluster11_black, c("C11orf53"))


# ~~~~~ 4. Good Ssample Genes
for(cluster in seq_along(clusters)){
	# Transpose
	expr0[[cluster]]$data <- as.data.frame(t(expr0[[cluster]]$data))
}
ggs = goodSamplesGenesMS(expr0) 
for(set in seq_along(clusters)){
	expr0[[set]]$data = expr0[[set]]$data[, ggs$goodGenes]
}
collectGarbage();

keep <- colnames(expr0[[1]]$data)[ggs$goodGenes]
cluster11_black <- cluster11_black[cluster11_black %in% keep]
print(paste("The final list for cluster11_black is", cluster11_black))

# ~~~~~ 5. Set up module labels that reflect the pathway membership
size <- checkSets(expr0)
nGenes = size$nGenes
pathGenes = cluster11_black
nPathGenes = length(pathGenes)
pathwayLabels = sample(c(0,1), nGenes, replace = TRUE)
names(pathwayLabels) <- colnames(expr0[[1]]$data)
pathwayLabels[cluster11_black] = 2


# Calculate module eigengenes of the pathway in each cluster
multiMEs = multiSetMEs(expr0, universalColors = pathwayLabels) 

KMEpathway = matrix(0, nPathGenes, length(clusters))

for (set in seq_along(clusters)) { 
	KMEpathway[, set] = bicor(expr0[[set]]$data[, pathGenes], multiMEs[[set]]$data[, 3], use = "p") 
} # If necessary, can plot a scatterplot of the module membeships if (FALSE) { sizeGrWindow(9,9); verbosePairs(KMEpathway, names = setNames, abline = TRUE) }

# If necessary, can plot a scatterplot of the module membeships 
#sizeGrWindow(9,9); verbosePairs(KMEpathway, names = setNames, abline = TRUE)

# Calculate module Preservation
multiColor = list()
for (set in seq_along(clusters)) {
	multiColor[[set]] = pathwayLabels
}
names(multiColor) = clusters
mp = modulePreservation(expr0, multiColor, networkType = "unsigned", corFnc = "bicor", referenceNetworks = c(1:length(clusters)), verbose = 4, maxModuleSize = 500, maxGoldModuleSize = 500, nPermutations = 200)
mp_signed = modulePreservation(expr0, multiColor, networkType = "signed", corFnc = "bicor", referenceNetworks = c(1:length(clusters)), verbose = 4, maxModuleSize = 500, maxGoldModuleSize = 500, nPermutations = 200)
save.image(file = paste(pathOut, "objects/MPupToindividualPathways-allgenes-mp.RData", sep = "/"));


#~~~~ 6. Calculate in group proportion
library(clusterRepro) 
cr = list()
set.seed(10)
for(set in seq_along(clusters)){
	printFlush(paste("Working on ", clusters[set]))
	# Leave out
	centr = as.matrix(multiMEs[[set]]$data[, c(2:3)])
	rownames(centr) = rownames(expr0[[set]]$data)
	colnames(centr) = c("Random", "cluster11_black")
	print(system.time({ cr[[set]] = clusterRepro(Centroids = centr, New.data = expr0[[set]]$data, Number.of.permutations = 2000); } ))
	save(cr, file = paste(pathOut, "objects/individualPathways-allgenes-cr.RData", sep = "/")) }
}



#~~~~~~ 7. Analysis and graphical representation of results - Doing this for the signed network first
nSets <- length(clusters)
Zsummary = matrix(NA, nSets, nSets)
Zdensity = matrix(NA, nSets, nSets)
Zconn = matrix(NA, nSets, nSets)
Zquality = matrix(NA, nSets, nSets)
for (ref in 1:nSets){
	for (test in 1:nSets){
		if (ref!=test) { 
			Z = c(as.matrix(mp_signed$preservation$Z[[ref]][[test]][-c(1:3), -1]))
			Zsummary[ref, test] = Z[1]; Zdensity[ref, test] = Z[2]
			Zconn[ref, test] = Z[3]; Zquality[ref, test] = median(c(as.matrix(mp_signed $quality$Z[[ref]][[test]][4, 2])))
		}
	}
}
Zs = list()
Zs[[1]] = Zquality
Zs[[2]] = Zsummary
Zs[[3]] = Zdensity
Zs[[4]] = Zconn
Znames = c("Zsummary.quality", "Zsummary.preservation", "Zdensity.preservation", "Zconnectivity.preservation")
setNames <- clusters

sizeGrWindow(10, 8)
layout(matrix(c(1:16), 4, 4, byrow = TRUE))
par(mar = c(5.3, 3.3, 4, 0.5))
par(mgp = c(1.9, 0.6, 0)) 
# Plot of summary Z statistics
ref=12
for (i in 2:4) { 
	pres = Zs[[i]][ref,-ref]
	labeledBarplot(pres, names = setNames[-ref], main = spaste(LETTERS[i-1], ". ", Znames[i], "\nReference set: cluster_11"), cex.names = 1.2, cex.axis = 1.2, ylab = Znames[i], cex.main = 1.4, cex.lab = 1.2, labels = setNames) 
}

# Plot of IGP calculated by clusterRepro 
labeledBarplot2(crPathway, names = setNames, main = "D. Observed IGP", ylab = "Actual.IGP", cex.names = 1.2, cex.axis = 1.2, cex.main = 1.4, cex.lab = 1.2)

# KME scatterplots - Again, cluster 11 (index 12) as the reference
par(mar = c(3.3, 3.3, 4, 0.5))
par(mgp = c(1.9, 0.6, 0))
ref = 12
ind = 5 
for (set in 1:nSets) {
	if (set!=ref) { 
		verboseScatterplot(KMEpathway[, ref], KMEpathway[, set], xlab = spaste("KME in ", setNames[ref]), ylab = spaste("KME in ", setNames[set]), main = spaste(LETTERS[ind], ". KME in ", setNames[set], "\nvs. ", setNames[ref], "\n"), cex.main = 1.4, cex.lab = 1.2, cex.axis = 1.2, abline = TRUE ) 
		ind = ind + 1; 
	} else {
		plot(c(0,1), type ="n", axes = FALSE, xlab = "", ylab ="")
	}
}

# Save all the results to a file
ref = 12
stats = list()
ind = 1
refNames = NULL
testNames = NULL
dropRows = c(1,3)
for (test in 1:nSets) {
	if (ref!=test) {
		stats[[ind]] = cbind(mp_signed$quality$observed[[ref]][[test]][-dropRows,, drop = FALSE],
		mp_signed$preservation$observed[[ref]][[test]][-dropRows, -1, drop = FALSE],
		mp_signed$referenceSeparability$observed[[ref]][[test]][-dropRows, -1, drop = FALSE],
		mp_signed$testSeparability$observed[[ref]][[test]][-dropRows, -1, drop = FALSE], 
		mp_signed$quality$Z[[ref]][[test]][-dropRows, -1, drop = FALSE], 
		mp_signed$quality$log.p[[ref]][[test]][-dropRows, -1, drop = FALSE], 
		mp_signed$quality$log.pBonf[[ref]][[test]][-dropRows, -1, drop = FALSE], 
		mp_signed$preservation$Z[[ref]][[test]][-dropRows,-1, drop = FALSE], 
		mp_signed$preservation$log.p[[ref]][[test]][-dropRows,-1, drop = FALSE], 
		mp_signed$preservation$log.pBonf[[ref]][[test]][-dropRows,-1, drop = FALSE], 
		mp_signed$referenceSeparability$Z[[ref]][[test]][-dropRows, -1, drop = FALSE], 
		mp_signed$referenceSeparability$log.p[[ref]][[test]][-dropRows, -1, drop = FALSE], 
		mp_signed$referenceSeparability$log.pBonf[[ref]][[test]][-dropRows, -1, drop = FALSE], 
		mp_signed$testSeparability$Z[[ref]][[test]][-dropRows, -1, drop = FALSE], 
		mp_signed$testSeparability$log.p[[ref]][[test]][-dropRows, -1, drop = FALSE], 
		mp_signed$testSeparability$log.pBonf[[ref]][[test]][-dropRows, -1, drop = FALSE]);
		ind = ind + 1
		refNames = c(refNames, setNames[ref])
		testNames = c(testNames, setNames[test])
		}
	}

table = NULL
nPlots = ind-1
for (p in 1:nPlots) {
	nMods = nrow(stats[[p]])
	modNames = rownames(stats[[p]])
	modStatus = rep("Cluster11_black", nMods)
	modStatus[modNames=="grey"] = "improper module" 
	modStatus[modNames=="0.1"] = "all-network sample" 
	description = cbind( referenceSet = rep(refNames[p], nMods), testSet = rep(testNames[p], nMods), module = rownames(stats[[p]]), moduleType = modStatus)
	if (is.null(table)) {
		table = cbind(description, stats[[p]]) 
		} else {
			table = rbind(table, cbind(description, stats[[p]])); 
		} 
}
x = as.matrix(table)
write.table(table, file = paste(pathOut, "MP_cluster11_in_others/allResults.csv", sep = "/"), sep = ",", quote = FALSE, row.names = FALSE);




#~~~~~~~ 8. Custom circle plot
#===============================================================================================
#
# grey2red
#
#===============================================================================================

grey2red = function(n, base, gamma)
{
  red = seq(from=base^gamma, to=1, length.out = n)^(1/gamma)
  green = blue = seq(from = base^gamma, to=0, length.out = n)^(1/gamma);
  col = rgb(red, green, blue, maxColorValue = 1); 
}


# Example of grey2red:

if (FALSE)
{
  par(mfrow = c(5,1))
  par(mar = c(1,3,1,1))
  n= 100
  barplot(rep(1, n), col = grey2red(n, 0, 1))
  barplot(rep(1, n), col = grey2red(n, 1, 1))
  barplot(rep(1, n), col = grey2red(n, 0.5, 1))
  barplot(rep(1, n), col = grey2red(n, 0.5, 0.2))
  barplot(rep(1, n), col = grey2red(n, 0.5, 5.0))
}



#===============================================================================================
#
# Circle plot for generating VisANT-like plots
#
#===============================================================================================

circlePlot = function(
  adjacency,
  labels,
  order,
  startNewPlot = TRUE,
  plotBox = c(-1, 1, -1, 1),
  center = c(0,0), 
  radii = c(0.8, 0.8),
  startAngle = 0,
  variable.cex.labels = TRUE,
  min.cex.labels = 1,
  max.cex.labels = 1.5,
  variable.cex.points = TRUE,
  min.cex.points = 1,
  max.cex.points = 3,
  variable.line.width = TRUE,
  min.line.width = 1,
  max.line.width = 5,
  lineColors = grey2red(50, 0.6, 1),
  pch = 21,
  labelColors = "black",
  pointColors = "black",
  pointBg = "black",
  xMargin = 1-radii[1],
  yMargin = 1-radii[2],
  xLabelOffset = 0.01,
  yLabelOffset = 0.01,
  variableLabelAngle = TRUE,
  
  ...)
{

  if (startNewPlot)
    plot(plotBox[1:2], plotBox[3:4], axes = FALSE, type = "n", xlab = "", ylab = "", ...);

  # plot(c(-1-xMargin,1+xMargin), c(-1-yMargin,1+yMargin), axes = FALSE, type = "n", xlab = "", ylab = "", ...) 
  checkAdjMat(adjacency, min = -1)
  n = length(labels);
  angles = seq(from = startAngle, to = startAngle + 2*pi * (1-1/n), length.out = n);
  x = center[1] + radii[1] * sin(angles);  # This is intentional; top should correspond to angle=0
  y = center[2] + radii[2] * cos(angles);

  adjx = adjacency
  adjx[is.na(adjx)] = 0;
  connectivity = apply(abs(adjx), 2, sum)-diag(adjx)
  minConn = min(connectivity, na.rm = TRUE);
  maxConn = max(connectivity, na.rm = TRUE);

  if (length(pch)==1) pch = rep(pch, n);
  if (length(labelColors)==1) labelColors = rep(labelColors, n);
  if (length(pointColors)==1) pointColors = rep(pointColors, n);
  if (length(pointBg)==1) pointBg = rep(pointBg, n);
  if (length(xLabelOffset)==1) xLabelOffset = rep(xLabelOffset, n);
  if (length(yLabelOffset)==1) yLabelOffset = rep(yLabelOffset, n);

  oLabs = labels[order]
  oLColors = labelColors[order];
  oPColors = pointColors[order];
  oPBg = pointBg[order];
  oConn = connectivity[order];
  oAdj = adjx[order, order];
  oPch = pch[order];

  actualCexPts = rep(0, n);
  for (node in 1:n)
  {
    cex = min.cex.points;
    if (variable.cex.points)
       cex = min.cex.points + (max.cex.points - min.cex.points) * 
                                  (oConn[node] - minConn)/(maxConn - minConn)
    actualCexPts[node] = cex
  }

  diag(oAdj) = 0;
  maxA = max(abs(oAdj));
  if (sum(oAdj < 0) > 0)
  {
     adjCol = numbers2colors(oAdj, signed = TRUE, lim = c(-maxA, maxA));
  } else {
     adjCol = numbers2colors(oAdj, signed = FALSE, lim = c(0, maxA));
  }


  ltA = oAdj;
  diag(ltA) = NA;
  ltA[upper.tri(ltA)] = NA;

  adjOrder = order(c(abs(ltA)))
  rows = row(oAdj)[adjOrder];
  cols = col(oAdj)[adjOrder];

  nLines = n*(n-1)/2;
  for (line in 1:nLines)
  {
    n1 = rows[line];
    n2 = cols[line];
    a = oAdj[n1, n2];
    normA = abs(a)/maxA;

    w = min.line.width;
    if (variable.line.width)
      w = min.line.width + (max.line.width - min.line.width) * normA;

    #pRadius1 = par("cxy") * actualCexPts[n1]/35;  # Emprical fudge factor..
    #pRadius2 = par("cxy") * actualCexPts[n2]/35;
    lineLen = sqrt( (x[n1] - x[n2])^2 + (y[n1] - y[n2])^2);
    x1 = x[n1] #+ pRadius1[1] * (x[n2] - x[n1]) / lineLen
    y1 = y[n1] #+ pRadius1[1] * (y[n2] - y[n1]) / lineLen
    x2 = x[n2] #+ pRadius2[1] * (x[n1] - x[n2]) / lineLen
    y2 = y[n2] #+ pRadius2[1] * (y[n1] - y[n2]) / lineLen

    lines(c(x1,x2),c(y1, y2), lwd = w, col = adjCol[n1, n2]);
  }

  for (node in 1:n)
    points(x[node], y[node], pch = oPch[node], cex = actualCexPts[node], bg = oPBg[node], col = oPColors[node]);

  for (node in 1:n)
  {
    cex = min.cex.labels;
    if (variable.cex.labels)
       cex = min.cex.labels + (max.cex.labels - min.cex.labels) *
                                 (oConn[node] - minConn)/(maxConn - minConn)
    textWidth = strwidth(oLabs[node], cex = cex);
    textHeight = strheight(oLabs[node], cex = cex);
    if (variableLabelAngle)
    {
      ang = angles[node]/pi * 180;
      if (ang < 180) 
      {
         dir = 1;
      } else {
         dir = -1;
         ang = ang - 180;
      }
      ang = (90 - ang)/2
      xDir = 1;
      yDir = 1;
      cosAng = cos(ang/180*pi);
      sinAng = sin(ang/180*pi);
    } else {
      ang = 0;
      xDir = x[node];
      yDir = y[node];
      cosAng = 1;
      sinAng = 1;
      dir = 1;
    }
    angRad = ang/180*pi;
    pRadius = par("cxy") * actualCexPts[node]/5  ;  # Emprical fudge factor..
    effPointRadius = sqrt(sum(c(cosAng^2, sinAng^2) * pRadius^2));
    rotMat = matrix( c(cosAng, sinAng, -sinAng, cosAng), 2, 2);
    labelShift = rotMat %*% as.matrix(c(textWidth, textHeight));
    text(x[node] + dir * xDir * (labelShift[1]/2 + cosAng * effPointRadius + xLabelOffset[node]), 
         y[node] + dir * yDir * (labelShift[2]/2 + sinAng * effPointRadius + yLabelOffset[node]), 
         labels = oLabs[node], adj = c(0.5, 0.5), 
         cex = cex, col = oLColors[node], srt = ang, xpd = TRUE);
  }

}

# Example of circle plot:

if (FALSE)
{

   sizeGrWindow(8,8)
   par(mfrow = c(1,1));
   nS = 100;
   nn = 30;
   mod = simulateModule(rnorm(nS), nn);
   adjacency = abs(cor(mod))^3;
   
   order = c(1:nn);
   
   labels = paste("Gene", c(1:nn));
   circlePlot(adjacency, labels, order, variable.cex.labels = FALSE, radii = c(0.6, 0.6));
   # Fix the overlaping labels - for now this requires semi-manual intervention.
   xOffset = rep(0.01, nn);
   xOffset[16] = -strwidth(labels[16])/3;
   circlePlot(adjacency, labels, order, xLabelOffset = xOffset);


   # Plot two circles in one plot

   circlePlot(adjacency, labels, order, variable.cex.labels = FALSE, center = c(-0.5, -0.5), 
              radii = c(0.35, 0.35));

   circlePlot(adjacency, labels, order, startNewPlot = FALSE, 
              variable.cex.labels = FALSE, center = c(0.5, 0.5), 
              radii = c(0.35, 0.35));

   
}

# Calculate adjacencies within the module
pathwayAdjs = list()
for (set in 1:nSets) { 
	printFlush(paste("Working on set", setNames[set]))
	bc = bicor(expr0[[set]]$data[, pathGenes], use = "p")
	#bc[bc<0] = 0
	pathwayAdjs[[set]] = abs(bc)^4 * sign(bc)
}

# We order the genes by a weighted average connectivity
conn = matrix(0, nPathGenes, nSets)
for (set in 1:nSets) {
	conn[, set] = apply(abs(pathwayAdjs[[set]]), 2, sum)-1
}
weights = c(3,1,5,1, 3,1,5,1)
wMat = matrix(weights, nPathGenes, nSets, byrow = TRUE)
wconn = apply(conn * wMat, 1, sum)
order = order(-wconn)
# use the gene names as lables
labels = colnames(pathwayAdjs[[ref]])
sizeGrWindow(12,6)
par(mfrow =c(2,4))
par(mar = c(0.3, 0.2, 1.5, 0.2)) 
for (set in 1:nSets) {
	circlePlot(pathwayAdjs[[set]], labels, order, main = setNames[set],
		variable.cex.labels = TRUE,
		radii = c(0.56, 0.62),
		center = c(0.1, 0.04),
		min.cex.labels = 1.2,
		max.cex.labels = 1.4,
		cex.main = 1.4);
}









#### SIMPLE EXPLORATION
# looking at the pairwise correlations between cluster 11 black genes  + C11orf53 with one another in each cluster
cluster11_black <- rownames(gI[gI$moduleColor == "black",])
# Could subset these here for trans-eQTL targets in this cluster?
trans <- c("LRMP", "SH2D6", "PSTPIP2", "HTR3E", "TRPM5", "HTR3C", "ALOX5", "OGDHL", "BMX", "MATK", "SH2D7", "PIK3CG", "PLCG2", "PTGS1", "IL17RB", "AZGP1", "GNG13", "CAMP", "ANKHD1", "GIN1", "SPAG6")
cluster11_black <- cluster11_black[cluster11_black %in% trans] 
# Add C11orf53
cluster11_black <- c(cluster11_black, c("C11orf53"))

expr0_trans <- vector("list", length = length(clusters))
cors <- vector("list", length = length(clusters))
ps <- vector("list", length = length(clusters))
names(expr0_trans) <- clusters
names(cors) <- clusters
names(ps) <- clusters
for(cluster in seq_along(clusters)){
	expr0_trans[[cluster]] <- expr0[[cluster]]$data[,colnames(expr0[[cluster]]$data) %in% cluster11_black]
	cors[[cluster]] <- matrix(NA, length(cluster11_black), length(cluster11_black))
	rownames(cors[[cluster]]) <- cluster11_black
	colnames(cors[[cluster]]) <- cluster11_black
	rownames(ps[[cluster]]) <- cluster11_black
	colnames(ps[[cluster]]) <- cluster11_black
	for(c in seq_along(cluster11_black)){
		for(r in seq_along(cluster11_black)){
			cors[[cluster]] <- cor.test(expr0_trans[[c]], expr0_trans[[r]], )
		}
	}
}