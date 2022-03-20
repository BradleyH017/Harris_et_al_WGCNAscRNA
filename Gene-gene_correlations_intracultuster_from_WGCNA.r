#### Bradley Jan 2022
# Looking at the correlation of z-scores between Peter trans-eQTL targets from WGCNA within Smillie clusters
# Simplifying the anaysis done so far (see WGCNA_within_clusters_summary_Peters.r)

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
set <- "Epithelial"
pathOut = paste("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/BH_analysis/Seurat/August_2021", set, "WGCNA_within_clusters", sep = "/")


# generate the cluster list
clusters <- seq(1:12)
clusters <- clusters-1
clusters <- paste("cluster_", clusters, sep = "")


# Load the expression data used in each analysis, extract that of the trans-eQTL targets only + C11orf53 (C11orf53 expression was removed to look for e=correlations with in the WGCNA. So need to re-append this onto the expression)
want <- v("C11orf53", "LRMP", "SH2D6", "HTR3E", "ALOX5", "OGDHL", "BMX", "MATK", "SH2D7", "PTGS1", "IL17RB", "AZGP1", "GNG13")
expr_list <- vector("list", length = length(clusters))
names(expr_list) <- clusters
for(cluster in seq_along(clusters)){
	pathOut = paste("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/BH_analysis/Seurat/August_2021", set, "WGCNA_within_clusters", sep = "/")
	load(paste(pathOut, clusters[cluster], "objects", "UpToSftThreshold.Rds", sep = "/"))
	expr_list[[cluster]] <- exprdata[rownames(exprdata) %in% want,]
	expr_list[[cluster]] <- as.data.frame(t(expr_list[[cluster]]))
	expr_list[[cluster]]$C11orf53 <- meta$C11orf53
}

# Calculate the pairwise correlations for each cluster - only use if the pvalue is < 0.05
cors <- vector("list", length = length(clusters))
names(cors) <- clusters
ps <- vector("list", length = length(clusters))
names(ps) <- clusters
for(cluster in seq_along(clusters)){
	cors[[cluster]] <- matrix(0, nrow = length(want), ncol = length(want))
	ps[[cluster]] <- matrix(0, nrow = length(want), ncol = length(want))
	use <- colnames(expr_list[[cluster]])
	rownames(cors[[cluster]]) <- want
	colnames(cors[[cluster]]) <- want
	for(gene1 in seq_along(use)){
		for(gene2 in seq_along(use)){
			temp <-  cor.test(
				expr_list[[cluster]][,which(colnames(expr_list[[cluster]]) == use[gene1])],
				expr_list[[cluster]][,which(colnames(expr_list[[cluster]]) == use[gene2])],
				method = "pearson",
				use = "complete.obs")
			ps[[cluster]][which(want == use[gene1]), which(want == use[gene2])] <- temp$p.value
			if(temp$p.value < 0.05)	{
				cors[[cluster]][which(want == use[gene1]), which(want == use[gene2])] <- temp$estimate
			}
		}
	}
}

# Plot this
my_fun <- function(cc, tt){
	corrplot::corrplot(cc, , is.corr = T , type = "upper",diag = F, na.label = "NA", na.label.col = "white", title = tt, mar=c(0,0,0.8,0))
}


plot_list <- vector("list", length = length(clusters))
names(plot_list) <- clusters
par(mfrow=c(3,4))


plot_list[[1]] <- my_fun(cc = cors[[1]], tt = clusters[1])
plot_list[[2]] <- my_fun(cc = cors[[2]], tt = clusters[2])
plot_list[[3]] <- my_fun(cc = cors[[3]], tt = clusters[3])
plot_list[[4]] <- my_fun(cc = cors[[4]], tt = clusters[4])
plot_list[[5]] <- my_fun(cc = cors[[5]], tt = clusters[5])
plot_list[[6]] <- my_fun(cc = cors[[6]], tt = clusters[6])
plot_list[[7]] <- my_fun(cc = cors[[7]], tt = clusters[7])
plot_list[[8]] <- my_fun(cc = cors[[8]], tt = clusters[8])
plot_list[[9]] <- my_fun(cc = cors[[9]], tt = clusters[9])
plot_list[[10]] <- my_fun(cc = cors[[10]], tt = clusters[10])
plot_list[[11]] <- my_fun(cc = cors[[11]], tt = clusters[11])
plot_list[[12]] <- my_fun(cc = cors[[12]], tt = clusters[12])



# 2.  Using WGCNA circle plot representation

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

for(cluster in seq_along(clusters)){
	cors[[cluster]] <- matrix(NA, nrow = length(want), ncol = length(want))
	ps[[cluster]] <- matrix(NA, nrow = length(want), ncol = length(want))
	use <- colnames(expr_list[[cluster]])
	rownames(cors[[cluster]]) <- want
	colnames(cors[[cluster]]) <- want
	for(gene1 in seq_along(use)){
		for(gene2 in seq_along(use)){
			temp <-  cor.test(
				expr_list[[cluster]][,which(colnames(expr_list[[cluster]]) == use[gene1])],
				expr_list[[cluster]][,which(colnames(expr_list[[cluster]]) == use[gene2])],
				method = "pearson",
				use = "complete.obs")
			ps[[cluster]][which(want == use[gene1]), which(want == use[gene2])] <- temp$p.value
			if(temp$p.value < 0.05)	{
				cors[[cluster]][which(want == use[gene1]), which(want == use[gene2])] <- temp$estimate
			}
		}
	}
}

sizeGrWindow(12,6)
par(mfrow =c(3,4))
par(mar = c(0.3, 0.2, 1.5, 0.2)) 

nSets <- length(clusters)
labels <- colnames(cors[[12]])
setNames=clusters
for (set in 1:nSets) {
	circlePlot(cors[[set]], order=seq(1:length(labels)), labels, main = setNames[set],
		variable.cex.labels = F,
		radii = c(0.56, 0.62),
		center = c(0.1, 0.04),
		min.cex.labels = 1.2,
		max.cex.labels = 1.4,
		cex.main = 1.4);
}




# SCRAP


my_fun <- function(cc, pp, tt){
	corrplot::corrplot(as.matrix(cc) , is.corr = T , type = "upper",diag = F, method = "circle",tl.cex = 1.0 ,cl.cex = 1.0, pch.cex = 10.0, title = tt,  mar=c(0,0,1,0), p.mat = pp, insig = "p-value", sig.level = -1, cl.pos = 'n') 
}









# Check enrichment results
res <- read.csv(paste(pathOut, "fgsea.trans_in_all_modules.GS.C11orf53.csv", sep = "/"), row.names=1)




