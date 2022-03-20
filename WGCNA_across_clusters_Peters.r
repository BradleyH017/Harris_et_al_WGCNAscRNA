# Bradley June 2021
# WGCNA in pseudo-bulked scRNA-sequencing data
# WGCNA_scrna.r

# Set up
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp")
source('analysis.r')
library(edgeR)
library(fastcluster)
library(dynamicTreeCut)
library(WGCNA)
options(stringsAsFactors=F)


# Now load in the seurat object 
from_file = F
z_score = T
by = "genes"
nom_trans = T

# Set up
set <- "Epithelial"
pathOut = paste("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/BH_analysis/Seurat/August_2021", set, sep = "/")


# Define functions
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
}


# First time running this, want to generate the pseudo-bulked matrix
if(from_file ==T){
	# Load the seurat object (any will do, just want the raw counts)
	load(paste(pathOut, "objects/Processed_seur.integrated.clusteredRds", sep = "/"))
	#Pseudo-bulk
	tmm <- tmm_PB(seur.integrated, 'RNA', z_score = z_score, by = by);
	# Manually construct the meta
	meta <- data.frame(batch = rep("", length(levels(factor(seur.integrated@meta.data$Sample)))), 
		Sex = c("F","F", "M","M", "F","F",  "M","M",  "M","M",  "M","M",  "F","F",  "F","F",  "F","F",  "M","M",  "F","F"), 
		Site = c("R","R", "R","R", "R","R", "R","R", "R","R", "R","R", "R","R",  "","", "L","L", "R","R", "",""), 
		Sample = levels(factor(seur.integrated@meta.data$Sample)))
	meta$Subject <- unlist(strsplit(meta$Sample, "\\."))[c(T,F)]
        for(x in seq(1:nrow(meta))){
               if(as.numeric(gsub("N", "", meta$Subject[x])) < 40) {
                       meta$batch[x] <- "train"
               } else {
                       meta$batch[x] <- "valid"
               }
        }
	rownames(meta) <- NULL
	rownames(meta) <- meta$Sample;
	# Binarise / make variables in the meta continous
	meta <- meta %>% mutate(batch=recode(batch, "train" = 0, "valid" = 1));
	meta <- meta %>% mutate(Sex=recode(Sex, "F" = 0, "M" = 1))
	meta$Site.R <- rep(0, nrow(meta))
	for(r in seq(1:nrow(meta))){
		if(meta$Site[r] == "R"){
			meta$Site.R[r] <- 1;
		}
		if(meta$Site[r] == "T"){
			meta$Site.R[r] <- 1
		}
		if(meta$Site[r] == ""){
			meta$Site.R[r] <- NA
		}
		if(meta$Site[r] == "L"){
			meta$Site.R[r] <- 0
		}
	}
	meta <- meta[,-which(colnames(meta) == "Site")]
	meta <- meta[,-which(colnames(meta) == "Subject")]
	#Save after removing C11orf53, COLCA1, COLCA2 from the tmm (don't want to look for correlations against itself)
	C53 <- data.frame(t(tmm[rownames(tmm) %in% c("C11orf53", "COLCA1", "COLCA2"),]));
	C53$Sample <- rownames(C53)	
	meta <- merge(meta, C53, by = "Sample");
	rownames(meta) <- meta$Sample
	meta <- meta[,-which(colnames(meta) == "Sample")]
	tmm <- tmm[-c(which(rownames(tmm) == "C11orf53"),which(rownames(tmm) == "COLCA1"), which(rownames(tmm) == "COLCA2")),];
	write.csv(tmm, paste(pathOut, "WGCNA_across_clusters/tables/tmm_from_rawRNA.csv", sep = "/"))
	write.csv(meta, paste(pathOut, "WGCNA_across_clusters/tables/meta_binarized_inc_C11orf53.csv", sep = "/"))
	exprdata <- tmm
	} else {
		exprdata <- read.csv(paste(pathOut, "WGCNA_across_clusters/tables/tmm_from_rawRNA.csv", sep = "/"), row.names=1);
		meta <- read.csv(paste(pathOut, "WGCNA_across_clusters/tables/meta_binarized_inc_C11orf53.csv", sep = "/"), row.names = 1)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~ Begin the analysis
if(nom_trans == T){
	# Load the trans-eQTLs
	if(from_file == F){
		trans <- read.csv("Peter_nom_trans/Peter_nom_sig_trans.csv")
		# Manually correct some genes
		colnames(trans)[which(colnames(trans) == "gene")] <- "probe"
		colnames(trans)[which(colnames(trans) == "EnsemblGeneID")] <- "gene"
		colnames(trans)[which(colnames(trans) == "GeneName")] <- "human_SYMBOL"
		trans[which(trans$gene == "ENSG00000271615"), which(colnames(trans) == "human_SYMBOL")] <- "ACTG1P22"
                trans[which(trans$gene == "ENSG00000282776"), which(colnames(trans) == "human_SYMBOL")] <- "IGHV1-45"
                trans[which(trans$gene == "ENSG00000213250"), which(colnames(trans) == "human_SYMBOL")] <- "RMBS2P1"
                trans[which(trans$gene == "ENSG00000285665"), which(colnames(trans) == "human_SYMBOL")] <- "TOP1P1"
                trans[which(trans$gene == "ENSG00000262001"), which(colnames(trans) == "human_SYMBOL")] <- "DLGAP1-AS2"
                trans[which(trans$gene == "ENSG00000283393"), which(colnames(trans) == "human_SYMBOL")] <- "RPL23AP63"
                trans[which(trans$gene == "ENSG00000232142"), which(colnames(trans) == "human_SYMBOL")] <- "RPS25P9"
                trans[which(trans$gene == "ENSG00000228415"), which(colnames(trans) == "human_SYMBOL")] <- "PTMAP1"
                trans[which(trans$gene == "ENSG00000215481"), which(colnames(trans) == "human_SYMBOL")] <- "BCRP3"
                trans[which(trans$gene == "ENSG00000270069"), which(colnames(trans) == "human_SYMBOL")] <- "MIR222HG"
                trans[which(trans$gene == "ENSG00000171570"), which(colnames(trans) == "human_SYMBOL")] <- "RAB4B-EGLN2"
                trans[which(trans$gene == "ENSG00000236493"), which(colnames(trans) == "human_SYMBOL")] <- "EIF2S2P3"
                trans[which(trans$gene == "ENSG00000259172"), which(colnames(trans) == "human_SYMBOL")] <- "SNRPA1-DT"
                trans[which(trans$gene == "ENSG00000225031"), which(colnames(trans) == "human_SYMBOL")] <- "EIF4BP7"
                trans[which(trans$gene == "ENSG00000248476"), which(colnames(trans) == "human_SYMBOL")] <- "BACH1-IT1"
                trans[which(trans$gene == "ENSG00000257723"), which(colnames(trans) == "human_SYMBOL")] <- "CHCHD3P2"
                trans[which(trans$gene == "ENSG00000231993"), which(colnames(trans) == "human_SYMBOL")] <- "EP300-AS1"
                trans[which(trans$gene == "ENSG00000221930"), which(colnames(trans) == "human_SYMBOL")] <- "DENND10P1"
                trans[which(trans$gene == "ENSG00000235847"), which(colnames(trans) == "human_SYMBOL")] <- "LDHAP7"
                trans[which(trans$gene == "ENSG00000246100"), which(colnames(trans) == "human_SYMBOL")] <- "LINC00100"
		# Subset for the nominally significant - Leaves us with 2035 genes. 1635 are unique
		trans <- trans[trans$pvalue < 0.01,]
		# For the purpose of this analysis, cannot use the Illumina IDs for the clustering as the expression data is done using gens.
		# Therefore need to subset for the human_SYMBOLs
		trans <- trans[!is.na(trans$human_SYMBOL),]
		# Take the most significant hit for each duplicated gene
		trans <- trans[!duplicated(trans$human_SYMBOL),]
		# Leaves us with 1634 genes
		# write to file
		write.csv(trans, paste(pathOut, "WGCNA_across_clusters_Peter/tables/trans_nom_significant.csv", sep="/"))
	} else {
		trans <- read.csv(paste(pathOut, "WGCNA_across_clusters_Peter/tables/trans_nom_significant.csv", sep="/"), row.names = 1)
	};
	exprdata <- exprdata[rownames(exprdata) %in% trans$human_SYMBOL,]
	# Leaves us with 273 genes being tested (if using p<0.01)
}

# Save trans included
write.csv(trans[trans$human_SYMBOL %in% rownames(exprdata),], paste(pathOut, "WGCNA_across_clusters_Peter/tables/trans_nom_significant_in_smillie.csv", sep = "/"))
		
#1. QC
# Rearrange the meta so that it matches the order of samples in the expression matrix
t_exprdata <- as.data.frame(t(exprdata))
meta <- meta[match(rownames(t_exprdata), rownames(meta)),]
y <- meta$C11orf53


# Redefine the path
pathOut <- paste(pathOut, "WGCNA_across_clusters_Peter", sep = "/")

pdf(file = paste(pathOut, "plots/QC/SampleTreeForOutliers.pdf", sep="/"))
plotClusterTreeSamples(datExpr=t_exprdata, y=y)
dev.off()



# This searches for expression data with missing values. If no missing values, it returns TRUE. If not, then remove 
#the offending genes and samples 
gsg=goodSamplesGenes(t_exprdata,verbose=3)
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(t_exprdata)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(t_exprdata)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  t_exprdata = t_exprdata[gsg$goodSamples, gsg$goodGenes]
}

# All genes and samples are good


sampleTree2 = hclust(dist(t_exprdata), method = "average")
# Convert traits to a colour, white = low
traitColors = numbers2colors(meta, signed = FALSE);

pdf(file = paste(pathOut, "plots/QC/Dendro_and_traits_heatmap.pdf", sep="/"))
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(meta),
                    main = "Sample dendrogram and datTraits heatmap")
dev.off()




#******** Soft threshold calculations - To pick - Only run the first time because very time consuming
powers = c(c(1:10), seq(from = 12, to=40, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(t_exprdata, powerVector = powers, verbose = 5)

print(paste("The power estimate by the authors function =", sft$powerEstimate))
write.csv(sft$fitIndices, paste(pathOut, "/tables/sft_threshold_fit_indices.csv", sep = "/"))


# Plot the results:
cex1 = 0.8;


pdf(file = paste(pathOut, "plots/QC/Soft_threshold_pick.pdf", sep="/"))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")
dev.off()


# Mean connectivity
pdf(file = paste(pathOut, "plots/QC/Mean_connectivity.pdf", sep="/"))
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


# Save everything
save.image(paste(pathOut, "/UpToSftThreshold.Rds", sep = "/"))





# Calculate adjacency and turn into a topological overlap matrix
# Pick soft threshold and modify the path
softPower=sft$powerEstimate;


adjacency = adjacency(t_exprdata, power = softPower, type = "signed");
TOM = TOMsimilarity(adjacency, TOMType = "signed");
dissTOM = 1-TOM;
save(TOM, file = paste(pathOut, "TOM.Rds", sep = "/"))



# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
pdf(file = paste(pathOut, "plots/results/Gene-TOM-clustering.pdf", sep = "/"))
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut: (Note : Example uses deepsplit of 2, using 3 here based on previous analysis. See other WGCNA2.r script)
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)


# Plot the dendrogram and colors underneath
pdf(file = paste(pathOut, "plots/results/DynamicTreeCut.pdf", sep = "/"))
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()



# Merging similar modules

MEList = moduleEigengenes(t_exprdata, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");


pdf(file = paste(pathOut, "plots/results/Module_eigengene_clustering.pdf", sep = "/"))
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25      #This is usually set at 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()


# Call an automatic merging function
merge = mergeCloseModules(t_exprdata, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

# To see what this has done, plot the dynamic vs the merged
pdf(file = paste(pathOut, "plots/results/Dynamic_vs_merged.pdf", sep = "/"))
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#Rename the module Colours
# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;


####### Quantifying module-trait associations
nGenes = ncol(t_exprdata);
nSamples = nrow(t_exprdata);

# Recalculate MEs with color labels - NOTE - This looks for associations with datTraits, 
#so might want to go back to this and alter what you are looking for - 
#Also remember that these NEED to be in numeric form



MEs0 = moduleEigengenes(t_exprdata, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, meta, use = "p"); #Do this with both datTraits and datTraits reverse
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


# Then plot this - uncorrected p-values
pdf(file = paste(pathOut, "plots/results/Module_trait_matrix.pdf", sep = "/"))
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 10, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(meta),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.0,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


## Corrected p-values
correct_pvals <- T
if(correct_pvals == T){
	BH_pvals <- p.adjust(moduleTraitPvalue, method = "BH")
	BH_pvals_df <- matrix(BH_pvals, ncol = ncol(moduleTraitPvalue), nrow = nrow(moduleTraitPvalue))
	rownames(BH_pvals_df) <- rownames(moduleTraitPvalue)
	colnames(BH_pvals_df) <- colnames(moduleTraitPvalue)
	pdf(file = paste(pathOut, "plots/results/Module_trait_matrix_BH_pvals.pdf", sep = "/"))
	# Will display correlations and their p-values
	textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
	                   signif(BH_pvals_df, 1), ")", sep = "");
	dim(textMatrix) = dim(moduleTraitCor)
	par(mar = c(6, 10, 3, 3));
	# Display the correlation values within a heatmap plot
	labeledHeatmap(Matrix = moduleTraitCor,
	               xLabels = names(meta),
	               yLabels = names(MEs),
	               ySymbols = names(MEs),
	               colorLabels = FALSE,
	               colors = blueWhiteRed(50),
	               textMatrix = textMatrix,
	               setStdMargins = FALSE,
	               cex.text = 1.0,
	               zlim = c(-1,1),
	               main = paste("Module-trait relationships"))
	dev.off()
}







######Then want to quantify the associations between our actual genes and the trait of interest 
#(in this case the genotype). Also want to define a quantitative measure of module membership (MM)
#MM is the correlation of an eigen gene with with the gene expression profile it is based upon. 
#This therefore qantifies the similarity of all genes on the array to every module
# Define variable genotype containing the genotype column of datTrait
C11orf53 = as.data.frame(meta$C11orf53);
names(C11orf53) = "C11orf53"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(t_exprdata, MEs, use = "p"));

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");

names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(t_exprdata, C11orf53, use = "p"));

GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(C11orf53), sep="");

names(GSPvalue) = paste("p.GS.", names(C11orf53), sep="");

save.image(paste(pathOut, "UpTo_geneTraitSignificance.Rds", sep = "/"))


geneInfo0 = data.frame(colnames(t_exprdata),
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, C11orf53, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)){
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.C11orf53));
geneInfo = geneInfo0[geneOrder, ]


write.csv(geneInfo, paste(pathOut, "tables/geneInfo.csv", sep = "/"))

# Plot the correlation between MM and GS.C11orf53 for blue module
int <- geneInfo[geneInfo$moduleColor == "blue",]
pdf(file=paste(pathOut, "plots/results/Cor_blue_MM.blue_GS.C11orf53.pdf", sep = "/"))
par(mar=c(5,5,4,5))
plot(int$MM.blue, int$GS.C11orf53, main = "Correlation between GS.C11orf53 and MM.blue", xlab="MM.blue", ylab="GS.C11orf53", pch=16, col="blue", cex.axis = 1.5, cex.lab = 2.0, cex.main=2.0)
abline(lm(int$GS.C11orf53~int$MM.blue), col = "red", lty="dashed")
text(0.3,0.8, paste("cor=", signif(cor.test(int$GS.C11orf53, int$MM.blue, method="pearson")$estimate,2), "\n", "p=", signif(cor.test(int$GS.C11orf53, int$MM.blue, method="pearson")$p.value,2), sep = ""), cex=1.5)
dev.off()



# Plot the interesting module
library(GGally)
library(network)
library(sna)
library(ggplot2)
library(intergraph)
module="blue"
probes = names(t_exprdata)
inModule = is.finite(match(moduleColors, module));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
adj <- modTOM
adj[adj > 0.3] <- 1
adj[adj != 1] = 0
for(r in 1:nrow(adj)){
        adj[r,r] <- 0
}
dimnames(adj) <- list(modProbes, modProbes)

adj <- adj[,colSums(adj) > 0]
adj <- adj[rowSums(adj) > 0,]

# Build network
network = network(adj)
# Subset
network <- delete.vertices(network, colnames(adj[,colSums(adj) == 0]))
# Extract the names to make colours correct
genes_in <- sapply(network$val, function(x) {x$vertex.names} ) 
# Add FDR status colour
trans$FDR_sig <- ifelse(trans$FDR < 0.05, "tomato", "steelblue") %>% as.factor() 
to_add <- trans$FDR_sig
names(to_add) <- trans$human_SYMBOL
to_add <- to_add[names(to_add) %in% genes_in]
to_add <- to_add[match(genes_in, names(to_add))]
network %v% "rs3087967_trans-eQTL" = as.character(to_add)
#colpal = RColorBrewer::brewer.pal(9, "Set1")[ c(1,2) ]

pdf(file = paste(pathOut, "plots/results/Network_plot_turquoise.pdf", sep = "/"))
par(mar=c(2,2.5,2,2))
ggnet2(network, label=T,  color = "rs3087967_trans-eQTL",  
	mode = "kamadakawai", legend.position = "top")
dev.off()


# Identifying hub genes
# calculating the intramodular connectivity (the degree of co-expression of a given gene with respect to the genes of a particular module)
# Generate colours, a vector of length nGenes which the colour label for each gene
want <- c("colnames.t_exprdata.", "moduleColor")
col_df <- geneInfo[, colnames(geneInfo) %in% want]
t_exprdata <- t_exprdata[, order(colnames(t_exprdata))]
col_df <- col_df[order(col_df$colnames.t_exprdata.),]
modcols <- col_df$moduleColor


kIM <- intramodularConnectivity.fromExpr(t_exprdata, colors = modcols, power = 14)
rownames(kIM) <- colnames(t_exprdata)
rownames(col_df) <- col_df$colnames.t_exprdata.
kIM <- cbind(kIM, col_df)
int_IM <- kIM[kIM$moduleColor == module,]
int_IM <- int_IM[order(-int_IM$kWithin),]
#write.csv(int_IM, paste(pathOut, "/tables/", module, "_intramodular_connectivity_scores.csv", sep = ""))
write.csv(kIM, paste(pathOut, "tables/all_intramodular_connectivity_scores.csv", sep = "/"))



# Now do an enrichment of  Peter's trans-eQTLs against the genes of each module, ranked by their kWithin
trans <- vector("list", length = 1)
names(trans) <- c("Vaughan_Shaw_HT12")
trans[[1]] <- c("LRMP", "SH2D6", "PSTPIP2", "HTR3E", "TRPM5", "HTR3C", "ALOX5", "OGDHL", "BMX", "MATK", "SH2D7", "PIK3CG", "PLCG2", "PTGS1", "IL17RB", "AZGP1", "GNG13", "CAMP", "ANKHD1", "EIF4EBP", "GIN1", "SPAG6")


# Do enrichment
library(gage, lib.loc = "/exports/igmm/eddie/dunlop-lab/BradleyH/R/x86_64-pc-linux-gnu-library/4.0.2")
library(fgsea, lib.loc = "/exports/igmm/eddie/dunlop-lab/BradleyH/R/x86_64-pc-linux-gnu-library/4.0.2")
library(ggplot2, lib.loc = "/exports/igmm/eddie/dunlop-lab/BradleyH/R/x86_64-pc-linux-gnu-library/4.0.2")
library(ggridges, lib.loc = "/exports/igmm/eddie/dunlop-lab/BradleyH/R/x86_64-pc-linux-gnu-library/4.0.2")

modules = levels(factor(geneInfo$moduleColor))
fgseaResults = vector("list", length = length(modules))
names(fgseaResults) <- modules
for(x in seq(1:length(modules))){
	geneList <- kIM[kIM$moduleColor == modules[x],]$kWithin;
	names(geneList) <- rownames(kIM[kIM$moduleColor == modules[x],]);
	fgseaResults[[x]] <- fgseaMultilevel(trans, geneList, minSize = 7, maxSize = 500, scoreType = "pos", eps=0);
	fgseaResults[[x]] <- as.data.frame(fgseaResults[[x]]);
	print(paste("fgseaRes for the", modules[x], "module"))
	fgseaResults[[x]]
}
fgseaResults <- do.call(rbind, fgseaResults)
fgseaResults <- fgseaResults[,-8]

write.csv(fgseaResults, paste(pathOut, "tables/fgseaRes_eQTLs_in_all_modules_kIM.csv", sep = "/"))


pdf(file = paste(pathOut, "plots/results/fgsea.kIM.RNASeq_trans_in_blue.pdf", sep = "/"))
plotEnrichment(trans[['RNA_Seq']], geneList)
dev.off()

# Do the same but ranking by MM
modules = levels(factor(geneInfo$moduleColor))
fgseaResults = vector("list", length = length(modules))
names(fgseaResults) <- modules
for(x in seq(1:length(modules))){
	geneList <- geneInfo[geneInfo$moduleColor == modules[x],which(colnames(geneInfo) == paste("MM.", modules[x], sep = ""))];
	names(geneList) <- geneInfo[geneInfo$moduleColor == modules[x],]$colnames.t_exprdata.;
	fgseaResults[[x]] <- fgseaMultilevel(trans, geneList, minSize = 7, maxSize = 500, scoreType = "std", eps=0);
	fgseaResults[[x]] <- as.data.frame(fgseaResults[[x]]);
	print(paste("fgseaRes for the", modules[x], "module"))
}
fgseaResults <- do.call(rbind, fgseaResults)
fgseaResults <- fgseaResults[,-8]

write.csv(fgseaResults, paste(pathOut, "tables/fgseaRes_eQTLs_in_all_modules_MM.csv", sep = "/"))

pdf(file = paste(pathOut, "plots/results/fgsea.MM.RNASeq_trans_in_blue.pdf", sep = "/"))
plotEnrichment(trans[['Vaughan_Shaw_HT12']], geneList)
dev.off()

########### Dividing samples by their relative expression of interesting module hub genes

## Plotting a heatmap of samples with z-scores for the genes hub genes. Starting with those with adjacency > 0.3
hub_blue <- c("PLCG2", "HCK", "MATK", "TRPM5", "ALOX5", "HTR3E", "SH2D6", "PSTPIP2", "BMX", "PTGS1", "GNG13", "SH2D7")
hub_k <- rownames(kIM[kIM$moduleColor == "blue" & kIM$kWithin > 0.7,])
hub_MM <- geneInfo[geneInfo$moduleColor == "blue" & geneInfo$MM.blue > 0.7,]$colnames.t_exprdata.
hub <- intersect(hub_blue, intersect(hub_k, hub_MM))

test <- exprdata[rownames(exprdata) %in% hub,]
# z scoring across hub genes
z <- apply(test, 1, FUN=scale, center = T, scale = T)
rownames(z) <- colnames(test)


# order the genes by intramodular connectivity
z <- data.frame(t(z))
z <- z[order(rownames(z)),]
use <- kIM[kIM$moduleColor == "blue",]
use <- use[rownames(use) %in% rownames(z),]
use <- use[order(-use$kWithin),]
z <- z[match(rownames(z), rownames(use)),]
z <- data.frame(t(z))

#heatmap(z, scale = "none")
col<- colorRampPalette(c("blue", "grey", "red"))(256)
# Use RColorBrewer color palette names
z <- as.matrix(z)
#heatmap(z, scale = "none", col =  col, cexCol=0.8) 

library(dendextend)
dist_m <- dist(z)
dend <- hclust(dist_m, method="complete")
row_dend = as.dendrogram(dend)
#row_dend = color_branches(row_dend, k = 2) # `color_branches()` returns a dendrogram object

library(ComplexHeatmap)
pdf(file = paste(pathOut, "plots/results/Blue_hub_gene_heatmap_tidy.pdf", sep = "/"))
Heatmap(z, name = "z_score", column_title = "Blue module hub genes", row_title = "Samples", row_dend_reorder=T,
	row_title_gp = gpar(fontsize=17), column_title_gp = gpar(fontsize=17),
	row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize=17), cluster_rows = row_dend)
dev.off()

library(ComplexHeatmap)
pdf(file = paste(pathOut, "plots/results/Turquoise_hub_gene_heatmap.pdf", sep = "/"))
Heatmap(z, name = "z_score_across_genes", column_title = "Blue_hub_genes", row_title = "Samples", row_names_gp = gpar(fontsize = 10))
dev.off()


# how robust is this clustering method? - Bootstrap using 10 randomly selected genes n times
# This is based on euclidian hierarchical clustering of the rows into these 2 groups
# Redefine this into a vector
set.seed(123)
dist_m <- dist(z)
dend <- hclust(dist_m, method="complete")
clust_turq <- cutree(dend, k=2)

boot_list <- vector("list", length = 10000)

for(x in 1:10000){
	# Subset the matrix
	test <- exprdata[sample(1:nrow(exprdata), length(hub), replace = F),]
	# z scoring across genes
        z <- apply(test, 1, FUN=scale, center = T, scale = T)
	dist_m <- dist(z)
        dend <- hclust(dist_m, method="complete")
        boot_list[[x]] <- cutree(dend, k=2)
}

# Sum up the overlap from each iteration
inters <- sapply(boot_list, function(x){
	low = which(x == 2)
	correct_low <- length(intersect(low, which(clust_turq == 2)))
	correct_low / length( which(clust_turq == 2))
	}
)
	

# What is the proportion that I am getting correct randomly
all_cor <- length(inters[inters == 1])
all_cor / 10000


# I am only getting this correct 505 times in 10000. Which means that the likelihood I get this result p = 0.0505, making it interesting that I am getting it with the genes I have
	

