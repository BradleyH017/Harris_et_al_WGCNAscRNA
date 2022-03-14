# Bradley August 2021
# WGCNA within clusters
# Agnostic identification of modules correlated with C11orf53 expression within clusters.
# Requires the selection of sft threshold picks etc for each module
# Want to run this as an array job, doing each cluster independent of one another.

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
pathOut = paste("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/BH_analysis/Seurat/August_2021", set, sep = "/")

# Read in the option for the cluster interested in: Need to first generate the text file with each cluster as a row (0 - 11)
option_list = list(make_option(c("-c", "--cluster"), action = "store", default = NA, type ="character", help="cluster number"))
opt = parse_args(OptionParser(option_list=option_list))
c = opt$c;
print(paste("The option for cluster to study is", c, sep = " "))
		

# Define pseudo-bulking function
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


# Define options for the analysis
z_score = T
SCT = F
subset_variable_genes = T
by = "genes"
pseudo_bulk = T
nom_trans = F


# Load the seurat object
load(paste(pathOut, "objects/Processed_seur.integrated.clusteredRds", sep = "/"))


# Subset for cluster
int_seur <- subset(x = seur.integrated, idents = as.character(c))
DefaultAssay(int_seur) <- "log2TP10K"

# Find the variable features based on log2TP10K expression. Initially did the top 3000. Want to explore the top 5000 so matches bulkRNASeq
nfeatures = 5000
int_seur <-  FindVariableFeatures(int_seur, selection.method="vst", nfeatures = nfeatures)
var.features <- int_seur@assays$log2TP10K@var.features
int_genes <- c("ITPRID1", "POU2F3", "BMX", "SH2D7", "CHAT", "SH2D6", "HTR3E", "AZGP1", "ACTG1P22", "TRPM5", "OGDHL", "AVIL", "C11orf53", "COLCA1", "COLCA2")
print(paste("The intersect of turquoise module hub genes and variable features in cluster", c, "is")) 
intersect(var.features, int_genes)


# Extract the expression matrix. NOTE: Need to make sure am following the options for whether z_score is T or not. First attempt using cluster 11 showed enormouse modules, very little C11orf53 expression
if(pseudo_bulk == T){
	# Then execute the function
	exprdata <- tmm_PB(object = int_seur, assay = "RNA", z_score = z_score, by = by)
}



#Make sure to extract C11orf53 before subsetting
temp <- data.frame(t(exprdata))
C53 <- data.frame(Sample = rownames(temp), C11orf53 = temp$C11orf53)
exprdata <- exprdata[-which(rownames(exprdata) == "C11orf53"),]


if(subset_variable_genes == T){
        exprdata <- exprdata[rownames(exprdata) %in% var.features,]
}


if(pseudo_bulk == T){
        # Manually construct the meta
        meta <- data.frame(batch = rep("", length(levels(factor(seur.integrated@meta.data$Sample)))),
                Sex = c("F","F", "M","M", "F","F",  "M","M",  "M","M",  "M","M",  "F","F",  "F","F",  "F","F",  "M","M",  "F","F"),
                Site = c("R","R", "R","R", "R","R", "R","R", "R","R", "R","R", "R","R",  "","", "L","L", "R","R", "",""),
                Sample = levels(factor(seur.integrated@meta.data$Sample)))
        # Add turq_hub_cluster
	group_df <- data.frame(Sample = levels(factor(seur.integrated@meta.data$Sample)), turq_hub_cluster = c(rep("High", 4), "Low", "Low", rep("High", 8), "Low", "Low", rep("High", 5), "Low"))
	meta <- merge(meta, group_df, by = "Sample")
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
	# Binarise
	meta <- meta %>% mutate(batch=recode(batch, "train" = 0, "valid" = 1));
        meta <- meta %>% mutate(Sex=recode(Sex, "F" = 0, "M" = 1))
        meta <- meta %>% mutate(turq_hub_cluster=recode(turq_hub_cluster, "Low" = 0, "High" = 1))
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
	# Subset the meta for the samples present within this cluster
	meta <- meta[meta$Sample %in% colnames(exprdata),]
	# Add C11orf53 expression as a variable in the meta column
	meta <- merge(meta, C53, by = "Sample")
        rownames(meta) <- meta$Sample
        meta <- meta[,-which(colnames(meta) == "Sample")]	
}


# ~~~~~~~~~~~~~~~~~ Now can begin the analysis
# 1. QC
# Make sure these match
t_exprdata <- as.data.frame(t(exprdata))
meta <- meta[match(rownames(t_exprdata), rownames(meta)),]
y <- meta$C11orf53

pathOut <- paste(pathOut, "/WGCNA_within_clusters/cluster_", c, sep = "") 

pdf(paste(pathOut, "plots/QC/Sample_dendrogram_y.pdf", sep ="/"))
plotClusterTreeSamples(datExpr=t_exprdata, y=y)
dev.off()

# Check samples and genes
gsg=goodSamplesGenes(t_exprdata,verbose=3)
gsg$allOK


if (!gsg$allOK)
{
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(t_exprdata)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  t_exprdata = t_exprdata[gsg$goodSamples, gsg$goodGenes]
}

# Print the difference between interesting genes and those included. Many genes will be removed for having no expression
int_genes <- c("ITPRID1", "POU2F3", "BMX", "SH2D7", "CHAT", "SH2D6", "HTR3E", "AZGP1", "ACTG1P22", "TRPM5", "OGDHL", "AVIL", "C11orf53", "COLCA1", "COLCA2", "FCGBP", "MUC2", "CLCA1")
print("Interesting genes not included downstream =")
setdiff(int_genes, colnames(t_exprdata))

print(paste("The number of genes remaining in the cluster to analyse = ", ncol(t_exprdata)))


# Now do the dendro and traits
sampleTree2 = hclust(dist(t_exprdata), method = "average")
# Convert traits to a colour, white = low
traitColors = numbers2colors(meta, signed = FALSE);


#******** Soft threshold calculations 
powers = c(c(1:10), seq(from = 12, to=40, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(t_exprdata, powerVector = powers, verbose = 5)

cex1 = 0.8;

pdf(file = paste(pathOut, "plots/QC/Soft_threshold.pdf", sep="/"))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")
dev.off()

pdf(file = paste(pathOut, "plots/QC/Mean_connectivity.pdf", sep="/"))
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


# Save everything after removing the seurat objects
rm(seur.integrated)
rm(int_seur)
save.image(paste(pathOut, "objects/UpToSftThreshold.Rds", sep = "/"))



# As using numerous runs independently, want to pick sft threshold automatically. 
# Use the power estimate
softPower=sft$powerEstimate
adjacency = adjacency(t_exprdata, power = softPower);
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM;

#save(TOM, file = paste(pathOut, "/objects/TOM_sft", softPower, ".Rds", sep = ""))
#Just this once (09/06/21)
#save(TOM, file = "TOM_test.Rds")
save(TOM, file = paste(pathOut, "objects/TOM.Rds", sep = "/"))

# Now run the module calculation and merging steps

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
write.csv(moduleTraitCor, paste(pathOut, "tables/moduleTraitCor.csv", sep = "/"))
write.csv(moduleTraitPvalue, paste(pathOut, "tables/moduleTraitPvalue.csv", sep = "/"))


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
BH_pvals <- p.adjust(moduleTraitPvalue, method = "BH")
BH_pvals_df <- matrix(BH_pvals, ncol = ncol(moduleTraitPvalue), nrow = nrow(moduleTraitPvalue))
rownames(BH_pvals_df) <- rownames(moduleTraitPvalue)
colnames(BH_pvals_df) <- colnames(moduleTraitPvalue)

# plot the ones with any FDR correlations < 0.1
moduleTraitCorSig <- moduleTraitCor
for(r in 1:nrow(BH_pvals_df)){
	if(sum(BH_pvals_df[r,]<0.1) == 0){
		BH_pvals_df[r, ] <- rep(NA,ncol(BH_pvals_df))
	}
}
BH_pvals_df <- BH_pvals_df[complete.cases(BH_pvals_df),]
moduleTraitCorSig <- moduleTraitCorSig[rownames(moduleTraitCorSig) %in% rownames(BH_pvals_df),]
modnames <- gsub("ME", "", rownames(BH_pvals_df))
modcols <- rownames(BH_pvals_df)
names(meta)[grep("cluster", names(meta))] <- "Blue_hub_gene_group"


pdf(file = paste(pathOut, "plots/results/Module_trait_matrix_BH_pvals.pdf", sep = "/"))
## Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCorSig, 2), "\n(",
                   signif(BH_pvals_df, 2), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCorSig)
par(mar = c(8, 10, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCorSig,
               xLabels = names(meta),
               yLabels = modcols,
               ySymbols = modcols,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.0,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


# EXPERIMENTAL
remove_grouping <- T
if(remove_grouping == T & c == 11){
	# Remove blue hub gene cluster
	moduleTraitCor_no <- moduleTraitCor[,-which(colnames(moduleTraitCor) == "turq_hub_cluster")]
	moduleTraitPvalue_no <- moduleTraitPvalue[,-which(colnames(moduleTraitPvalue) == "turq_hub_cluster")]
	# p-value adjust
	BH_pvals_no <- p.adjust(moduleTraitPvalue_no, method = "BH")
	BH_pvals_df_no <- matrix(BH_pvals_no, ncol = ncol(moduleTraitPvalue_no), nrow = nrow(moduleTraitPvalue_no))
	rownames(BH_pvals_df_no) <- rownames(moduleTraitPvalue_no)
	colnames(BH_pvals_df_no) <- colnames(moduleTraitPvalue_no)
	# plot the ones with any FDR correlations < 0.1
	moduleTraitCorSig_no <- moduleTraitCor_no
	for(r in 1:nrow(BH_pvals_df_no)){
		if(sum(BH_pvals_df_no[r,]<0.1) == 0){
			BH_pvals_df_no[r, ] <- rep(NA,ncol(BH_pvals_df_no))
		}
	}
	BH_pvals_df_no <- BH_pvals_df_no[complete.cases(BH_pvals_df_no),]
	moduleTraitCorSig_no <- moduleTraitCorSig_no[rownames(moduleTraitCorSig_no) %in% rownames(BH_pvals_df_no),]
	modnames_no <- gsub("ME", "", rownames(BH_pvals_df_no))
	modcols_no <- rownames(BH_pvals_df_no)
	pdf(file = paste(pathOut, "plots/results/Module_trait_matrix_BH_pvals_no_grouping.pdf", sep = "/"))
	## Will display correlations and their p-values
	textMatrix_no = paste(signif(moduleTraitCorSig_no, 2), "\n(",
	                   signif(BH_pvals_df_no, 2), ")", sep = "");
	dim(textMatrix_no) = dim(moduleTraitCorSig_no)
	par(mar = c(8, 10, 3, 3));
	# Display the correlation values within a heatmap plot
	labeledHeatmap(Matrix = moduleTraitCorSig_no,
	               xLabels = names(meta)[-grep("cluster", names(meta))],
	               yLabels = modcols_no,
	               ySymbols = modcols_no,
	               colorLabels = FALSE,
	               colors = blueWhiteRed(50),
	               textMatrix = textMatrix_no,
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
C11orf53 <- meta$C11orf53

turq = as.data.frame(meta$turq_hub_cluster);
names(turq) = "turq_hub_cluster"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(t_exprdata, MEs, use = "p"));

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");

names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(t_exprdata, turq, use = "p"));

GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(turq), sep="");

names(GSPvalue) = paste("p.GS.", names(turq), sep="");

save.image(paste(pathOut, "objects/UpTo_geneTraitSignificance.Rds", sep = "/"))




#Can use the GS and MM to caluculate the weighting of genes and thier module membership in
#interesting modules
#module = "maroon" #Picked this based on the result of above heatmap. Change to look at different modules
#column = match(module, modNames);
#moduleGenes = moduleColors==module;


#pdf(file = paste(pathOut, "/plots/results/", module, "_scattergraph.pdf", sep = ""))
#verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
#                   abs(geneTraitSignificance[moduleGenes, 1]),
#                   xlab = paste("Module Membership in", module, "module"),
#                   ylab = "Gene significance for C11orf53 expression",
#                   main = paste("Module membership vs. gene significance\n"),
#                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
#dev.off()


geneInfo0 = data.frame(colnames(t_exprdata),
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, C11orf53, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.turq_hub_cluster));
geneInfo = geneInfo0[geneOrder, ]


write.csv(geneInfo, paste(pathOut, "tables/geneInfo.csv", sep = "/"))



# Scatter plot of genes MM with GS.C11orf53
if(c == 11){
	
	library(ggplot2)
	library(ggrepel)
	black <- geneInfo[geneInfo$moduleColor == "black",]
	Peter_trans <- c("LRMP", "SH2D6", "PSTPIP2", "HTR3E", "TRPM5", "HTR3C", "ALOX5", "OGDHL", "BMX", "MATK", "SH2D7", "PIK3CG", "PLCG2", "PTGS1", "IL17RB", "AZGP1", "GNG13", "CAMP", "ANKHD1", "EIF4EBP", "GIN1", "SPAG6")
	black$to_label <- ifelse(rownames(black) %in% Peter_trans, rownames(black), "")
	cor_res <- cor.test(black$MM.black, black$GS.C11orf53)
	cor_label <- paste("cor=", signif(cor_res$estimate,2), "\n", "  p=", signif(cor_res$p.value,2), sep = "")

	pdf(file = paste(pathOut, "plots/results/MM.black_GS.C11orf53.pdf",sep = "/"))
	ggplot(black, aes(x=MM.black, y=GS.C11orf53, label=to_label)) + 	
		geom_point(color = ifelse(black$to_label == "", "black", "red")) +
	#	geom_text(aes(label=to_label)) + 
		geom_text_repel(size=4, max.overlaps=Inf, label.size=NA, fill=NA, seed=123, box.padding = 1.35) + 
		xlim(0,1) + 
		ylim(0,1) + 
		geom_hline(yintercept=0.5, col = "blue", lty="dashed") +
		ggtitle("11q23.1 trans-eQTL targets in black module") +
		theme_bw() + 
		theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size = 12),  axis.title.y = element_text(size = 16), axis.title.x = element_text(size = 16), plot.title=element_text(face="bold", size = 16), legend.position = "none") + 
		annotate(geom="text", x=0.08, y=0.95, label=cor_label,
              color="black", size=5)
	dev.off()
}





