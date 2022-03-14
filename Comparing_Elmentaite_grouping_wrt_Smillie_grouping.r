# Bradley Jan 2022
# Comparing the expression of the blue module hub genes from Smillie et al WGCNA in Elmentaite et al

# Set up
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp")
source('analysis.r')
library(edgeR)
library(fastcluster)
library(dynamicTreeCut)
library(WGCNA)
options(stringsAsFactors=F)


# load the Smillie et al results for pan cluster WGCNA
from_file = F
z_score = T
by = "genes"
seurat_variable = F
nfeatures = 5000
MAD = F
nom_trans = T

# Load
set <- "Epithelial"
pathOut = paste("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/BH_analysis/Seurat/August_2021", set, sep = "/")

# Define function
tmm_PB <- function(object, assay, z_score, by, split.by){
        splitObj <- SplitObject(object, split.by = split.by);
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


# Load kIM and GI
pathOut <- paste(pathOut, "WGCNA_across_clusters_Peter", sep = "/")
kIM <- read.csv(paste(pathOut, "tables/all_intramodular_connectivity_scores.csv", sep = "/"), row.names=1)
geneInfo <- read.csv(paste(pathOut, "tables/geneInfo.csv", sep = "/"), row.names=1)

# Load the rest (need the expression we used)
load(paste(pathOut, "UpTo_geneTraitSignificance.Rds", sep = "/"))

# 1. ~~~~~ Regenerate the sample division results, just to check
hub_blue <- c("PLCG2", "HCK", "MATK", "TRPM5", "ALOX5", "HTR3E", "SH2D6", "PSTPIP2", "BMX", "PTGS1", "GNG13", "SH2D7")
hub_k <- rownames(kIM[kIM$moduleColor == "blue" & kIM$kWithin > 0.7,])
hub_MM <- geneInfo[geneInfo$moduleColor == "blue" & geneInfo$MM.blue > 0.7,]$colnames.t_exprdata.
hub <- intersect(hub_blue, intersect(hub_k, hub_MM))

# Extract the hub genes
hub_expr <- exprdata[rownames(exprdata) %in% hub,]
# z scoring across hub genes
z <- apply(hub_expr, 1, FUN=scale, center = T, scale = T)
rownames(z) <- colnames(test)
# order the genes by intramodular connectivity
z <- data.frame(t(z))
z <- z[order(rownames(z)),]
use <- kIM[kIM$moduleColor == "blue",]
use <- use[rownames(use) %in% rownames(z),]
use <- use[order(-use$kWithin),]
z <- z[match(rownames(z), rownames(use)),]
z <- data.frame(t(z))

col<- colorRampPalette(c("blue", "grey", "red"))(256)
# Use RColorBrewer color palette names
z <- as.matrix(z)
library(dendextend)
dist_m <- dist(z)
dend <- hclust(dist_m, method="complete")
row_dend = as.dendrogram(dend)

# Plot
library(ComplexHeatmap)
Heatmap(z, name = "z_score", column_title = "Blue module hub genes", row_title = "Samples", row_dend_reorder=T,
        row_title_gp = gpar(fontsize=17), column_title_gp = gpar(fontsize=17),
        row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize=17), cluster_rows = row_dend)

# Pseudo-bulking of the Smillie data was done by z-scoring across genes
# Need to load the raw Smillie data to compare with Elmentaite later
load(paste(pathOut, "../objects/Processed_seur.integrated.clusteredRds", sep = "/"))

# pseudo-bulk this WITHOUT z-scoring, so can compare with Elmentaite
tmm_noz <- tmm_PB(seur.integrated, 'RNA', z_score = F, by = "genes", split.by = "Sample");

# Extract the hub genes
tmm_noz_hub <- tmm_noz[rownames(tmm_noz) %in% hub,]

# 2. ~~~~~~~~ Now pseudo-bulk the Elmentaite data 
# Read in the object
load("Elmentaite/results/objects/Processed_seur.integrated.clusteredRds")

# Pseudo-bulk without z scoring
tmm_noz_elm <- tmm_PB(seur.integrated, 'RNA', z_score = F, by = "genes", split.by = "orig.ident");

# Extract the hub genes
tmm_noz_elm_hub <- tmm_noz_elm[rownames(tmm_noz_elm) %in% hub,]

# 3. ~~~~~~~~ Compare the expression of hub genes
# combine 
common_hub <- intersect(rownames(tmm_noz_elm_hub), rownames(tmm_noz_hub))
tmm_noz_elm_hub <- tmm_noz_elm_hub[rownames(tmm_noz_elm_hub) %in% common_hub,]
tmm_noz_hub <- tmm_noz_hub[rownames(tmm_noz_hub) %in% common_hub,]
tmm_noz_elm_hub <- tmm_noz_elm_hub[match(rownames(tmm_noz_hub), rownames(tmm_noz_elm_hub)),]
all_hub <- cbind(tmm_noz_elm_hub, tmm_noz_hub)

# Reformat to make boxplots of each hub gene across samples
all_hub <- as.data.frame(t(all_hub))



# 4. ~~~~~~~~ Alternative method - Using seurat to calculate the average expression of each hub gene, within samples (across all clusters)
# Load the Smillie data
load("BH_analysis/Seurat/August_2021/Epithelial/objects/Processed_seur.integrated.clusteredRds")
DefaultAssay(seur.integrated) <- 'RNA'

# Divide this into individual samples (all clusters) and calculate the average expression of the hub genes
splitObj <- SplitObject(seur.integrated, split.by ="Sample")
Sample <- levels(factor(seur.integrated@meta.data$Sample))
var_res <- vector("list", length = length(Sample))
names(var_res) <- Sample
for(samp in seq_along(Sample)){
	print(Sample[samp])
	splitObj[[samp]] <- FindVariableFeatures(splitObj[[samp]], selection.method="vst", nfeatures = 5000)
	# Grab results
	var_res[[samp]] <- HVFInfo(object = splitObj[[samp]])
	# Extract hub genes
	var_res[[samp]] <- var_res[[samp]][rownames(var_res[[samp]]) %in% hub,]
	var_res[[samp]]$Sample <- rep(Sample[samp], nrow(var_res[[samp]]))
	var_res[[samp]]$gene <- rownames(var_res[[samp]])
	print(head(var_res[[samp]]))
}

# Combine
var_res_all <- do.call(rbind, var_res)
rownames(var_res_all) <- NULL















