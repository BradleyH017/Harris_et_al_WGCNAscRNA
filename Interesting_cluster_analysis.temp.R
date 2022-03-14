# Bradley August 2021
# scRNA-Seq - tidied pipeline, analysis of interesting clusters, gene variability etc
# Interesting_cluster_analysis.r

# Set up
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp")
source('analysis.r')
library('Seurat')
library(ggsci)
library(rstatix)
library(optparse)
library(patchwork)
library(gage)
library(fgsea)
library(data.table)
library(limma)
library(qusage)
library(farver)
library(labeling)

# Choose set and define path - Editting to use te unmerged clusters for epithelial (as not clusters were merged), but merged for Immune
set <- "Epithelial"
pathOut = paste("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/BH_analysis/Seurat/August_2021", set, sep = "/")

# Read in the object
load(paste(pathOut, "objects/Processed_seur.integrated.clusteredRds", sep = "/"))

#1. GSEA of the trans-eQTLs (Peter's and mine) against the markers of each cluster
trans <- read.csv("/exports/igmm/eddie/dunlop-lab/BradleyH/Human_RNA-Seq/data/2021PEER_FineMapping/N219_counts_nolog_TMM_INT_AGE_BMI_Sex_Site_Batch/eQTL/trans_eQTL.csv")
trans <- trans[,-1]
trans <- trans[trans$snp == "rs3087967",]

# Convert to gene names
library(AnnotationDbi)
library(org.Hs.eg.db)

MgeneList2entrez <- AnnotationDbi::mapIds(keys = as.character(trans$gene),  x = org.Hs.eg.db, keytype = "ENSEMBL", column = "SYMBOL",
                                          multiVals = "first");
DF_MgeneList2entrez <- data.frame("gene" = names(MgeneList2entrez), "human_SYMBOL" = MgeneList2entrez);
        rownames(DF_MgeneList2entrez) <- NULL;
        head(DF_MgeneList2entrez);
trans <- merge(trans, DF_MgeneList2entrez, by = "gene")
trans[trans$gene == "ENSG00000271615", which(colnames(trans) == "human_SYMBOL")] <- "ACTG1P22"
trans <- trans[!duplicated(trans$human_SYMBOL),]
trans <- trans[!is.na(trans$human_SYMBOL),]

trans_list <- vector("list", length = 3)
names(trans_list) <- c("nom_sig_RNASeq_rs3087967", "FDR_sig_RNASeq_rs3087967", "FDR_sig_Vaughan_Shaw_HT12")

trans_list[[1]] <- trans[trans$pvalue < 0.05,]$human_SYMBOL
trans_list[[2]] <- trans[trans$FDR < 0.1,]$human_SYMBOL
trans_list[[3]] <- c("LRMP", "SH2D6", "PSTPIP2", "HTR3E", "TRPM5", "HTR3C", "ALOX5", "OGDHL", "BMX", "MATK", "SH2D7", "PIK3CG", "PLCG2", "PTGS1", "IL17RB", "AZGP1", "GNG13", "CAMP", "ANKHD1", "EIF4EBP", "GIN1", "SPAG6")

# Now read in the markers of mine
if(set == "Epithelial"){
	markers_all <- read.csv(paste(pathOut, "tables/markers_unmerged_clusters/MAST_fully_processed_seur.csv", sep = "/"))
	}
if(set == "Immune"){
	markers_all <- read.csv(paste(pathOut, "tables/markers_merged_clusters/MAST_fully_processed_seur.csv", sep = "/"))
	}
markers_all <- markers_all[,-1]

# Split the markers into their clusters 
markers <- vector("list", length = length(levels(factor(markers_all$cluster))))
names(markers) <- levels(factor(markers_all$cluster))
for(x in seq(1:length(markers))){
        markers[[x]] <- markers_all[markers_all$cluster == names(markers)[x],]
}


# Now turn this into a ranked gene list (by logFC)
geneList <- vector("list", length = length(markers))
names(geneList) <- names(markers)
for(x in seq(1:length(markers))){
        geneList[[x]] <- markers[[x]]$avg_log2FC;
        names(geneList[[x]]) <- markers[[x]]$gene
}

# Now do enrichment for FDR and nom sig and manually BH correct
fgsea_list <- vector("list", length= length(geneList))
names(fgsea_list) <- names(geneList)
for(x in seq(1:length(fgsea_list))){
        fgsea_list[[x]] <- fgseaMultilevel(trans_list, geneList[[x]], minSize = 0, maxSize = 500, scoreType = "pos")
        fgsea_list[[x]]$cluster <- rep(names(geneList)[x], nrow(fgsea_list[[x]]))
}

fgseaRes <- do.call(rbind, fgsea_list)
fgseaRes

# Save these results
fgseaRes <- data.frame(fgseaRes)
fgseaRes <- fgseaRes[,-which(colnames(fgseaRes) == "leadingEdge")]
if(set == "Epithelial"){
	write.csv(fgseaRes, paste(pathOut, "tables/markers_unmerged_clusters/fgsea_trans_in_clusters.csv", sep ="/"))
	}
if(set == "Immune"){
	write.csv(fgseaRes, paste(pathOut, "tables/markers_merged_clusters/fgsea_trans_in_clusters.csv", sep ="/"))
	}


# Plot this 
pdf(file = paste(pathOut, "plots/results/Enrichment_of_my_FDR_trans_in_cluster_11.pdf", sep = "/")
plotEnrichment(trans_list[['FDR_sig_RNASeq_rs3087967']], geneList[['18']])
dev.off()





#2. Now analyse the variability in the expression of my FDR-significant trans-eQTLs across clusters and plot


# First need to define the pseudo-bulking function
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
                        if( by == "Samples"){
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


# Then perform the pseudo-bulking (and z-scoring) within clusters
if(set == "Epithelial"){
	split.by = "seurat_clusters"
}
if(set == "Immune"){
	split.by = "seurat_clusters_merged"
}

cluster_list <- SplitObject(seur.integrated, split.by = split.by)
for(x in seq(1:length(cluster_list))){ 
	cluster_list[[x]] <- tmm_PB(object = cluster_list[[x]], assay = "RNA", z_score = T, by = "Samples")
}


# Now plot the variability of the FDR significant trans-eQTLs, in addition to COLCA1/2. ACTGP122, ITPRID1 and CHAT not included in the analysis
genes.to.test <- c("C11orf53", "COLCA1", "COLCA2", "LRMP", "SH2D6", "PSTPIP2", "HTR3C", "ALOX5", "OGDHL", "MATK", "PIK3CG", "PLCG2", "PTGS1", "IL17RB", "GNG13", "CAMP", "ANKHD1", "GIN1", "SPAG6", "SH2D7", "BMX", "POU2F3", "HTR3E", "AZGP1", "TRPM5", "OGDHL", "AVIL")
genes.to.test <- intersect(genes.to.test, rownames(seur.integrated@assays$RNA))

rbind.named.fill <- function(x) {
    nam <- sapply(x, names)
    unam <- unique(unlist(nam))
    len <- sapply(x, length)
    out <- vector("list", length(len))
    for (i in seq_along(len)) {
        out[[i]] <- unname(x[[i]])[match(unam, nam[[i]])]
    }
    setNames(as.data.frame(do.call(rbind, out), stringsAsFactors=FALSE), unam)
}

int <- cluster_list
for(x in seq(1:length(int))){
	int[[x]] <- int[[x]][rownames(int[[x]]) %in% genes.to.test,]
}

# TRY 1
z_genes <- vector("list", length = length(genes.to.test))
names(z_genes) <- genes.to.test
for(y in seq(1:length(z_genes))){
        z_genes[[y]] <- vector("list", length = length(int));
        names(z_genes[[y]]) <- names(int);
        for(x in seq(1:length(z_genes[[y]]))){
                z_genes[[y]][[x]] <- int[[x]][which(rownames(int[[x]]) == genes.to.test[y]),]
        };
        z_genes[[y]] <-  as.data.frame(t(rbind.named.fill(z_genes[[y]])));
        colnames(z_genes[[y]]) <- names(int);
        z_genes[[y]] <- z_genes[[y]][,order(as.numeric(colnames(z_genes[[y]])))]
        z_genes[[y]]$Sample <- factor(rownames(z_genes[[y]]))
        z_genes[[y]] <- melt(z_genes[[y]])
        colnames(z_genes[[y]])[which(colnames(z_genes[[y]]) == "variable")] <- "seurat_cluster";
        colnames(z_genes[[y]])[which(colnames(z_genes[[y]]) == "value")] <- "zlogcpm";
}


library(ggplot2)
boxes <- vector("list", length = length(z_genes))
for(x in seq(1:length(z_genes))){
        boxes[[x]] <- ggboxplot(z_genes[[x]], "seurat_cluster", y = "zlogcpm") +
              geom_point(data=z_genes[[x]], aes(x=factor(seurat_cluster), y = zlogcpm, color = Sample), position = position_jitter(width = 0.0)) +
                        theme_classic() + ylab("z-logcpm expression") + xlab("Cluster") + ggtitle(names(z_genes)[x]) + theme(plot.title=element_text(face="bold", size = 22))
        }

library(ggpubr)
p2 <- ggarrange(boxes[[1]] + rremove("ylab") + rremove("xlab"), boxes[[2]] + rremove("ylab") + rremove("xlab"), boxes[[3]] + rremove("ylab") + rremove("xlab"),  boxes[[4]] + rremove("ylab") + rremove("xlab"),  boxes[[5]] + rremove("ylab") + rremove("xlab"),
        boxes[[6]] + rremove("ylab") + rremove("xlab"),  boxes[[7]] + rremove("ylab") + rremove("xlab"),  boxes[[8]] + rremove("ylab") + rremove("xlab"),
        boxes[[9]] + rremove("ylab") + rremove("xlab"),  boxes[[10]] + rremove("ylab") + rremove("xlab"), boxes[[11]] + rremove("ylab") + rremove("xlab"), boxes[[12]] + rremove("ylab") + rremove("xlab"),
         boxes[[13]] + rremove("ylab") + rremove("xlab"),  boxes[[14]] + rremove("ylab") + rremove("xlab"), boxes[[15]] + rremove("ylab") + rremove("xlab"), boxes[[16]] + rremove("ylab") + rremove("xlab"), 
          boxes[[17]] + rremove("ylab") + rremove("xlab"),  boxes[[18]] + rremove("ylab") + rremove("xlab"), boxes[[19]] + rremove("ylab") + rremove("xlab"), boxes[[20]] + rremove("ylab") + rremove("xlab"),
           boxes[[21]] + rremove("ylab") + rremove("xlab"),  boxes[[22]] + rremove("ylab") + rremove("xlab"), boxes[[23]] + rremove("ylab") + rremove("xlab"), boxes[[24]] + rremove("ylab") + rremove("xlab"), 
        boxes[[25]] + rremove("ylab") + rremove("xlab"),
        ncol = 5, nrow = 5, common.legend = TRUE, legend = "none")
p2 <- annotate_figure(p2, left = text_grob("z-logcpm expression", rot = 90, size = 24),
                    bottom = text_grob("Cluster", size = 24))

pdf(file = paste(pathOut, "plots/results/trans_eQTL_variability_by_cluster.pdf", sep = "/"))
p2
dev.off()


# 3. Do the statistical tests to investigate variability
tests <- c("max", "min", "range", "median", "mean", "sd", "mad")

summary.stats <- function(bulk, z_score, int_gene){
        if( z_score == T ){
                for(x in seq(1:length(bulk))){
                        genes <- rownames(bulk[[x]]);
                        bulk[[x]] <- apply(bulk[[x]], 2, FUN=scale, center = T, scale = T);
                        rownames(bulk[[x]]) <- genes
                }
        } else {};
        sum.table <- matrix(nrow = length(bulk), ncol = length(tests));
        rownames(sum.table) <- names(bulk);
        colnames(sum.table) <- tests;
                for(r in seq(1:length(bulk))){
                        sum.table[r, which(colnames(sum.table) == "max")] <- max(bulk[[r]][grep(int_gene, rownames(bulk[[r]])),]);
                        sum.table[r, which(colnames(sum.table) == "min")] <- min(bulk[[r]][grep(int_gene, rownames(bulk[[r]])),]);
                        sum.table[r, which(colnames(sum.table) == "range")] <- sum.table[r, which(colnames(sum.table) == "max")] - sum.table[r, which(colnames(sum.table) == "min")];
                        sum.table[r, which(colnames(sum.table) == "median")] <- median(bulk[[r]][grep(int_gene, rownames(bulk[[r]])),]);
                        sum.table[r, which(colnames(sum.table) == "mean")] <- mean(bulk[[r]][grep(int_gene, rownames(bulk[[r]])),]);
                        sum.table[r, which(colnames(sum.table) == "sd")] <- sd(bulk[[r]][grep(int_gene, rownames(bulk[[r]])),])
                        sum.table[r, which(colnames(sum.table) == "mad")] <- mad(bulk[[r]][grep(int_gene, rownames(bulk[[r]])),]);
                }
        sum.table <- as.data.frame(sum.table);
        sum.table <- sum.table[order(-sum.table$mad),]
        return(sum.table)
}

genes.to.test <- c(c("C11orf53", "COLCA1", "COLCA2"), c("LRMP", "SH2D6", "PSTPIP2", "HTR3E", "TRPM5", "HTR3C", "ALOX5", "OGDHL", "BMX", "MATK", "SH2D7", "PIK3CG", "PLCG2", "PTGS1", "IL17RB", "AZGP1", "GNG13", "CAMP", "ANKHD1", "EIF4EBP", "GIN1", "SPAG6" ))
z.sum.list = vector("list", length = length(genes.to.test))

# A. z-scoring
names(z.sum.list) <- genes.to.test
for(x in seq(1:length(z.sum.list))){
        z.sum.list[[x]] <- summary.stats(cluster_list, z_score = T, int_gene = genes.to.test[x]);
        write.csv(z.sum.list[[x]], paste(pathOut, "/tables/int_gene_variability_by_cluster/", genes.to.test[x], "_variability_analysis_zTMM.csv", sep = ""))
}


# B. Investigating the variability of gene expression in the absence of z-scoring --> Gives us a direct method of comparison with the other dataset
cluster_list_noz <- SplitObject(seur.integrated, split.by = split.by)
for(x in seq(1:length(cluster_list_noz))){ 
	cluster_list_noz[[x]] <- tmm_PB(object = cluster_list_noz[[x]], assay = "RNA", z_score = F, by = "orig.ident")
}

# now calculate statistics of variability
sum.list = vector("list", length = length(genes.to.test))
names(sum.list) <- genes.to.test
for(x in seq(1:length(sum.list))){
        sum.list[[x]] <- summary.stats(cluster_list_noz, z_score = F, int_gene = genes.to.test[x]);
        write.csv(sum.list[[x]], paste(pathOut, "/tables/int_gene_variability_by_cluster/", genes.to.test[x], "_variability_analysis_TMM.csv", sep = ""))
}


# C.  Looking at variability of gene expression at the single cell level, within clusters
# Divide into clusters
meta <- seur.integrated@meta.data
clusters <- levels(factor(meta$seurat_clusters))
cl_tog <- vector("list", length = length(clusters))
names(cl_tog) <- clusters
var_res <- vector("list", length = length(clusters))
names(var_res) <- clusters
Idents(seur.integrated) <- "seurat_clusters"
for(cluster in seq_along(clusters)){
	print(clusters[cluster])
	# Split
	cl_tog[[cluster]] <- subset(seur.integrated, idents = clusters[cluster])
	DefaultAssay(cl_tog[[cluster]]) <- "RNA"
	# Calculate gene variability
	cl_tog[[cluster]] <- FindVariableFeatures(cl_tog[[cluster]], selection.method="vst", nfeatures = 5000)
	# Grab results
	var_res[[cluster]] <- HVFInfo(object = cl_tog[[cluster]])
	# Subset for our genes
	var_res[[cluster]] <- var_res[[cluster]][rownames(var_res[[cluster]]) %in% genes.to.test, ]
	var_res[[cluster]] <- var_res[[cluster]][order(-var_res[[cluster]]$variance.standardized),]
	print(var_res[[cluster]])
	var_res[[cluster]]$cluster <- rep(clusters[cluster], nrow(var_res[[cluster]]))
}

# Collect results per gene and save
var_res_all <- do.call(rbind, var_res)
var_res_gene <- vector("list", length = length(genes.to.test))
names(var_res_gene) <- genes.to.test
for(gene in seq_along(genes.to.test)){
	var_res_gene[[gene]] <- var_res_all[grep(genes.to.test[gene], rownames(var_res_all)),] 
	var_res_gene[[gene]] <- var_res_gene[[gene]][order(-var_res_gene[[gene]]$variance),]
	write.csv(var_res_gene[[gene]], paste(pathOut, "/tables/int_gene_variability_by_cluster/", genes.to.test[gene], "_variability_analysis_seur_single_cells.csv", sep = ""))
}



# OLD# Stats for expression of the genes to test
#sum.sc = vector("list", length = length(genes.to.test))
#names(sum.sc) <- genes.to.test
#for(x in seq(1:length(sum.sc))){
#		print(genes.to.test[[x]])
#        sum.sc[[x]] <- summary.stats(cl_tog, z_score = F, int_gene = genes.to.test[x]);
#        sum.sc[[x]] <- sum.sc[[x]][order(-sum.sc[[x]]$sd),]
#        print(sum.sc[[x]])
#        write.csv(sum.sc[[x]], paste(pathOut, "/tables/int_gene_variability_by_cluster/", genes.to.test[x], "_variability_analysis_single_cells.csv", #sep = ""))
#}

### Use seurat to calculate gene variability at the single cell level, within clusters
#seur.sc <- vector("list", length = length(genes.to.test))
#names(sum.sc) <- genes.to.test



#4. Is the rank of samples based on C11orf53 expression limited to cluster 11, or shared across clusters?
# Extract the gene of interest and rank patients by the expression of this gene
# Define the function

library(purrr)
gene_rank <- function(x, gene, keep_gene_vals){
        # Generate the empty list for temp
        gene_rank_list <- vector("list", length = length(x));
        names(gene_rank_list) <- names(x);
        # Extract C11orf53 per cluster and rank
        for(y in seq(1:length(x))){
                gene_rank_list[[y]] <- x[[y]][grep(gene, rownames(x[[y]])),];
                gene_rank_list[[y]] <- data.frame(Sample = names(gene_rank_list[[y]]), ph = gene_rank_list[[y]]);
                colnames(gene_rank_list[[y]])[which(colnames(gene_rank_list[[y]]) == "ph")] <- gene;
                gene_rank_list[[y]] <- gene_rank_list[[y]][order(-gene_rank_list[[y]][,which(colnames(gene_rank_list[[y]]) == gene)]),];
                gene_rank_list[[y]]$rank <- seq(1:nrow(gene_rank_list[[y]]));
                colnames(gene_rank_list[[y]])[which(colnames(gene_rank_list[[y]]) == "rank")] <- paste("cluster.", names(x)[y], ".rank", sep = "")
        if(keep_gene_vals == T){
                colnames(gene_rank_list[[y]])[which(colnames(gene_rank_list[[y]]) == gene)] <- paste("cluster.", names(x)[y], ".", gene, sep = "")
                } else {
                gene_rank_list[[y]] <- gene_rank_list[[y]][,-grep(gene, colnames(gene_rank_list[[y]]))]
                }
        };
        # Combine into a single dataframe
        gene_ranks <- gene_rank_list %>% purrr::reduce(left_join, by = "Sample");
        rownames(gene_ranks) <- gene_ranks$Sample;
        gene_ranks <- gene_ranks[,-which(colnames(gene_ranks) == "Sample")];
        return(gene_ranks)
}

z <- cluster_list
for(x in seq(1:length(z))){
                        genes <- rownames(z[[x]]);
                        z[[x]] <- apply(z[[x]], 2, FUN=scale, center = T, scale = T);
                        rownames(z[[x]]) <- genes
                }

C53_gene_ranks <- gene_rank(z, "C11orf53", keep_gene_vals = F)


# Save the expression values
int_cluster = "11"

C53_gene_expr <- gene_rank(z, "C11orf53", keep_gene_vals = T)
C53_gene_expr <- C53_gene_expr[,-grep("rank", colnames(C53_gene_expr))]
colnames(C53_gene_expr) <- unlist(strsplit(colnames(C53_gene_expr), "\\."))[c(F,T,F)]
C53_gene_expr <- C53_gene_expr[,order(as.numeric(colnames(C53_gene_expr)))]
C53_gene_expr <- C53_gene_expr [order(-C53_gene_expr[which(colnames(C53_gene_expr) == int_cluster)]),]
write.csv(C53_gene_expr, paste(pathOut, "tables/markers_unmerged_clusters/C11orf53_zlogcpm_within_clusters.csv", sep = "/"))

# Save the gene ranks
colnames(C53_gene_ranks) <- unlist(strsplit(colnames(C53_gene_ranks), "\\."))[c(F,T,F)]
C53_gene_ranks <- C53_gene_ranks[,order(as.numeric(colnames(C53_gene_ranks)))]
C53_gene_ranks <- C53_gene_ranks [order(C53_gene_ranks[which(colnames(C53_gene_ranks) == int_cluster)]),]
write.csv(C53_gene_ranks, paste(pathOut, "tables/markers_unmerged_clusters/C11orf53_ranks_zlogcpm_within_clusters.csv", sep = "/"))


# Compute the pairwise spearman rank correlations between cluster ranks of patients
# First generate an empty matrix to fill with the values
mat <- matrix(nrow = ncol(C53_gene_ranks), ncol = ncol(C53_gene_ranks))

dimnames = colnames(C53_gene_ranks)
colnames(mat) <- dimnames
rownames(mat) <- dimnames


for(r in seq(1:nrow(mat))){
        for(c in seq(1:ncol(mat))){
                mat[r,c] <- cor(C53_gene_ranks[,colnames(mat)[c]],
                        C53_gene_ranks[,rownames(mat)[r]],
                        method = "spearman",
                        use = "complete.obs")
                }
        }

mat <- mat[order(as.numeric(rownames(mat))), order(as.numeric(colnames(mat)))]



# Now want to plot this, and save every thing
library(corrplot, lib='/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/temp_package_dir')

pdf(file = paste(pathOut, "plots/results/C11orf53_Spearmann_rank_correlations_within_clusters.pdf", sep = "/"))
corrplot::corrplot(as.matrix(mat) , is.corr = T , type = "upper",diag = F,method = "color",tl.cex = 1.0 ,cl.cex = 1)
dev.off()

write.csv(mat, paste(pathOut, "tables/markers_unmerged_clusters/C11orf53_Spearmann_rank_correlations_within_clusters.csv", sep = "/"))



# If want to incorporate the pvalue into this visualisation
mat2 <- matrix(nrow = ncol(C53_gene_ranks), ncol = ncol(C53_gene_ranks))
colnames(mat2) <- dimnames
rownames(mat2) <- dimnames

for(r in seq(1:nrow(mat2))){
        for(c in seq(1:ncol(mat2))){
                mat2[r,c] <- cor.test(C53_gene_ranks[,colnames(mat)[c]],
                        C53_gene_ranks[,rownames(mat)[r]],
                        method = "spearman",
                        use = "complete.obs")$p.value
                }
        }

mat2 <- mat2[order(as.numeric(rownames(mat2))), order(as.numeric(colnames(mat2)))]

# mat3 contains only the significant ones
mat3 <- mat
for(r in seq(1:nrow(mat2))){
        for(c in seq(1:ncol(mat2))){
                if(mat2[r,c] > 0.05){
                        mat3[r,c] <- 0
                } else {}
        }
}


pdf(file = paste(pathOut, "plots/results/C11orf53_Spearmann_rank_correlations_within_clusters_sig_only.pdf", sep = "/"))
corrplot::corrplot(as.matrix(mat3) , is.corr = T , type = "upper",diag = F,method = "color",tl.cex = 1.0 ,cl.cex = 1)
dev.off()


# Are the expression of cis-eQTLS (C11orf53, COLCA1 and COLCA2) correlated with one another? Both within and accross clusters
# 1.  Across all cells
# Re-pseudo-bulkd the expression, but without z-scoring
tmm_all <- tmm_PB(object = seur.integrated, assay = "RNA", z_score = T, by = "genes")

# Extract the cis-eQTLs
tmm_cis <- tmm_all[rownames(tmm_all) %in% c("C11orf53", "COLCA1", "COLCA2"), ]
tmm_cis <- data.frame(t(tmm_cis))

# Compute the pairwise correlations and pvals in expression
mat_cis <- matrix(nrow=ncol(tmm_cis), ncol = ncol(tmm_cis))
dims = colnames(tmm_cis)
rownames(mat_cis) <- dims
colnames(mat_cis) <- dims
for(r in seq(1:nrow(mat_cis))){
        for(c in seq(1:ncol(mat_cis))){
                mat_cis[r,c] <- cor(tmm_cis[,c],
                        tmm_cis[,r],
                        method = "pearson")
                }
        }

mat_cis_p <- matrix(nrow=ncol(tmm_cis), ncol = ncol(tmm_cis))
rownames(mat_cis_p) <- dims
colnames(mat_cis_p) <- dims
for(r in seq(1:nrow(mat_cis))){
        for(c in seq(1:ncol(mat_cis))){
                mat_cis_p[r,c] <- cor.test(tmm_cis[,c],
                        tmm_cis[,r],
                        method = "pearson")$p.value
                }
        }




# Save this as a plot
corrplot::corrplot(as.matrix(mat_cis) , is.corr = T , type = "upper",diag = F, method = "circle",tl.cex = 1.0 ,cl.cex = 1, title = "All_clusters",  mar=c(0,0,1,0), p.mat = mat_cis_p,  insig = "p-value", sig.level = -1) 

# 2. Now do this within clusters 
cluster_list <- SplitObject(seur.integrated, split.by = split.by)
for(x in seq(1:length(cluster_list))){ 
	cluster_list[[x]] <- tmm_PB(object = cluster_list[[x]], assay = "RNA", z_score = T, by = "genes")
}

cluster_list <- cluster_list[order(as.numeric(names(cluster_list)))]


cis_list <- vector("list", length = length(cluster_list))
cor_list <- vector("list", length = length(cluster_list))
p_list <- vector("list", length = length(cluster_list))
for(cluster in seq_along(cluster_list)){
	cis_list[[cluster]] <- cluster_list[[cluster]][rownames(cluster_list[[cluster]]) %in% c("C11orf53", "COLCA1", "COLCA2"),]
	cis_list[[cluster]] <- data.frame(t(cis_list[[cluster]]))
	cor_list[[cluster]] <- matrix(nrow=ncol(cis_list[[cluster]]), ncol = ncol(cis_list[[cluster]]))
	p_list[[cluster]] <- matrix(nrow=ncol(cis_list[[cluster]]), ncol = ncol(cis_list[[cluster]]))
	dims <- colnames(cis_list[[cluster]])
	rownames(cor_list[[cluster]]) <- dims
	colnames(cor_list[[cluster]]) <- dims
	rownames(p_list[[cluster]]) <- dims
	colnames(p_list[[cluster]]) <- dims
	for(r in seq(1:nrow(cor_list[[cluster]]))){
		for(c in seq(1:ncol(cor_list[[cluster]]))){
			cor_list[[cluster]][r,c] <- cor(cis_list[[cluster]][,c],
                        cis_list[[cluster]][,r],
                        method = "pearson")
            p_list[[cluster]][r,c] <- cor.test(cis_list[[cluster]][,c],
                        cis_list[[cluster]][,r],
                        method = "pearson", use = "complete.obs")$p.value
			}
		}
	}

titles <- names(cluster_list)


# Now plot together
my_fun <- function(cc, pp, tt){
	corrplot::corrplot(as.matrix(cc) , is.corr = T , type = "upper",diag = F, method = "circle",tl.cex = 1.0 ,cl.cex = 1.0, pch.cex = 10.0, title = tt,  mar=c(0,0,1,0), p.mat = pp, insig = "p-value", sig.level = -1, cl.pos = 'n') 
}


par(mfrow=c(4,3))
plot_list[[1]] <- my_fun(cc = cor_list[[1]], pp = p_list[[1]], tt = paste("Cluster", titles[1]))
plot_list[[1]] <- my_fun(cc = cor_list[[2]], pp = p_list[[2]], tt = paste("Cluster", titles[2]))
plot_list[[1]] <- my_fun(cc = cor_list[[3]], pp = p_list[[3]], tt = paste("Cluster", titles[3]))
plot_list[[1]] <- my_fun(cc = cor_list[[4]], pp = p_list[[4]], tt = paste("Cluster", titles[4]))
plot_list[[1]] <- my_fun(cc = cor_list[[5]], pp = p_list[[5]], tt = paste("Cluster", titles[5]))
plot_list[[1]] <- my_fun(cc = cor_list[[6]], pp = p_list[[6]], tt = paste("Cluster", titles[6]))
plot_list[[1]] <- my_fun(cc = cor_list[[7]], pp = p_list[[7]], tt = paste("Cluster", titles[7]))
plot_list[[1]] <- my_fun(cc = cor_list[[8]], pp = p_list[[8]], tt = paste("Cluster", titles[8]))
plot_list[[1]] <- my_fun(cc = cor_list[[9]], pp = p_list[[9]], tt = paste("Cluster", titles[9]))
plot_list[[1]] <- my_fun(cc = cor_list[[10]], pp = p_list[[10]], tt = paste("Cluster", titles[10]))
plot_list[[1]] <- my_fun(cc = cor_list[[11]], pp = p_list[[11]], tt = paste("Cluster", titles[11]))
plot_list[[1]] <- my_fun(cc = cor_list[[12]], pp = p_list[[12]], tt = paste("Cluster", titles[12]))





