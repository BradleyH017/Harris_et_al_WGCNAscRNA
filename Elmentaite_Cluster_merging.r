# Bradley October 2021
# scRNASeq preprocessing of Elmentaite (Epithelial) cells, aligned using 10X 
# Testing the merging of clusters

# Set up
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp")
source('analysis.r')
library(optparse)
setwd("Elmentaite/results")

# Read in the object
load("objects/Processed_seur.integrated.clusteredRds")

# read in the results of the roc test for each comparison
res <- list.files("tables/pairwiseDEonSeuratClusters/roc")
res <- paste("tables/pairwiseDEonSeuratClusters/roc", res, sep = "/")
roc <- lapply(res, read.table, sep = '\t')
roc_raw <- do.call(rbind, roc)

# Using the non-like_smil approach ends up merging an enormous number of clusters, for this reason have left it the way it is
like_smil = T
if(like_smil == T){
        for(x in seq(1:length(roc))){
                roc[[x]]$gene <- rownames(roc[[x]])
                roc[[x]] <-  roc[[x]][roc[[x]]$power > 0.2,]
                }
        } else { for(x in seq(1:length(roc))){
                roc[[x]]$gene <- rownames(roc[[x]])
                roc[[x]] <-  roc[[x]][roc[[x]]$power > 0.45,];
                roc[[x]] <- roc[[x]][(abs(roc[[x]]$pct.1 - roc[[x]]$pct.2) > 0.3),];
        }
}




# Exploring the distribution of total DE genes based on this
counts <- vector("list", length = length(roc))
for(x in seq(1:length(roc))){
        counts[[x]] <- nrow(roc[[x]])
}
counts <- do.call(rbind, counts)

pdf(file = "tables/pairwiseDEonSeuratClusters/roc/hist_of_DE.pdf")
hist(counts, main = "Number of DE genes between clusters after filtering")
dev.off()

cuts <- data.frame(cutoff = seq(from = 0, to = 200, by = 10), nmerge = rep(0, length(seq(from = 0, to = 200, by = 10))))
for(r in seq(1:nrow(cuts))){
        cuts[r,2] <- length(counts[counts > cuts[r,1]])
}
write.csv(cuts, "tables/pairwiseDEonSeuratClusters/roc_cutoffs_after_filtering.csv")



# Generate a matrix, colour by the number of DE genes based on this subsetting
mat <- matrix(nrow = length(levels(factor(seur.integrated@meta.data$seurat_clusters))), ncol = length(levels(factor(seur.integrated@meta.data$seurat_clusters))))

dimnames = seq(1:(length(levels(factor(seur.integrated@meta.data$seurat_clusters)))))
dimnames = dimnames - 1

colnames(mat) <- dimnames
rownames(mat) <- dimnames


# Fill the matrix
roc_df <- do.call(rbind, roc)
roc_df$cond_pair <- rep(0,nrow(roc_df))
for(r in seq(1:nrow(roc_df))){
        if(as.numeric(roc_df$cond1[r]) > as.numeric(roc_df$cond2[r])){
                roc_df$cond_pair[r] = paste(roc_df$cond1[r], roc_df$cond2[r], sep = ".")
        } else {
                roc_df$cond_pair[r] = paste(roc_df$cond2[r], roc_df$cond1[r], sep = ".")
        }
}

for(c in seq(1:ncol(mat))){
        for(r in seq(1:nrow(mat))){
                if(as.numeric(colnames(mat)[c]) == as.numeric(rownames(mat)[r])){
                        mat[r,c] <- NA
                } else {
                        if(as.numeric(colnames(mat)[c]) > as.numeric(rownames(mat)[r])) {
                                mat[r,c] <- nrow(roc_df[roc_df$cond1 == colnames(mat)[c] & roc_df$cond2 == rownames(mat)[r],]);
                        } else {
                                mat[r,c] <- nrow(roc_df[roc_df$cond1 == colnames(mat)[r] & roc_df$cond2 == rownames(mat)[c],]);
                        }
                }
        }
}


# For visualisation, may want to reduce values so the max is 40
set_cap <- F

if(set_cap == T){
        minV <- 0
        maxV <- 40
        # Limit matrix
        mat <- apply(mat, c(1, 2), function(x) min(max(x,minV),maxV))
} else {}

library(corrplot, lib='/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/temp_package_dir')

pdf(file = "tables/pairwiseDEonSeuratClusters/pairwiseDE_after_filter_heatmap_maxed.pdf")
corrplot::corrplot(as.matrix(mat) , is.corr = FALSE , type = "upper",diag = F,method = "color",tl.cex = 1,cl.cex = 1)
dev.off()



# Now explore the number of clusters which would be merged for a range of values
cut_off <- 30
comp <- levels(factor(roc_df$cond_pair))
comp_test <- vector("list", length = length(comp))
for(x in seq(1:length(comp))){
        names(comp_test)[x] <- comp[x];
        comp_test[[x]] <- nrow(roc_df[roc_df$cond_pair == comp[x],])
}
comp_vec <- do.call(rbind, comp_test)
comp_vec <- comp_vec[,1]
names(comp_vec) <- names(comp_test)
comp_vec <- comp_vec[comp_vec < cut_off]

#••••••• Using the authors roc value threshold (power=0.2), but an increased cut off (30 DEGs), no clusters being merged - Have already do cluster marker analysis of these, so don't need to worry about recalculating










