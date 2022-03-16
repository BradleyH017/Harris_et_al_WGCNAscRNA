# Bradley August 2021
# Getting cluster proportions for DA testing in R

# Set up
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp")
source('analysis.r')
set <- "Epithelial"
pathOut = paste("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/BH_analysis/Seurat/August_2021", set, sep = "/")

# Load the object
load(paste(pathOut, "objects/Processed_seur.integrated.clusteredRds", sep = "/"))

# Generate the cluster counts matrix
freq <- data.frame(sample = factor(seur.integrated@meta.data$Sample), seurat_clusters = factor(seur.integrated@meta.data$seurat_clusters))
library(plyr)
tally <- ddply(freq, .(freq$sample, freq$seurat_cluster), nrow, .drop=F)
colnames(tally)<- c("orig.ident", "seurat_cluster", "freq")

#~~~ Reformat
samps <- levels(factor(tally$orig.ident))
indx <- seq(1:length(samps))
per_samp <- vector("list", length = length(samps))
names(per_samp) <- samps
for(x in indx){per_samp[[x]] = tally[tally$orig.ident == samps[x],c(2:3)];
rownames(per_samp[[x]]) <- per_samp[[x]][,1];
per_samp[[x]] <- per_samp[[x]][,-1];
colnames(per_samp[[x]]) <= levels(factor(tally$seurat_cluster))}
freq <- data.frame(matrix(unlist(per_samp), nrow = length(per_samp), byrow = T))
rownames(freq) <- names(per_samp)
colnames(freq) <- levels(factor(tally$seurat_cluster))

# Now convert to proportions
for(r in seq(1:nrow(freq))){
        freq[r,] <- (freq[r,]/rowSums(freq[r,]))*100
}


# Manually add the grouping (defined by the WGCNA across clusters)
colnames(freq) <- paste("cluster", colnames(freq), sep="_")
freq$turq_hub_cluster_grouping = c(rep("High", 4), rep("Low", 2), rep("High", 8), rep("Low", 2), rep("High", 5), "Low")

# Stats
library(rstatix)
library(ggpubr)
stat.test <- freq %>% wilcox_test(cluster_11 ~ turq_hub_cluster_grouping, paired = F)
stat.test$y.position = 4.2

bxp <- ggboxplot(freq, x = "turq_hub_cluster_grouping", y = "cluster_11", fill = "turq_hub_cluster_grouping", outlier.shape = NA, add = "jitter", palette = c("#DC0000FF", "#00A087FF")) +
        theme(legend.position="none") +
        scale_x_discrete(limits = c("High", "Low")) +
        xlab("turq_hub_cluster_grouping") + ylab("Cluster proportion (%)") +
        ggtitle("Cluster 11") + 
        theme(strip.background = element_blank()) +
        theme(axis.text.x = element_text(angle=45, hjust=1, size=14), axis.text.y = element_text(size = 14),  axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 0), plot.title = element_text(face = "bold", size = 18)) + 
        theme_classic()

bxp <- bxp + stat_pvalue_manual(stat.test, label = "p")


pdf(file = paste(pathOut, "plots/results/cluster11_abundance_across_turq_hub_grouping.pdf", sep = "/"))
bxp
dev.off()

# Cor between C11orf53 pseudo-bulk expression and cluster 11 abundance
# Pseudo-bulk the expression of each sample as in the WGCNA within clusters
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
        
        
 tmm <- tmm_PB(object = seur.integrated, assay = "RNA", z_score = T, by = "genes")
 
 # Extract C11orf53
tmm <- data.frame(t(tmm))
int <- tmm[, colnames(tmm) %in% c("C11orf53", "COLCA1", "COLCA2")]
names(C53) <- rownames(tmm)

# Add onto freq
int <- int[match(rownames(freq), rownames(int)),]
freq <- cbind(freq, int)

# Subset for what we want and melt
want <- c("cluster_11", "C11orf53", "COLCA1", "COLCA2")
sub <- freq[,colnames(freq) %in% want]
library(reshape)
sub$cluster_11 <- as.factor(sub$cluster_11)
melt <- melt(sub)
colnames(melt)[which(colnames(melt) == "variable")] <- "ciseQTL"

# Plot
ggscatter(freq, x = "C11orf53", y = "cluster_11",
          add = "reg.line", conf.int = TRUE,    
          add.params = list(fill = "lightgray"),
          )+
  stat_cor(method = "pearson", 
           label.x = -2, label.y = 3) 
           
melt$cluster_11 <- as.numeric(as.character(melt$cluster_11))
b <- ggplot(melt, aes(x = value, y = cluster_11)) + geom_point() +
  geom_smooth(method = "lm") + facet_wrap(~ciseQTL) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + stat_cor(method = "pearson", 
           label.x = -2, label.y = 3) + xlab("Pseudo-bulked C11orf53 expression") + theme_bw() + theme(strip.text.x = element_text(size = 12), axis.title.y = element_text(size = 14), axis.title.x = element_text(size=14))
	
pdf(file = paste(pathOut, "plots/results/C11orfs_cor_cluster11.pdf", sep = "/"))
b
dev.off()
	
	
	
	
	
	

