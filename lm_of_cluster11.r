# Bradley October 2021
# Linear modelling of the abundance of cluster 11 cells within each sample of the Smillie et al paper
# Agnostically identifying features (genes) associated with the abundance of this cluster
# Need to generate a model (will use elastic) to identify the significant associations

#~~~~ Overall question: What are the expression predictors of the abundance of the tuft cell resembling cluster? ~~~~~~

# Set up
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp")
source('analysis.r')
set="Epithelial"
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


# Load in the linear model packages - *Not already installed*
library(minfi)

# Pseudo-bulk the expression data: Want to determine which gene's expression on the whole are predictors of this
# Define the function
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

# Pseudo-bulk - without z-scoring so is more comparable
tmm <- tmm_PB(seur.integrated, 'RNA', z_score = F, by = "genes")

# Generate the meta, add the abundance of cluster 11  onto this
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
# Now add cluster 11
meta <- meta[match(rownames(freq), rownames(meta)),]
meta$cluster_11 <- freq$'11'
meta <- meta[,-which(colnames(meta)  == "Sample")]


#~~~~~~~~~~ 1. Do an independent linear regression between C11orf53 expression and the abundance of cluster 11 (our interesting cluster)
pathOut <- c("BH_analysis/Seurat/August_2021/Epithelial/lm_cluster11")
library("broom")

cluster11 <- meta$cluster_11
if(all(colnames(tmm) == rownames(meta)) == F){
	tmm <- tmm[,match(rownames(meta), colnames(tmm))]
}
tmm <- data.frame(t(tmm))

# Generate the model
lm_c11 <- lm(cluster11 ~ tmm$C11orf53)
lm_c11



write.csv(lm_c11, paste(pathOut, "tables/lm_C11orf53_vs_cluster11.csv", sep = "/"))


# Plot this
pdf(file=paste(pathOut, "plots/lm_C11orf53_vs_cluster11.pdf", sep = "/"))
par(mar=c(5,5,4,5))
plot(tmm$C11orf53, cluster11, xlab = "C11ORF53 expression", ylab = "Cluster 11 (%)", pch = 16, 
	cex.axis = 1.5, cex.lab = 2.0, cex.main = 2.0,
	main="Linear model of Cluster 11 \nand C11ORF53 expression")
abline(lm_c11, lty="dashed", col="red")
dev.off()


# Test the hypothesis
library("broom")
tidy(lm_c11)


#~~~~~~~~ 2. Now do this for every feature, and identify the top correlated features
# Generate the design matrix
library("limma")
design <- model.matrix(~cluster11)
dim(design)

# Fit a model on all features
tmm <- data.frame(t(tmm))
fit <- lmFit(tmm, design = design)

# perform the pooled estimation of standard errors that results in the moderated t-statistics and resulting p-values.
fit <- eBayes(fit)

# Extract the results
toptab <- topTable(fit, coef = 2, number = nrow(fit))
head(toptab, n=35)
write.csv(toptab, paste(pathOut, "tables/lm_all_tmm_vs_cluster11.csv", sep = "/"))


# Plot this
plot(toptab$logFC, -log10(toptab$adj.P.Val),
    xlab = "Effect size", ylab = bquote(-log[10](adj.p-value)),
    pch = 19
)
abline(h = -log10(0.05), lty = "dashed", col = "red")
sig_only <- toptab[toptab$adj.P.Val < 0.05,]
text(sig_only$logFC, -log10(sig_only$adj.P.Val), row.names(sig_only), cex=1, pos=3, col="red")


# Alternative plot
#sig_only <- toptab[toptab$adj.P.Val < 0.05,]
# Add diffexpressed column
toptab$diffexpressed <- ifelse(toptab$adj.P.Val < 0.05 & toptab$logFC > 1, "Yes", "No")
toptab$delabel <- NA
for(r in 1:nrow(toptab)){
	if(toptab$diffexpressed[r] == "Yes"){
		toptab$delabel[r] <- rownames(toptab)[r]
	}
}
# Add tanseQTLs
Peter_trans <- c("LRMP", "SH2D6", "PSTPIP2", "HTR3E", "TRPM5", "HTR3C", "ALOX5", "OGDHL", "BMX", "MATK", "SH2D7", "PIK3CG", "PLCG2", "PTGS1", "IL17RB", "AZGP1", "GNG13", "CAMP", "ANKHD1", "EIF4EBP", "GIN1", "SPAG6")
toptab$Peter <- ifelse(rownames(toptab) %in% Peter_trans, "Yes", "No")
toptab$Peter_label <- NA
toptab$Peter_diff <- NA
toptab$Peter_diff_label <- NA
for(r in 1:nrow(toptab)){
	if(toptab$Peter[r] == "Yes"){
		toptab$Peter_label[r] <- rownames(toptab)[r]
	}
	if(toptab$Peter[r] == "Yes" & toptab$diffexpressed[r] == "Yes"){
		toptab$Peter_diff[r] <- "Yes"
		toptab$Peter_diff_label[r] <- rownames(toptab)[r]
	}
}


# Plot
library(ggplot2)
library(ggrepel)
ggplot(data=toptab, aes(x=logFC, y=-log10(adj.P.Val), col= diffexpressed, label= Peter_diff_label)) + 
	geom_point() + 
	theme_bw() + 
	scale_color_manual(values=c("black", "red")) +
	geom_hline(yintercept=-log10(0.05), col="red", lty="dashed") + 
	geom_vline(xintercept=1, col="red", lty="dashed") + 
	geom_text_repel(size=5, max.overlaps=Inf, label.size=NA, fill=NA, seed=1234, nudge_x = .15, nudge_y = 0, box.padding = 0.7, xlim=c(1.2,NA)) + 
	ggtitle("Cluster 11 (%) associated genes") + 
	theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size = 20),  axis.title.y = element_text(size = 24), axis.title.x = element_text(size = 24), plot.title=element_text(face="bold", size = 28), legend.position = "none")
	
ggplot(data=toptab, aes(x=logFC, y=-log10(adj.P.Val), col= diffexpressed, label= Peter_diff_label)) + 
	geom_point() + 
	theme_bw() + 
	scale_color_manual(values=c("black", "red")) +
	geom_hline(yintercept=-log10(0.05), col="red", lty="dashed") + 
	geom_text_repel(size=5, label.size=NA, fill=NA, seed=134)
	



# Enrichment of our trans-eQTLs in nominally signifcantly associated genes
trans <- vector("list", length = 2)
names(trans) <- c("Vaughan_Shaw_HT12", "blue_module")
trans[[1]] <- c("LRMP", "SH2D6", "PSTPIP2", "HTR3E", "TRPM5", "HTR3C", "ALOX5", "OGDHL", "BMX", "MATK", "SH2D7", "PIK3CG", "PLCG2", "PTGS1", "IL17RB", "AZGP1", "GNG13", "CAMP", "ANKHD1", "EIF4EBP", "GIN1", "SPAG6")
trans[[2]] <- c("MYOM1", "MFSD7", "FANK1", "P2RY2", "KLK4", "LINC00261", "TRPM5", "PSTPIP2", "SH2D6", "ALOX5", "BMX", "GNG13", "SH2D7", "HCK", "PLCG2", "MATK")


# Do enrichment
library(gage,)
library(fgsea)
library(ggplot2)
library(ggridges)

# Using all genes
geneList <- toptab$logFC
names(geneList) <- rownames(toptab)
res <- fgseaMultilevel(trans, geneList, minSize = 7, maxSize = 500, scoreType = "pos", eps=0);
res <- as.data.frame(res[,-8])
Trans <- plotEnrichment(trans[[1]], geneList) + theme_bw() + ggtitle("Trans-eQTLs in linear model") + xlab("Rank") + ylab("Enrichment score") + theme(axis.text.x = element_text(size=19), axis.text.y = element_text(size = 20),  axis.title.y = element_text(size = 24), axis.title.x = element_text(size = 24), plot.title=element_text(face="bold", size = 27)) 
Blue <- plotEnrichment(trans[[2]], geneList) + theme_bw() + ggtitle("Blue hub genes in linear model") + xlab("Rank") + ylab("Enrichment score") + theme(axis.text.x = element_text(size=19), axis.text.y = element_text(size = 20),  axis.title.y = element_text(size = 24), axis.title.x = element_text(size = 24), plot.title=element_text(face="bold", size = 27)) 
write.csv(res, paste(pathOut, "tables/fgsea_res_trans_in_all_genes_lm_cluster11.csv", sep = "/"))

library(patchwork)
both <- ggarrange(Trans + rremove("ylab") + rremove("xlab"), Blue + rremove("xlab") + rremove("ylab"), ncol=1, nrow=2)
both <- annotate_figure(both, left = text_grob("Enrichment score", rot = 90, size = 22), bottom = text_grob("Rank", size=22))


# using sig only genes
geneList <- toptab[toptab$adj.P.Val < 0.05,]$logFC
names(geneList) <- rownames(toptab[toptab$adj.P.Val < 0.05,])







