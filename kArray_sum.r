# Bradley August 2021
# scRNA-Seq - tidied pipeline, Summarising the results of kArray test
# kArray_sum.r

# Set up
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp")
source('analysis.r')
library('Seurat')
library(ggsci)
library(rstatix)



# Choose set and define path
set <- "Immune"
pathOut = paste("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/BH_analysis/Seurat/August_2021", set, "tables/markers", sep = "/")



# Read in all of the fgsea results as a list
setwd(pathOut)
temp = list.files(pattern="*.fgsea_vs_authors.csv")
temp.name <- temp
temp.name <- gsub(".fgsea_vs_authors.csv", "", temp.name)
temp.name <- gsub("k", "", temp.name)
res = lapply(temp, read.csv, row.names = 1)
names(res) <- temp.name


# Divide each element of the list into a list of each cluster
res.l <- vector("list", length = length(res))
names(res.l) <- temp.name
for(k in 1:length(res.l)){
	res.l[[k]] <- vector("list", length = length(levels(factor(res[[k]]$cluster))))
	for(cluster in seq(1:length(levels(factor(res[[k]]$cluster))))){
		res.l[[k]][[cluster]] <- res[[k]][res[[k]]$cluster == levels(factor(res[[k]]$cluster))[cluster],]
	}
	names(res.l[[k]]) <- levels(factor(res[[k]]$cluster))
}



#1. ~~~~~~~~~~~~ Summarising the number of confident hits: Where these is only one significant hit, or the top hit is 1000x more significant than the next, or the NES is 1.5* larger than that of the next
# Generate the matrix

res.conf <- res.l
for(k in 1:length(res.conf)){
	for(cluster in 1:length(res.l[[k]])){
		res.conf[[k]][[cluster]] <- res.conf[[k]][[cluster]][order(res.conf[[k]][[cluster]]$padj),]
		if(res.conf[[k]][[cluster]]$padj[2]/res.conf[[k]][[cluster]]$padj[1] > 1000){
			res.conf[[k]][[cluster]] <- res.conf[[k]][[cluster]][1,]
		} else {
			if(nrow(res.conf[[k]][[cluster]][res.conf[[k]][[cluster]]$padj < 0.05,]) == 1){
				res.conf[[k]][[cluster]] <- res.conf[[k]][[cluster]][1,]
			} else {
				res.conf[[k]][[cluster]] <- res.conf[[k]][[cluster]][order(-res.conf[[k]][[cluster]]$NES),]
				if(res.conf[[k]][[cluster]]$NES[1]/res.conf[[k]][[cluster]]$NES[2] > 1.5 & res.conf[[k]][[cluster]]$padj[1] < 0.05){
					res.conf[[k]][[cluster]] <- res.conf[[k]][[cluster]][1,]
				} else {
					if(nrow(res.conf[[k]][[cluster]]) > 1){
						res.conf[[k]][[cluster]] <- NA
							}
						}
					}
				}
			}
		res.conf[[k]] <- do.call(rbind, res.conf[[k]]);
		# Remove NA
		res.conf[[k]] <- res.conf[[k]][complete.cases(res.conf[[k]]),];
		# Remove those top of multiple clusters
		res.conf[[k]] <- res.conf[[k]][!duplicated(res.conf[[k]]$pathway),]
		res.conf[[k]]$k <- rep(names(res.conf)[k], nrow(res.conf[[k]]))
}


# Plotting the results
res.conf2 <- do.call(rbind, res.conf)
want <- c("cluster", "k")
res.conf2 <- res.conf2[,colnames(res.conf2) %in% want]
freq <- table(res.conf2$k)
df <- data.frame(k = names(freq), hits = data.frame(freq)$Freq)
df$k <- as.numeric(df$k)

library("ggsci")
df <- df[order(as.numeric(df$k)),]
plot(df$k, df$hits)
p<-ggplot(df, aes(x=k, y=hits, label=k)) + geom_line(size = 0.75) + geom_point() + geom_text(vjust = -0.5, nudge_x = 15) + theme_minimal() + ylab("Number of confident hits") + xlab("kNN paramater") + theme(axis.text=element_text(size=12),  axis.title=element_text(size=14,face="bold")) + theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

pdf(file = "../../kArray/Confident_hits_per_k.pdf")
p
dev.off()




#2. ~~~~~~~~~~~~ Summarise the number of times each of the authors clusters was a confident hit
# Tally the pathways for each kNN paramaeter. Also add those never found confidently
total <- do.call(rbind, res.l[[1]])
total <- levels(factor(total$pathway))


res.sum <- do.call(rbind, res.conf)
freq <- table(res.sum$pathway)
df <- data.frame(cluster = names(freq), hits = data.frame(freq)$Freq)
setdiff(total, df$cluster)
dfnot <- data.frame(cluster = total, hits = rep(0, length(total)))
df <- rbind(df, dfnot)
df <- df[order(-df$hits),]

p<-ggplot(df, aes(x=reorder(cluster, -hits), y=hits, fill=cluster)) + geom_bar(stat="identity") + theme_minimal() + theme(legend.position = "none") + theme(axis.text.x = element_text(angle=45, hjust=1, size=12), axis.text.y = element_text(size = 12), axis.title=element_text(size=14,face="bold")) + theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) + xlab("Authors cluster")+ ylab("Frequency found as confident hit")


pdf(file = "../../kArray/Confident_hits_per_authors_cluster.pdf")
p
dev.off()

