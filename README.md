# Harris_et_al_WGCNAscRNA
Repository for code used in the WGCNA/scRNAseq analysis of colonic epithelium: https://www.nature.com/articles/s41598-022-17887-5#Sec21


## Pipeline order
To replicate the analysis performed in this paper, download the datasets into an accessible directory (see Methods), along with the analysis.r script (which will need to be adjusted for your own directories/R environment). Scripts can then be ran as follows

### Preprocessing and clustering - Smillie
  1. Pre_processing.r
  2. Initial_clustering_dimred.r
  3. kArray_clustering_gsea_vs_authors.r
  4. kArray_sum.r
  5. Doublet.detection.r
  6. Fully_processed_clustering.r
  7. DE_between_cluster.r
  8. Cluster_merging.r
  9. Interesting_cluster_analysis.temp.R

### WGCNA - Smillie
  1. WGCNA_across_clusters_Peters.r
  2. WGCNA_within_clusters.r
  3. Pairwise_module_preservation_WGCNA.r
  4. Pariwise_module_preservation_WGCNA_summary.r
  5. Gene-gene_correlations_intracultuster_from_WGCNA.r

### Abundance testing - Smillie
  1. lm_of_cluster11.r
  2. DA_milo.R

### Preprocessing and clustering - Elmentaite
  1. Elmentaite_alignment.txt
  2. Elmentaite_Preprocessing.r
  3. Elmentaite_Initial_clustering.r
  4. Elmentaite_kArray.r
  5. Elmentiate_kArray_sum.r
  6. Elmentaite_Doublet_removal.r
  7. Elmentaite_Fully_processed_clustering.r
  8. Elmentaite_DE_between_clusters.r
  9. Elmentaite_cluster_merging.r
  10. Elmentaite_Interesting_cluster_analysis.r

### Genotyping and gene variability analysis
  1. Freebayes_variant_calling.bash
  2. samtools_variant_calling.bash
  3. Gene_variability_summary_both_datasets.R

Any questions about the code or analysis please contact bradley.harris@ed.ac.uk
