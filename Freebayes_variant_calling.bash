# Bradley October 2021
# Variant calling using FreeBayes
# Following https://github.com/freebayes/freebayes

# This is to be ran after aligning the fast files to get counts using cell_ranger_counts.sh


# Load it with conda
#conda create --prefix /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/Elmentaite/alignment/alignment_conda freebayes
source activate /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/Elmentaite/alignment/alignment_conda

###### Need to change the bam file @RG bit for samples Human_colon_16S8002627 and Human_colon_16S8002630
# Make a temp copy, treat this, so don't have to run into issues later
cd /gpfs/igmmfs01/eddie/CCGG-tumour-WGS/BradTemp/Elmentaite/data/Human_colon_16S8002630/outs
cp possorted_genome_bam.bam test.bam
# View the RG bit of header
module load igmm/apps/samtools/1.6
samtools view -H test.bam | grep '^@RG'


###### test with picard - takes about 40 mins
cd /gpfs/igmmfs01/eddie/CCGG-tumour-WGS/BradTemp/Elmentaite/data/Human_colon_16S8002627outs
cp possorted_genome_bam.bam test.bam
java -jar /exports/igmm/eddie/CCGG-mouse/c11mouseWGS/WGS_download_buffer/picard.jar AddOrReplaceReadGroups \
       I=test.bam \
       O=testp.bam \
       RGID=Human_colon_16S8002627:0:1:1101:10311 \
       RGLB=0.1 \
       RGPL=ILLUMINA \
       RGPU=Human_colon_16S8002627:0:1:1101:10311 \
       RGSM=Human_colon_16S8002627
# Check - validation step takes about 10 mins
samtools view -H testp.bam | grep '^@RG'
samtools view -H testp.bam
samtools index testp.bam
java -jar /exports/igmm/eddie/CCGG-mouse/c11mouseWGS/WGS_download_buffer/picard.jar  ValidateSamFile I=testp.bam MODE=SUMMARY

# Now for the other sample
cd /gpfs/igmmfs01/eddie/CCGG-tumour-WGS/BradTemp/Elmentaite/data/Human_colon_16S8002630/outs
java -jar /exports/igmm/eddie/CCGG-mouse/c11mouseWGS/WGS_download_buffer/picard.jar AddOrReplaceReadGroups \
       I=test.bam \
       O=testp.bam \
       RGID=Human_colon_16S8002630:0:1:1101:10311 \
       RGLB=0.1 \
       RGPL=ILLUMINA \
       RGPU=Human_colon_16S8002630:0:1:1101:10311 \
       RGSM=Human_colon_16S8002630
# Check - validation step takes about 10 mins
samtools view -H testp.bam | grep '^@RG'
samtools view -H testp.bam
samtools index testp.bam
java -jar /exports/igmm/eddie/CCGG-mouse/c11mouseWGS/WGS_download_buffer/picard.jar  ValidateSamFile I=testp.bam MODE=SUMMARY
########



######### Alternative ----- Replacing ID and PU using sed
samtools view -H test.bam | sed "s/ID:Human_colon_16S8002630:0:1:1181:1031 1/ID:Human_colon_16S8002630:0:1:1181:10311/" | samtools reheader - test.bam > test2.bam;
samtools view -H test2.bam | sed "s/PU:Human_colon_16S8002630:0:1:1181:1031 1/PU:Human_colon_16S8002630:0:1:1181:10311/" | samtools reheader - test2.bam > test3.bam;
# Check
samtools view -H test3.bam
samtools view -H test3.bam | grep '^@RG'
#index (takes 5-10mins)
 samtools index test3.bam

# Now do the other sample
cd ../../Human_colon_16S8002627/outs/
cp possorted_genome_bam.bam test.bam
samtools view -H test.bam | sed "s/ID:Human_colon_16S8002627:0:1:2049:1031 1/ID:Human_colon_16S8002627:0:1:1101:10311/" | samtools reheader - test.bam > test2.bam;
samtools view -H test2.bam | sed "s/PU:Human_colon_16S8002627:0:1:2049:1031 1/PU:Human_colon_16S8002627:0:1:1101:10311/" | samtools reheader - test2.bam > test3.bam;
# Check
samtools view -H test3.bam
samtools view -H test3.bam | grep '^@RG'
# index
 samtools index test3.bam
#########

# Now can run freebayes on all samples - NOTE: For array 6 and 7 (the replaced header samples), need to run this on test3.bam
# Want to run this on all samples simultaneously
# Do this in an array script. Limited to just the region surrounding rs3087967 (11:111,156,830-111,156,840 ... SNP is 111,156,836)
freebayes -f $REFasta -r $region "$samp"/outs/testp.bam > temp/$samp.vcf



# bgzip all of the resulting vcfs and combine into one vcf
cd /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/Elmentaite/data/temp
module load igmm/apps/bcftools/1.9
for f in *.vcf 
do 
bgzip $f 
tabix -p vcf $f.gz
done

# Make sure all the zipped vcf files have been done so with bgzip!!!!
bcftools merge *.vcf.gz -Oz -o rs3087967.vcf.gz

gunzip rs3087967.vcf.gz

# Reading into R to save as a csv
module load igmm/apps/R/4.0.2
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/Elmentaite/data/temp")
lb="/exports/igmm/eddie/dunlop-lab/BradleyH/R/x86_64-pc-linux-gnu-library/4.0.2"
install.packages('vcfR', lib=lb)
library(vcfR, lib=lb)
vcf <- read.vcfR("rs3087967.vcf", verbose = FALSE )
write.csv(vcf@gt, "rs3087967.csv")














######### SCRAP
java -jar /exports/igmm/eddie/CCGG-mouse/c11mouseWGS/WGS_download_buffer/picard.jar AddOrReplaceReadGroups \
	I=possorted_genome_bam.bam \
	O=possorted_genome_bam.bam \
	RGID=Human_colon_16S8002630:0:1:1181:1031 \
	RGLB=0.1 \
	RGPL=ILLUMINA \
	RGPU=Human_colon_16S8002630:0:1:1181:1031 \
	RGSM=Human_colon_16S8002630




