# Using sam tools to variant call from scRNAseq reads
# Recommended from https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1863-4 


#~~~~~~~ Ran this as script 
# Load the module
module load igmm/apps/bcftools/1.9
. /etc/profile.d/modules.sh
# Load and change dir
module load igmm/apps/bcftools/1.9
cd /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/Elmentaite/data

# Define
REFasta=ref/hg19/fasta/genome.fa
samples=samples.txt
samp=`head -n $SGE_TASK_ID $samples | tail -n 1 | awk '{ print $1 }'`
# Run
echo "Running on"
echo $samp
samtools mpileup -Q 30 -A -x -Ov -v -r 11  --fasta-ref $REFasta\
        "$samp"/outs/possorted_genome_bam.bam \
        > ../alignment/samtools_align/"$samp".chr11.vcf
echo "Done alignment"
# Change dir
cd ../alignment/samtools_align
# Process
bgzip "$samp".chr11.vcf

#~~~~~~~~~~~~~~~~~~~~~~~~~

# Combine and extract
cd /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/Elmentaite/alignment/samtools_align
for f in *.vcf.gz
do
tabix -p vcf $f
done

module load igmm/apps/bcftools/1.9
bcftools merge *.vcf.gz -Oz -o chr11.all.vcf.gz
tabix chr11.all.vcf.gz
tabix chr11.all.vcf.gz 11:111,156,830-111,156,840 > rs3087967.vcf






