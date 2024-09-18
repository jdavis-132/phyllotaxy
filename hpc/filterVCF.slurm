#!/usr/bin/bash
#SBATCH --job-name=filterVCF
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=25:00
#SBATCH --partition=jclarke,batch,guest
#SBATCH --output=/work/schnablelab/jdavis132/phyllotaxy2024/log/%x_%a.out
#SBATCH --error=/work/schnablelab/jdavis132/phyllotaxy2024/log/%x_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jdavis132@huskers.unl.edu

ml vcftools
ml bcftools

cd $WORK/phyllotaxy2024/genotype

#vcftools --gzvcf SAP.preprocess.vcf.gz --out SAP_phyllotaxyGenotypes --keep phyllotaxyGenotypes.tsv --recode
#vcftools --vcf SAP_phyllotaxyGenotypes.recode.vcf --out SAP_MAF_filter --maf 0.05 --recode
#bcftools filter SAP_MAF_filter.recode.vcf --exclude 'F_PASS(GT=="het") > 0.1' -o SAP_phyllotaxy_filtered.vcf
bcftools view -m2 -M2 -v snps SAP_phyllotaxy_filtered.vcf > SAP_phyllotaxy_SNPs.vcf
bcftools view --exclude-types indels SAP_phyllotaxy_SNPs.vcf > SAP_phyllotaxy_noSV.vcf
vcftools --vcf SAP_phyllotaxy_noSV.vcf --remove-indels --out SAP_phyllotaxy_biallelic.vcf --recode --recode-INFO-all

