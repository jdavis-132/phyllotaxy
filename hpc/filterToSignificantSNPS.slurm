#!/usr/bin/bash
#SBATCH --job-name=filterToSignificantSNPS
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30G
#SBATCH --time=2:05:00
#SBATCH --partition=jclarke,schnablelab,batch,guest
#SBATCH --output=/work/schnablelab/jdavis132/phyllotaxy2024/log/%x_%a.out
#SBATCH --error=/work/schnablelab/jdavis132/phyllotaxy2024/log/%x_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jdavis132@huskers.unl.edu

ml vcftools

cd $WORK/phyllotaxy2024/genotype

vcftools --vcf SAP_phyllotaxy_biallelic.vcf.recode.vcf --snps sigSNPS.tsv --recode --out sigSNPSMed
