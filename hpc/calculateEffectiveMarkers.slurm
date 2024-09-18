#!/usr/bin/bash
#SBATCH --job-name=calculateEffectiveMarkers
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=2:05:00
#SBATCH --partition=jclarke,batch,guest
#SBATCH --output=/work/schnablelab/jdavis132/phyllotaxy2024/log/%x_%a.out
#SBATCH --error=/work/schnablelab/jdavis132/phyllotaxy2024/log/%x_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jdavis132@huskers.unl.edu

ml plink2
ml anaconda
conda activate GWAS

cd $WORK/phyllotaxy2024/genotype

plink2 --vcf SAP_phyllotaxy_biallelic.vcf.recode.vcf --make-bed --out SAP_phyllotaxy_filtered

java -Xmx5g -jar ../gec.jar --effect-number --plink-binary SAP_phyllotaxy_filtered --genome --out SAP_phyllotaxy_gec
