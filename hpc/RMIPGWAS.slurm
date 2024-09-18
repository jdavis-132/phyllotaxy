#!/usr/bin/bash
#SBATCH --job-name=RMIPGWAS
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=24G
#SBATCH --time=55:00:00
#SBATCH --partition=jclarke,batch,guest
#SBATCH --output=/work/schnablelab/jdavis132/phyllotaxy2024/log/%x_%a.out
#SBATCH --error=/work/schnablelab/jdavis132/phyllotaxy2024/log/%x_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jdavis132@huskers.unl.edu

ml anaconda
conda activate GWAS

cd /scratch
cp $WORK/phyllotaxy2024/phyllotaxyPhenotypes.csv phenotypes.csv
cp $WORK/phyllotaxy2024/ResamplingFarmCPU_reseq.R RMIPGWAS.R
cp $WORK/phyllotaxy2024/genotype/SAP_phyllotaxy_biallelic.vcf.recode.vcf SAP_phyllotaxy_biallelic.vcf

R CMD BATCH RMIPGWAS.R

cp /scratch/*.csv $WORK/phyllotaxy2024/signals
cp /scratch/*.Rout $WORK/phyllotaxy2024/log
