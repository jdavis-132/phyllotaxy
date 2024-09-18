#!/usr/bin/bash
#SBATCH --job-name=calculateLD
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30G
#SBATCH --time=2:05:00
#SBATCH --partition=jclarke,schnablelab,batch,guest
#SBATCH --output=/work/schnablelab/jdavis132/phyllotaxy2024/log/%x_%a.out
#SBATCH --error=/work/schnablelab/jdavis132/phyllotaxy2024/log/%x_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jdavis132@huskers.unl.edu

ml plink
ml anaconda
conda activate GWAS

cd $WORK/phyllotaxy2024/genotype

plink --bfile SAP_phyllotaxy_filtered --chr 5 --r2 --ld-snp Chr05_12109370 --ld-window-kb 99999 --ld-window 99999 --ld-window-r2 0 --out LD_Chr05_12109370
plink --bfile SAP_phyllotaxy_filtered --chr 5 --r2 --ld-snp Chr05_65733791 --ld-window-kb 99999 --ld-window 99999 --ld-window-r2 0 --out LD_Chr05_65733791
plink --bfile SAP_phyllotaxy_filtered --chr 6 --r2 --ld-snp Chr06_41390777 --ld-window-kb 99999 --ld-window 99999 --ld-window-r2 0 --out LD_Chr06_41390777
