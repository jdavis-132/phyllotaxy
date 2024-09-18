library(rMVP)
library(readr)
library(data.table)
library(dplyr)
library(tidyverse)

TOTAL_MARKERS <- 4693810 
EFFECTIVE_MARKERS <- 1088251.19
effective_ratio <- EFFECTIVE_MARKERS/TOTAL_MARKERS
## regular phenotype input

# load rmvp formatted data
MVP.Data(fileVCF="SAP_phyllotaxy_biallelic.vcf",
         filePhe="phenotypes.csv",
         sep.phe=",",
         fileKin=TRUE,
         filePC=TRUE,
         priority="memory",
         maxLine=10000,
         out="mvp.hmp"
)

genotype <- attach.big.matrix("mvp.hmp.geno.desc")
phenotype <- read.table("mvp.hmp.phe",head=TRUE)
map <- read.table("mvp.hmp.geno.map" , head = TRUE)
Kinship <- attach.big.matrix("mvp.hmp.kin.desc")
Covariates_PC <- bigmemory::as.matrix(attach.big.matrix("mvp.hmp.pc.desc"))

#FarmCPU bootstrapFarmCPU_signals
# args=commandArgs(TRUE)# receive argument x from terminal
# x=as.numeric(args) # x is the no. of bootstrap from 1 to 100

for(x in 1:100){
  phe1=phenotype # make a copy of phenotype
  nline = nrow(phe1)
  phe1[sample(c(1:nline), as.integer(nline*0.1)),2:ncol(phe1)]=NA  # randomly choose 10% phenotype to be NA 
  colnames(phe1)=paste0(colnames(phenotype),x)  # rename the phenotype by attaching bootstrap number
  for(i in 2:ncol(phe1)){
    imMVP<-MVP(phe = phe1[,c(1,i)], geno = genotype, map = map, K=Kinship, CV.FarmCPU=Covariates_PC,
               nPC.FarmCPU = 3, maxLoop = 10, method = "FarmCPU", priority = 'memory', p.threshold = (0.05/EFFECTIVE_MARKERS), 
               threshold = (0.05/effective_ratio),
               file.output = 'pmap.signal')
  }
}

# traits=c('am', 'am2', 'blue_glm', 'cmed3', 'med3', 'med2', 'cam3', 'cam3_pam', 'cam4_pam')


# get.support=function(trait){ # write a function to summarise the occurrence of signals, trait is what i have in the rmvp output filenames, disregarding the number of bootstrap
#   files = list.files(pattern = paste0(trait,".*FarmCPU_signals.csv"))
#   if (length(files)>=1){  
#     signals <-
#     files %>%
#     map_df(~read.csv(.,skip=1,header=F,colClasses = c("factor","integer","integer","factor","factor","numeric","numeric","numeric")))
#     header = read.csv(paste0(trait,"1.FarmCPU.csv"))
#     colnames(signals)=colnames(header)
#     colnames(signals)[8]<-"pvalue"
#     signals=signals %>%
#       group_by(SNP,CHROM,POS) %>%
#       summarise(P=mean(pvalue), RMIP = n()/100) #%>% ## if {{trait}} doesnot work otherwise change this name to something else 
#   #separate(SNP, c("CHR","BP"),remove=F)
#     write.table(signals, file=paste0("Z", trait, "signals.csv"), quote = F,row.names = F,sep=",")
#   }
#   else{
#       print(paste0("file not found", trait))
#     }
# }
# 
# for(x in traits){get.support(x)}

