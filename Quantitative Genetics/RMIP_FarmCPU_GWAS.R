library(rMVP)
library(data.table)
setwd("~/Documents/figuringout_GWAS/")
## regular phenotype input

# load rmvp formatted data
MVP.Data(fileVCF ="multispeqv3.vcf",
         filePhe="filtered_multispeq_cutoff.csv",
         #sep.hmp="\t",
         sep.phe=",",
         SNP.effect="Add",
         fileKin=TRUE,
         filePC=TRUE,
         priority="memory",
         #maxLine=10000,
         out="mvp.vcf"
)

genotype <- attach.big.matrix("mvp.vcf.geno.desc")
phenotype <- read.table("spatially_multispeq_blues_filtered_ordered_extremed_with_deviceid.csv", header = TRUE, sep = ",")
map <- read.table("mvp.vcf.geno.map" , head = TRUE)
Kinship <- attach.big.matrix("mvp.vcf.kin.desc")
Covariates_PC <- bigmemory::as.matrix(attach.big.matrix("mvp.vcf.pc.desc"))

#FarmCPU bootstrapFarmCPU_signals
# args=commandArgs(TRUE)# receive argument x from terminal
# x=as.numeric(args) # x is the no. of bootstrap from 1 to 100

for(x in 1:100){
  phe1=phenotype # make a copy of phenotype
  nline = nrow(phe1)
  phe1[sample(c(1:nline), as.integer(nline*0.1)), 2:ncol(phe1)]=NA  # randomly choose 10% phenotype to be NA
  colnames(phe1)=paste0(colnames(phenotype),x)  # rename the phenotype by attaching bootstrap number
  for(i in 2:ncol(phe1)){
    imMVP <- MVP(phe = phe1[,c(1,i)], geno = genotype, map = map, K=Kinship, CV.FarmCPU=Covariates_PC, file.output="pmap.signal",
                 nPC.FarmCPU = 3, maxLoop = 10, method = "FarmCPU", priority = 'memory',threshold=0.125, p.threshold=1.06e-8)
  }
}

traits=c('Relative_Chlorophyll', 'qL', 'PS1_Active.Centers', 'vH.', 'NPQt',
         'FvP_over_FmP', 'PhiNO', 'gH.', 'Phi2', 'PhiNPQ', 'ECS_tau',
         'PS1_Oxidized.Centers', 'PS1_Over.Reduced.Centers', 'PS1_Open.Centers')

get.support=function(trait){ # write a function to summarise the occurrence of signals, trait is what i have in the rmvp output filenames, disregarding the number of bootstrap
  files = list.files(pattern = paste0(trait,".*FarmCPU_signals.csv"))
  if (length(files)>=1){  
    signals <-
      files %>%
      map_df(~read.csv(.,skip=1,header=F,colClasses = c("factor","factor","integer","factor","factor","numeric","numeric","numeric")))
    header <- c("SNP","CHROM","POS","REF","ALT","Effect","SE","pvalue")
    colnames(signals)=header
    signals=signals %>%
      group_by(SNP,CHROM,POS) %>%
      summarise(P=mean(pvalue), support = n()/100) #%>% ## if {{trait}} doesnot work otherwise change this name to something else 
    #separate(SNP, c("CHR","BP"),remove=F)
    write.table(signals, file=paste0("Z", trait, "signals.csv"), quote = F,row.names = F,sep=",")
  }
  else{
    print(paste0("file not found", trait))
  }
}

for(x in traits){get.support(x)}
