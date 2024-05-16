suppressPackageStartupMessages(library("Rtsne"))
suppressPackageStartupMessages(require(VariantAnnotation))
suppressPackageStartupMessages(require(reshape))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(plotly))

#setwd("~/Desktop/HIPAA/super_auto/rapid_aut/")
#setwd("~/Desktop/HIPAA/super_auto/CLL_patients/251185/")
#setwd("~/Desktop/HIPAA/super_auto/")

args = commandArgs(trailingOnly=TRUE)
patient = args[1]
infile = args[2]
outfile = args[3]
vcf=readVcf(infile, "hg38")


#vcf=readVcf(paste("CLL_patients/",patient,"/",patient,".somatic.clustered.vcf",sep=""), "hg38")

print(patient)

#vcf=readVcf("251185.somatic.lichee_clustered.vcf", "hg38")
#vcf=readVcf("rapid_aut.subclones.vcf", "hg38")

AO = as.data.frame(geno(vcf)$AO); # Number of alt allele oberservations
RO = geno(vcf)$RO; # Number of ref allele overservations
DP = geno(vcf)$DP; # Read depth
SAF = as.vector(info(vcf)$SAF); # Alt obs on forward strand
SAR = as.vector(info(vcf)$SAR); # Alt obs on reverse strand
CLU = as.vector(info(vcf)$AFCLU);

selector=complete.cases(DP);

AF=matrix(rep(0, NROW(AO)*NCOL(AO)), ncol=NCOL(AO), byrow=T);
for(i in 1:NROW(AO)) {
  for(j in 1:NCOL(AO)) {
    AF[i, j] = AO[[i,j]] / (AO[[i,j]] + RO[[i, j]]);
  }
}

colnames(AF)=colnames(AO)
rownames(AF)=rownames(AO)

topo.colors(2)

strictSomatic=selector
print("Number of somatic variants")
print(sum(strictSomatic))
correctAF=AF

plotVariantCrossSamples=function(variant_vector=NA, samples,col,new_plot=TRUE){
  if (new_plot){
    plot(correctAF[strictSomatic,samples][1,],type='l', lwd=5, col='white',ylim=c(0,1),ylab='AF',xlab="r=0,grn=1,blu=2,blk=3,o=4,cy=5,prp=6,y=7,gry=8,pi=9,dg=10")
  }
  
  for(i in 1:sum(strictSomatic)){
    if (CLU[i] == "0"){
      lines(correctAF[strictSomatic,samples][i,], col="red")
    } else if (CLU[i] == "1"){
      lines(correctAF[strictSomatic,samples][i,], col="green")
    } else if (CLU[i] == "2"){
      lines(correctAF[strictSomatic,samples][i,], col="blue")
    } else if (CLU[i] == "3"){
      lines(correctAF[strictSomatic,samples][i,], col="black")
    } else if (CLU[i] == "4"){
      lines(correctAF[strictSomatic,samples][i,], col="orange")
    } else if (CLU[i] == "5"){
      lines(correctAF[strictSomatic,samples][i,], col="cyan1")
    } else if (CLU[i] == "6"){
      lines(correctAF[strictSomatic,samples][i,], col="purple")
    } else if (CLU[i] == "7"){
      lines(correctAF[strictSomatic,samples][i,], col="yellow")
    } else if (CLU[i] == "8"){
      lines(correctAF[strictSomatic,samples][i,], col="gray")
    } else if (CLU[i] == "9"){
      lines(correctAF[strictSomatic,samples][i,], col="pink")
    } else if (CLU[i] == "10"){
      lines(correctAF[strictSomatic,samples][i,], col="darkgreen")
    }
  }
}

pdf(file=outfile, width = 5, height = 5)
#pdf(file="251185.lichee.lines.pdf", width = 5, height = 5)
plotVariantCrossSamples(samples=c(1:length(AO)),col="Black")
dev.off()
