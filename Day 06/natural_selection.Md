# 1. Analyze SNP density along chromosome assemblies using sliding window analysis in VCFtools  
Modules: `bioinformatics/vcftools/0.1.14`  
Command: `vcftools --vcf <INPUT_VCF> --SNPdensity <WINDOW_SIZE>`  
## What happens when window size is changed?  

# 2. Plot distribution using ggplot2  
## Write R script
`library(ggplot2)`  
`library(ggplot2)`  
`siskin <- read.table("out.snpden", header = TRUE)`  
`pdf("siskin_density.pdf")`  
`qplot(siskin$BIN_START, siskin$VARIANTS.KB, data=siskin, geom= "point", colour=siskin$CHROM,main = "Siskin SNP Density", ylab="Density",xlab="ChromPosition") + labs(colour="Chromosome")`  
`dev.off()`  
Modules: `tools/R/3.4.1`  
Command: `Rscript <script_name.R>`  

# Download siskin_density.pdf using sftp/scp/etc

vcftools --vcf gVCF --LROH --chr <CHROMOSOME NAME> # Show LROHs per individual
vcftools --vcf gVCF --recode --out <INDIVIDUAL> --indv <INDIVIDUAL> # Make individual VCF
vcftools --vcf <INDIVIDUAl>.recode.vcf --recode --out <INDIVIDUAL>_hetsites --max-maf 0.9 --maf 0.1 # Get Heterozygous sites
vcftools --vcf <INDIVIDUAL>_hetsites.recode.vcf --SNPdensity 100000 --out <INDIVIDUAL>_hetsites
R

#In R
library(ggplot2)
siskin <- read.table("<INDIVIDUAL>_hetsites.snpden", header = TRUE)
pdf("<INDIVIDUAL>_heterozygosity.pdf")
qplot(siskin$BIN_START, siskin$VARIANTS.KB, data=siskin, geom= "point", colour=siskin$CHROM,main = "Siskin Heterozygosity", ylab="Heterozygosity",xlab="ChromPosition") + labs(colour="Chromosome")
dev.off()
q()
