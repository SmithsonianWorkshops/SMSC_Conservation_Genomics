#COMMANDS FOR ConGen Course
qrsh
module load bioinformatics/vcftools/0.1.14
vcftools --vcf gVCF --SNPdensity 100000 # Can play with this window size in course to show differences
module load tools/R/3.4.1
R

#In R
library(ggplot2)
siskin <- read.table("out.snpden", header = TRUE)
pdf("siskin_density.pdf")
qplot(siskin$BIN_START, siskin$VARIANTS.KB, data=siskin, geom= "point", colour=siskin$CHROM,main = "Siskin SNP Density", ylab="Density",xlab="ChromPosition") + labs(colour="Chromosome")
dev.off()
q()

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
