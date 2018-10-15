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

## Execute R script in job file  
Modules: `tools/R/3.4.1`  
Command: `Rscript <script_name.R>`  

## Download siskin_density.pdf using sftp/scp/etc  

## Do you observe any regions of low/high SNP density?  

# 3. Identify individual-specific ROH using VCFtools  
Command: `vcftools --vcf <INPUT_VCF> --LROH --chr <CHROMOSOME NAME> # Show LROHs per individual`

## Are any ROH identified in any individual? How long?  

# 4. Calculate and plot individual heterozygosity  
## Make individual VCF  

Command: `vcftools --vcf <INPUT_VCF> --recode --out <INDIVIDUAL> --indv <INDIVIDUAL>`  
## Get Heterozygous sites 
Command: `vcftools --vcf <INDIVIDUAl>.recode.vcf --recode --out <INDIVIDUAL>_hetsites --max-maf 0.9 --maf 0.1`  

## Get heterozygous densities  
Command: `vcftools --vcf <INDIVIDUAL>_hetsites.recode.vcf --SNPdensity <WINDOW_SIZE> --out <INDIVIDUAL>_hetsites`  

## Make R script and visualize density plot  

# 4. Identify individual-specific ROH using BCFtools  
Modules: `bioinformatics/bcftools/1.4.1`  
Command: `bcftools roh  -s <SAMPLE> -O r <INPUT_VCF.gz> > ROH.txt`  

BCFtools ROH requires a gzipped compressed VCF file.  Use `gzip` to compress your VCF.  

## Are the ROH identified the same as with VCFtools?  

