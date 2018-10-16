# 1. Analyze SNP density along chromosome assemblies using sliding window analysis in VCFtools  
Modules: `bioinformatics/vcftools/0.1.14`  
Command: `vcftools --vcf <INPUT_VCF> --SNPdensity <WINDOW_SIZE>`  
## What happens when window size is changed?  

# 2. Plot distribution using ggplot2  
Modules: `tools/R/3.4.1`  

## Write R script
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

# 5. Predict variant effects using SnpEff  

# Change default data directory in SnpEff.config file to /home/<YOURNAME>  

# Add siskin.genome to SnpEff.config file in the list of genomes.  

# Build siskin genome  
Modules: `java/1.8.0_45`  
Command: `mkdir genomes`  
Command: `mkdir siskin`  
Command: `mv <REFERENCE.fa> genomes/siskin.fa`  
Command: `mv <REFERENCE.gff> siskin/genes.gff`  
Command: `java -jar snpEff.jar build -gff3 -v siskin`  

# Compute variant effects using SnpEff  
Command: `java -jar snpEff.jar siskin <INPUT_VCF> > <OUT.vcf>`  

# Download snpEff_summary.html and inspect output  

# 6. Extract putatively selected genes from your VCF using SnpSift  
Modules: `java/1.8.0_45`  
Command: `java -jar SnpSift.jar filter "ANN[*].EFFECT has <missense_variant/stop_gained/etc>" <OUT.vcf> > <filtered.vcf>`  

SnpSift filter uses || and && for OR and AND respectively, allowing more stringent filtering.  
e.g. to get both missense and stop-gained sites, the filter line is: `"ANN[*].EFFECT has 'missense_variant' || ANN[*].EFFECT has 'stop_gained'"`  

# Compile list of possibly selected genes  
Command: `java -jar SnpSift.jar extractFields "ANN[*].GENE" <filtered.vcf> > genes.txt`  

This list will be non-unique, but can be manipulated in Excel or other tools to get a final gene list.  
