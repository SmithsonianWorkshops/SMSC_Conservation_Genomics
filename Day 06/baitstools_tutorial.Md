# 1. Obtain BaitsTools from GitHub and check the help
## Obtain BaitsTools and make executable  
Module: `tools/ruby/2.2.2/`  
Command: `git clone https://github.com/campanam/BaitsTools`  
Command: `chmod +x BaitsTools/*.rb`  

## Check available subcommands  
Command: `ruby baitstools.rb`  

## Check subcommand-specific help  
Command: `ruby baitstools.rb <subcommand> -h`  

# 2. tilebaits and basic sequence filtration  
## Generate 120 bp baits with 60 bp tiling to produce a 2× coverage of the input FASTA/FASTQ sequences.  
Command: `ruby baitstools.rb tilebaits -i <INPUT.fa> -L 120 -O 60 -B`  

BaitsTools will generate a fasta file "out-baits.fa" that contains the unfiltered candidate bait sequences and a BED file "out-baits.bed" listing the bait sequence coordinates.  

## Use the `-c` argument to remove candidate baits shorter than the requested bait length. Remove gaps and Ns with `-G exclude` and `-N` respectively. Use the `-w` command to output a TSV format table of parameters.  
Command: `ruby baitstools.rb tilebaits -i <INPUT.fa> -L 120 -O 60 -B -w -c -G exclude -N`  

baitstools will generate a fasta file out-filtered-baits.fa" and a "out-filtered-baits.bed" containing the filtered bait sequences and coordinates respectively. It will also output a TSV-format file "out-filtered-params.txt" containing the results of bait filtration.  

## Open the out-filtered-params.txt in a spreadsheet program or text editor.  
The columns from left to right are:  
* Chromosome (Haplotype): Coordinates: The reference sequence, the haplotype identification (for alternate alleles or alignment data), and the bait coordinates.  
* BaitLength: The length of the generated bait  
* GC%: Bait GC content in percent  
* Tm: Bait melting temperature  
* Masked%: Percent of bait masked  
* MeanQuality: Mean Phred-like base quality of the generated bait  
* MinQuality: Minimum Phred-like base quality of the generated bait  
* Ns: Whether the bait included Ns  
* Gaps: Whether the bait included gaps  
* Kept: Whether the bait was retained in the final filtered set. 

## Filter the previously-made baits by target GC content and melting temperature.  
To use a default value, simply use the argument flag without a specified value. Use the default minimum (`-n`) and maximum (`-x`) GC contents of 30% and 50% respectively. Change the minimum (`-q`) and maximum (`-z`) melting temperatures to 80°C and 140°C.   
Command: `ruby baitstools.rb checkbaits -i out-filtered-baits.fa -w -n -x -q 80 -z 140`  

## Inspect the parameters file to see which baits were retained after filtration.  

# 3. annot2baits and baiting of specific genomic features  
annot2baits generates baits from a GFF/GTF annotation file and the corresponding regions in a FASTA or FASTQ file.  

Command: `ruby baitstools.rb annot2baits -i <INPUT.gff> -r <Input.fa> -U <DESIRED_FEATURES> -B`  

annot2baits will give the standard output and a FASTA file of the extracted regions (in this case "out-regions.fa").  

# 4. vcf2baits and baits from genomic variants  
vcf2baits selects variants from a vcf file and then generates baits from a reference sequence.  

## Specify the input VCF using `-i` and the reference sequence with `-r`.  
In addition to the common output and quality-control options, vcf2baits can filter VCF variants by their QUAL score (specified using `-V`) to help prevent the inclusion of sequencing errors in the bait set. The total number of desired variants is controlled with `-t` across all sequences in the reference file. The maximum number of variants per contig or scaffold is controlled with `-m`. The minimum distance distance between variants is controlled with `-d`. First, we will generate baits for 40 variants with QUAL score >= 30 and a minimum distance of 20,000 bp apart. Note that both `-t` and `-m` are set to 40 here since there is only one reference sequence in the file.   
Command: `ruby baitstools.rb vcf2baits -i <INPUT.vcf> -r <INPUT.fa> -V 30 -t 40 -m 40 -d 20000`  

This will produce a set of 120 bp baits with the selected variants at the 61st base position in each bait.  

## Bait length is controlled with `-L`. Use `-b` to specify the number of bases before the variant within the bait.    
Command: `ruby baitstools.rb vcf2baits -i <INPUT.vcf> -r <INPUT.fa> -V 30 -t 40 -m 40 -d 20000 -L 80 -b 20`  

## Baits can be tiled across the selected variants. Use the `-O` argument to specify the base pair offset between tiled baits and the `-k` argument to specify the number of baits per variant. Use the `-b` argument to determine the starting position for tiling. Ggenerate 80 bp baits from with 3× variant coverage and a 15 bp offset between baits and a starting postion 79 bp upstream of the variant, enter the following command:  
Command: `ruby baitstools.rb vcf2baits -i <INPUT.vcf> -r <INPUT.fa> -V 30 -t 40 -m 40 -d 20000 -L 80 -b 79 -O 15 -k 3`  

## Use `-p` to select variants without generating corresponding baits. No reference sequence is needed if `-p` is specified.   
Command: `ruby baitstools.rb vcf2baits -i <INPUT.vcf> -V 30 -t 40 -m 40 -d 20000 -p`

## Use `-e` to generate baits from every variant within a VCF file. Generate 120 bp baits from the variants selected in the previous step, enter the following command:  
Command: `ruby baitstools.rb vcf2baits -i out-selected.vcf -r <INPUT.fa> -e`  

## Use `-j` to scale the maximum number of selected variants per contig by individual contig length. This argument overrides `-m`.  
Command: `ruby baitstools.rb vcf2baits -i <INPUT.vcf> -r <INPUT.fa> -V 30 -j`  

vcf2baits will select up to 30,000 variants that are spaced minimally 10,000 bp apart.  

## Finally, use `-a` to apply alternate alleles to baits to generate a balanced bait set.   
Command: `ruby baitstools.rb vcf2baits -i <INPUT.vcf> -r <INPUT.fa> -V 30 -j -a`  

vcf2baits will select up to 30,000 variants that are spaced minimally 10,000 bp apart and apply alternate alleles to the bait sequences.  
