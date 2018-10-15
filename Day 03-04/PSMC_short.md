###1. PSMC
* Reminder to take notes, share info here: [Etherpad](https://pad.carpentries.org/CuuMC5spi7)
* For this part of the tutorial, we're going to work with a sorted bam file, from which we already removed the duplicates
	+ BAM file:  ```/scratch/genomics/tsuchiyam/varcall_psmc/siskin_all_nodup.bam```
	+ Here is the reference: ```/scratch/genomics/dikowr/siskin_SMSC/siskin_Contig3141_pilon.fasta```

*  First we will get stats and average coverage from our file:
	+ **module**: ```bioinformatics/samtools/1.6```
	+ **command**: ```samtools depth siskin_all_nodup.bam > final_depth.out```
```awk '{sum+=$3 ; n++} END {if(n>0) print sum/n;}' final_depth.out > averagecov.txt```
	+ Later you will need mindepth: 1/3 of the averagecov.txt and maxdepth: 2X averagecov.txt

* Now we will run a series of steps that will get us back to a fastq file.
	+ **module**: ```bioinformatics/samtools/1.6```
	+ **command**: ```samtools mpileup -C50 -uf <YOUR_SCAFFOLD.fasta> <nodup_YOUR_SORTED.bam> > <OUTPUT.bcf>```
	+ **module**: ```bioinformatics/bcftools```
	+ **command**: ```bcftools call -c -O v YOUR_FILE.bcf > OUTPUT.vcf ```
```vcfutils.pl vcf2fq -d <mindepth> -D <maxdepth> <YOUR_FILE.vcf> > OUTPUT_consensus.fq```

* Finally, we will run PSMC
	+ **module**: ```module load bioinformatics/psmc```
	+ **module**: ```module load gnuplot```
	+ **command**: ```fq2psmcfa -q20 <outconcensus>.fq > <outconsensus>.psmcfa```
```psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o output.psmc input.psmcfa```
```parallel -j$NSLOTS "psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" input.psmcfa -o round-{}.psmc" :::: <(seq 100)```
```psmc_plot.pl <output_prefix> input.psmc```
