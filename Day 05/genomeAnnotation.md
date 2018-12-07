#Genome Annotation - a mini guide

This is a short guide on how to annotate an eukaryotic genome using AUGUSTUS (Stanke et al. 2006). To begin with, you will need a genome assembly in FASTA format. 


##Summary of steps:
1. Filter out small scaffolds
2. Rename scaffolds sequentially
3. (optional) Create a name map
4. Run BUSCO using the --long flag
5. Run RepeatMasker and create an AUGUSTUS hints file (src=RM)
6. Run BLAT and create an AUGUSTUS hints file (src=E)
7. Edit the extrinsic file to add only the evidence you have
8. Run AUGUSTUS


##Filter out the short scaffolds and rename them using bioawk


###Step 1. Filter out the scaffolds below 500bp (or 1000bp, your choice :))

`bioawk -c fastx '{ if(length($seq) > 499) { print ">"$name; print $seq }}' input.fasta > filtered.fasta`

Observation: this bioawk command will remove any scaffolds of the specified size an below. Since I want to keep the 500bp scaffolds, I used 499 instead of 500.


###Step 2. Rename the scaffolds sequentially.

`bioawk -c fastx '{ print ">scaffold_" ++i"\n"$seq }' < filtered.fasta > ordered.fasta`

**Why are we doing this?**
For the AUGUSTUS run, we will create a *job array* to make our annotation more efficient. Job arrays consist on multiple tasks inside the same job, and those tasks have to be sequential (or not random at least). If we had just filtered out the small scaffolds and did not perform this step, we would end up with lots of of "empty" tasks. Like this:

- Imagine that we kept scaffolds 1 until 10, then 50 until 60. if we don't do the renaming step, the job array would look for the tasks 11 until 49 and would find nothing.

- If we had only 100 scaffolds, that would probably be a just a minor inconvenience. But since most projects will have thousands, (or dozens of thousands) scaffolds, that mess can become difficult to clean and filter out - and it also wastes computer resources.

###Step 3 (optional). Create a name map of the renaming step. 
In case you had already run other analyses, or you want to be able to revert to the original scaffold names, you can obtain the name map using this command:

`paste <(grep ">" filtered.fasta) <(grep ">" ordered.fasta) | sed 's/>//g' > name_mapping.txt`

The name map file will look like this:

```
scaffold_1	scaffold_1
scaffold_2	scaffold_2
scaffold_4	scaffold_3
scaffold_6	scaffold_4
scaffold_7	scaffold_5
scaffold_10	scaffold_6
scaffold_11	scaffold_7
scaffold_13	scaffold_8
scaffold_14	scaffold_9
scaffold_18	scaffold_10

```

First column shows the original scaffold name (from the filtered file) and the second the new scaffold id (after the renaming).


###Step 4. Run BUSCO using the --long flag.

BUSCO (Simão et al. 2015; Waterhouse et al. 2017) assesses completeness by searching the genome for a selected set of single copy orthologous genes. There are several databases that can be used with BUSCO and they can be downloaded from here: [https://buscos.ezlab.org]().

According to the BUSCO manual, the `--long` flag turns on Augustus optimization mode for self-training. It can be used as a training set for AUGUSTUS.

##### *** Before running BUSCO, copy the file augustus/config folder to a place where you have writing privileges (e.g. /pool/genomics/username/augustus) ***

####Job file

```
# /bin/sh
# ----------------Parameters---------------------- #
#$ -S /bin/sh
#$ -pe mthread 4
#$ -q mThM.q
#$ -l mres=16G,h_data=16G,h_vmem=16G,himem
#$ -cwd
#$ -j y
#$ -N job_name
#$ -o job_name.log
#
# ----------------Modules------------------------- #
module load bioinformatics/busco/3.0
#
# ----------------Your Commands------------------- #
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
echo + NSLOTS = $NSLOTS
#
export AUGUSTUS_CONFIG_PATH="/pool/genomics/user_id/augustus/config"
#
run_BUSCO.py --long -o output -i assembly.fa -l database_name -c $NSLOTS -m genome
#
echo = `date` job $JOB_NAME done
```

#####Explanation:
```
--long: turn on Augustus optimization mode for self-training.
-o: name of the output folder and files
-i: input file (FASTA)
-l: path to the folder containing the database of BUSCOs (the one you downloaded from the BUSCO website).
-c: number of CPUs
-m: mode (options are genome, transcriptome, proteins

```


#####4a. Using BUSCO output for the AUGUSTUS run
Copy the folder retraining\_parameters from run\_BUSCO/augustus\_output to your augustus/config/species. Rename the folder with the name of the run (you can find it by looking at the file prefix inside the folder). In this case, I'd rename the folder BUSCO\_siskin\_busco\_2684346740

```
BUSCO_siskin_busco_2684346740_exon_probs.pbl
BUSCO_siskin_busco_2684346740_igenic_probs.pbl
BUSCO_siskin_busco_2684346740_intron_probs.pbl
BUSCO_siskin_busco_2684346740_metapars.cfg
BUSCO_siskin_busco_2684346740_metapars.cgp.cfg
BUSCO_siskin_busco_2684346740_metapars.utr.cfg
BUSCO_siskin_busco_2684346740_parameters.cfg
BUSCO_siskin_busco_2684346740_weightmatrix.txt

```

You can also modify all filenames to match just the species ("siskin"). But in this case, you need to rename all files AND replace the current species id (BUSCO\_siskin\_busco\_2684346740) by siskin in all files using *sed*.


###Step 5. Masking and annotating repetitive elements
From the RepeatMasker manual (Smit, 2013-2015):
> Repeatmasker is a program that screens DNA sequences for interspersed repeats and low complexity DNA sequences.The output of the program is a detailed annotation of the repeats that are present in the query sequence as well as a modified version of the query sequence in which all the annotated repeats have been masked.

####Job file

```
# /bin/sh 
# ----------------Parameters---------------------- #
#$ -S /bin/sh
#$ -pe mthread 24-64
#$ -q mThC.q
#$ -l mres=4G,h_data=4G,h_vmem=4G
#$ -cwd
#$ -j y
#$ -N job_name
#$ -o job_name.log
#
# ----------------Modules------------------------- #
module load bioinformatics/repeatmasker
#
# ----------------Your Commands------------------- #
#
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
echo + NSLOTS = $NSLOTS
#
RepeatMasker -species carnivore -pa $NSLOTS -xsmall assembly.fa
#
echo = `date` job $JOB_NAME done

```
#####Explanation:
```
-species: species/taxonomic group repbase database (browse available species here: https://www.girinst.org/repbase/update/browse.php)
-pa: number of cpus
-xsmall: softmasking (instead of hardmasking with N)

```

#####Output files:
- assembly.fa.tbl: summary information about the repetitive elements
- assembly.fa.masked: hardmasked assembly
- assembly.fa.softmasked: softmasked assembly
- assembly.fa.out: detailed information about the repetitive elements, including coordinates, repeat type and size.


####Creating hints files from RepeatMasker
In this step, we will use the .out file from the RepeatMasker run to create a hint file for AUGUSTUS.

##### 5a. Use the script `rmOutToGFF3.pl` to convert your .out file into GFF3

`rmOutToGFF3.pl <input.out> > <output>.gff3`

##### 5b. Use the script `gff2hints.pl` convert the gff3 into a hints file. 
This script can be found here: [http://iubio.bio.indiana.edu/gmod/genogrid/scripts/gff2hints.pl]()

`perl gff2hints.pl --in=input.gff3 --source=RM --out=hints_RM.out`

###STEP 6. Running BLAT
BLAT (BLAST-like Alignment Tool, Kent 2002) is a tool that aligns DNA (as well as 6-frame translated DNA or proteins) to DNA, RNA and proteins across different species. The output of this program will also be used as hints for AUGUSTUS.

#####Job file
```
# /bin/sh
# ----------------Parameters---------------------- #
#$ -S /bin/sh
#$ -q mThC.q
#$ -l mres=6G,h_data=6G,h_vmem=6G
#$ -cwd
#$ -j y
#$ -N job_name
#$ -o job_name.log
#
# ----------------Modules------------------------- #
module load bioinformatics/blat
#
# ----------------Your Commands------------------- #
#
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
#
blat -t=dna -q=rna assembly.fa \
transcriptome.fna \
assembly_output.psl
#
echo = `date` job $JOB_NAME done
```

#####Explanation:
```
-t: target database type (DNA, 6-frame translated DNA or protein)
-q: query database type (DNA, RNA, protein, 6-frame translated DNA or RNA).

```
#### Creating hints files from BLAT
In this step, we will use the .psl file from the BLAT run to create a hint file for AUGUSTUS.

##### 6a. Sort the .psl file

`cat blat.psl | sort -n -k 16,16 | sort -s -k 14,14 > blat_srt.psl`

##### 6b. Use the script `blat2hints.pl` from the Augustus 3.3 module

`blat2hints.pl --in=blat_srt.psl --out=blat_hints.out`

###7. Edit the extrinsic file 
Check which columns of evidence you have in the extrinsic file (augustus/config/extrinsic) and remove the ones you don't need.

#####Example of extrinsic file with RepeatMasker (RM) and BLAT (E) evidence.

```
# extrinsic information configuration file for AUGUSTUS
# 
# include with --extrinsicCfgFile=filename
# date: 15.4.2015
# Mario Stanke (mario.stanke@uni-greifswald.de)


# source of extrinsic information:
# M manual anchor (required)
# P protein database hit
# E EST/cDNA database hit
# C combined est/protein database hit
# D Dialign
# R retroposed genes
# T transMapped refSeqs
# W wiggle track coverage info from RNA-Seq

[SOURCES]
M RM E

#
# individual_liability: Only unsatisfiable hints are disregarded. By default this flag is not set
# and the whole hint group is disregarded when one hint in it is unsatisfiable.
# 1group1gene: Try to predict a single gene that covers all hints of a given group. This is relevant for
# hint groups with gaps, e.g. when two ESTs, say 5' and 3', from the same clone align nearby.
#
[SOURCE-PARAMETERS]


#   feature        bonus         malus   gradelevelcolumns
#		r+/r-
#
# the gradelevel colums have the following format for each source
# sourcecharacter numscoreclasses boundary    ...  boundary    gradequot  ...  gradequot
# 

[GENERAL]
      start      1          1  M    1  1e+100  RM  1     1    E 1    1
       stop      1          1  M    1  1e+100  RM  1     1    E 1    1
        tss      1          1  M    1  1e+100  RM  1     1    E 1    1
        tts      1          1  M    1  1e+100  RM  1     1    E 1    1
        ass      1      1 0.1  M    1  1e+100  RM  1     1    E 1    1
        dss      1      1 0.1  M    1  1e+100  RM  1     1    E 1    1
   exonpart      1  .992 .985  M    1  1e+100  RM  1     1    E 1  1e2
       exon      1          1  M    1  1e+100  RM  1     1    E 1  1e4
 intronpart      1          1  M    1  1e+100  RM  1     1    E 1    1
     intron      1        .34  M    1  1e+100  RM  1     1    E 1  1e6
    CDSpart      1     1 .985  M    1  1e+100  RM  1     1    E 1    1
        CDS      1          1  M    1  1e+100  RM  1     1    E 1    1
    UTRpart      1     1 .985  M    1  1e+100  RM  1     1    E 1    1
        UTR      1          1  M    1  1e+100  RM  1     1    E 1    1
     irpart      1          1  M    1  1e+100  RM  1     1    E 1    1
nonexonpart      1          1  M    1  1e+100  RM  1     1.15 E 1    1
  genicpart      1          1  M    1  1e+100  RM  1     1    E 1    1

#
# Explanation: see original extrinsic.cfg file
#
```

The values in the columns 3 correspond respectively to bonus and malus (penalty) for each annotation element. Example from the AUGUSTUS manual:

```
Example: 
  CDS     1000  0.7  ....
means that, when AUGUSTUS is searching for the most likely gene structure,
every gene structure that has a CDS exactly as given in a hint gets
a bonus factor of 1000. Also, for every CDS that is not supported the
probability of the gene structure gets a malus of 0.7. Increase the bonus to
make AUGUSTUS obey more hints, decrease the malus to make AUGUSTUS predict few 
features that are not supported by hints.

```

Starting on column 4, it's possible to change how each source of evidence affects the annotation process. For example: 

```
Examples:

M 1 1e+100
means for the manual hint there is only one score class, the bonus for this
type of hint is multiplied by 10^100. This practically forces AUGUSTUS to obey
all manual hints.

T    2       1.5 1 5e29
For the transMap hints distinguish 2 classes. Those with a score below 1.5 and
with a score above 1.5. The bonus if the lower score hints is unchanged and
the bonus of the higher score hints is multiplied by 5x10^29.
```
Check here the original paper discussing the use of extrinsic evidence for annotation with AUGUSTUS (Stanke et al. 2006, [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1810548/]())

###8. Run AUGUSTUS

For the AUGUSTUS job, we need the following input files:

1. assembly (masked)
2. hints file (merged RM and E hints, see below)
3. extrinsic file
4. training parameters (from BUSCO)

#####To merge the hints files, use this script:
`cat hints_RM.out blat_hints.out | sort -n -k 4,4 | sort -n -k 5,5 > hints_RM_E.gff3`

####AUGUSTUS job:

```
# /bin/sh 
# ----------------Parameters---------------------- #
#$ -S /bin/sh
#$ -q mThM.q
#$ -l mres=20G,h_data=20G,h_vmem=20G,himem
#$ -cwd
#$ -j y
#$ -N augustus
#$ -o augustus.log
#
# ----------------Modules------------------------- #
module load bioinformatics/augustus/3.3
#
# ----------------Your Commands------------------- #
#
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
echo + NSLOTS = $NSLOTS
#
export AUGUSTUS_CONFIG_PATH="/pool/genomics/user_id/augustus/config"
#
augustus --strand=both --singlestrand=true \
--hintsfile=hints_RM.E.gff3 \
--extrinsicCfgFile=extrinsic.M.RM.E.cfg \
--alternatives-from-evidence=true \
--gff3=on \
--uniqueGeneId=true \
--softmasking=1 \
--species=training_parameters \
input_masked.fa > output.gff
#
echo = `date` job $JOB_NAME done
#
```
AUGUSTUS will run serially, one scaffold at a time. In order to speed up the process, we can break the assembly into scaffolds and process them in paralel. To do so, we will use a script from EVM to split the assembly and the hints file, and create job arrays for AUGUSTUS.

####8a. Partition the assembly into scaffolds. 
EVM (Evidence Modeller, Haas et al. 2008) is a program that combines *ab initio* gene predictions and protein and transcript alignments into weighted consensus gene structures. We will use an EVM script that splits the assembly into folders, with one scaffold per folder plus its corresponding hints file (in gff).

We don't have EVM installed as a module on Hydra, but you can download ([https://github.com/EVidenceModeler/EVidenceModeler/archive/v1.1.1.tar.gz]()) and extract it in your Hydra folder. The script is in the folder EVmutils. This script runs fast, so we can use the interactive queue for it.

```
perl /PATH/TO/EVM/EvmUtils/partition_EVM_inputs.pl --genome genome.fasta \
     --gene_predictions hints_RM_E.gff3 \
     --segmentSize 100000 --overlapSize 10000 \
     --partition_listing partitions_list.out
```

Now you should have many (thousands) of folders, each one with one scaffold and its corresponding hints file. They all retained the same name or the original file, and the folders are identified as scaffold\_1, scaffold\_2... scaffold\_n.

Now, let's create a job array to run the augustus job. A job array is a job file with multiple tasks. We can create a job array with 50 tasks or 50 separate job files. There are multiple ways of creating a job array, and this is just an example. 

```
# /bin/sh 
# ----------------Parameters---------------------- #
#$ -S /bin/sh
#$ -q mThM.q
#$ -l mres=20G,h_data=20G,h_vmem=20G,himem
#$ -t 1-100 -tc 20
#$ -cwd
#$ -j y
#$ -N augustus
#$ -o augustus_$TASK_ID.log
#
# ----------------Modules------------------------- #
module load bioinformatics/augustus/3.3
#
# ----------------Your Commands------------------- #
#
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
echo + NSLOTS = $NSLOTS
#
export AUGUSTUS_CONFIG_PATH="/pool/genomics/user_id/augustus/config"
#
augustus --strand=both --singlestrand=true \
--hintsfile=scaffold_$SGE_TASK_ID/hints_RM.E.gff3 \
--extrinsicCfgFile=extrinsic.M.RM.E.cfg \
--alternatives-from-evidence=true \
--gff3=on \
--uniqueGeneId=true \
--softmasking=1 \
--species=training_parameters \
scaffold-$SGE_TASK_ID/input_masked.fa > \
augustus_output/augustus_out_$SGE_TASK_ID.gff
#
echo = `date` job $JOB_NAME done
#
```
The job is submitted from the directory where all the scaffold\_n folders are located, and not from inside each folder.
Also, I'm saving all output files in a separate directory (augustus\_output), to facilitate port-processing.

#####Important points:
- \#$ -t 1-100 -tc 20:
	- the -t flag identifies the id number of the tasks to be run by this job. Here, we are running the tasks numbers 1 until 100.
	- the -tc flag indicates the number of consecutive tasks to be run. 
- \#$ -o augustus\_$TASK_ID.log
	- each task will have its own log file, identified by the variable $TASK\_ID.
- \--hintsfile=scaffold\_$SGE\_TASK\_ID/hints\_RM.E.gff3
	-  The variable $SGE\_TASK\_ID is used to identify each folder. In our example, this value would be anything between 1 and 100. For example, for task 10, the $SGE\_TASK\_ID value would equal to 10. 
	-  The same applies to scaffold-$SGE_TASK_ID/input_masked.fa and the output file.

####Combining the results.

Use the script `join_aug_pred.pl` from AUGUSTUS:

1. Concatenate all output files from augustus in numerical order:
`cat $(find . name "augustus_*.gff" | sort -V) > augustus.concat`

2. Join the results using the `join_aug_pred.pl`
`cat augustus.concat < join_aug_pred.pl > augustus_all.gff3`

3. Convert the Augustus GFF3 to EVM GFF3 (linear, without the sequences):
`/PATH/TO/EVM/EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl augustus_all.gff3 > augustus_final.gff3`

###References:
Haas et al. Automated eukaryotic gene structure annotation using EVidenceModeler and the Program to Assemble Spliced Alignments. Genome Biology 2008, 9:R7doi:10.1186/gb-2008-9-1-r7.

Kent, W. J. (2002). BLAT—the BLAST-like alignment tool. Genome research, 12(4), 656-664.

M. Stanke , O. Schöffmann , B. Morgenstern, S. Waack (2006)
Gene prediction in eukaryotes with a generalized hidden Markov model that uses hints from external sources 
BMC Bioinformatics 7, 62.

Simão, F. A., Waterhouse, R. M., Ioannidis, P., Kriventseva, E. V., & Zdobnov, E. M. (2015). BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs. Bioinformatics, 31(19), 3210-3212.

Smit, AFA, Hubley, R & Green, P. RepeatMasker Open-4.0.
2013-2015 <http://www.repeatmasker.org>.

Waterhouse, R. M., Seppey, M., Simão, F. A., Manni, M., Ioannidis, P., Klioutchnikov, G., ... & Zdobnov, E. M. (2017). BUSCO applications from quality assessments to gene prediction and phylogenomics. Molecular biology and evolution, 35(3), 543-548.

