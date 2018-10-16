## _De novo_ transcriptome assembly with Trinity


####Trinity _de novo_ assembly

Now we will start the Trinity run.

First open the [QSubGen web application](https://hydra-3.si.edu/tools/QSubGen).

- Choose medium time limit and reserve 6GB of memory.
- Select 'multi-thread' and choose 10 CPUs.
- In the modules field, type ```trinity``` and the select the module ```bioinformatics/trinity/2.1.1```.
- In job specific commands, type 
```
Trinity --seqType fq  \
          --left data/wt_SRR1582649_1.fastq,data/wt_SRR1582651_1.fastq,data/wt_SRR1582650_1.fastq,data/GSNO_SRR1582648_1.fastq,data/GSNO_SRR1582646_1.fastq,data/GSNO_SRR1582647_1.fastq \
          --right data/wt_SRR1582649_2.fastq,data/wt_SRR1582651_2.fastq,data/wt_SRR1582650_2.fastq,data/GSNO_SRR1582648_2.fastq,data/GSNO_SRR1582646_2.fastq,data/GSNO_SRR1582647_2.fastq \
          --max_memory 2G --min_contig_length 150 --CPU $NSLOTS --full_cleanup
 ```
- Choose a descriptive job name then click on ```Check if OK```
- If it passes, either save it and upload it to Hydra, or copy the text and paste it directly into your favority text editor.

Now save your text file into your ```/pool/genomics/<username>/RNAseq_workshop``` directory as ```trinity.job```.

Now submit your job with the command: ```qsub trinity.job```

Soon your transcriptome assembly will be finished!