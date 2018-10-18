## _De novo_ transcriptome assembly with Trinity


####Trinity _de novo_ assembly

Now we will start the Trinity run.

First open the [QSubGen web application](https://hydra-4.si.edu/tools/QSubGen).

- Choose medium time limit and reserve 6GB of memory.
- Select 'multi-thread' and choose 10 CPUs.
- In the modules field, type ```trinity``` and the select the module ```bioinformatics/trinity/2.6.6```.
- In job specific commands, type:
          
      	Trinity --seqType fq 
      	--left data/RNA_Eye_1_val_1.fastq \
      	--right data/RNA_Eye_2_val_2.fastq \
      	--max_memory 60G --CPU $NSLOTS --full_cleanup

- Choose a descriptive job name then click on 
- ```Check if OK```
- If it passes, either save it and upload it to Hydra, or copy the text and paste it directly into your favority text editor.

Now save your text file into your ```/pool/genomics/<username>/RNAseq_SMSC``` directory as ```trinity.job```.

Now submit your job with the command: ```qsub trinity.job```

Soon your transcriptome assembly will be finished!

```
Results of all tutorials can be found here:
/data/genomics/workshops/smsc/RNA_Seq/SMSC_results.tar.gz
```