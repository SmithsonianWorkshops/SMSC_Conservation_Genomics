## Differential Gene Expression with edgeR

_Note: First we need to install one more ```R``` package. To do so follow these instructions:_

Once you are logged in to hydra:

```
$ module load bioinformatics/trinity/2.6.6
$ R
```

Now that you are in ```R```, you need to load the biocLite function again:

```
> source("http://bioconductor.org/biocLite.R")
> biocLite()
> biocLite('qvalue')
> biocLite('fastcluster')
```

For both the last two commands, you will be prompted to update your other libraries. Respond with ```n```.

Then exit

```
> quit()
```

When promted to "Save workspace image? [y/n/c]:" reply `n`. 

####Create sample text file

You will need to create a tab delimited text file containing information for your different samples. In the first column, you will have the name of the condition and, in the second column, you will enter the name of the sample. Start a new text file, samples.txt, with nano.

```
$ nano samples.txt
```

Enter the following (keep in mind that the condition and names should be separated by tabs!).

```
GSNO	GSNO_SRR1582648.RSEM
GSNO	GSNO_SRR1582647.RSEM
GSNO	GSNO_SRR1582646.RSEM
WT	wt_SRR1582651.RSEM
WT	wt_SRR1582649.RSEM
WT	wt_SRR1582650.RSEM
```

Since software can be very picky about whether you specified config files correctly, it is sometimes good to check that you did, indeed, enter the correct characters. You can view special characters with:

```
$ cat -te samples.txt
```

If your file was specified correctly, it should look like this:

```
GSNO^IGSNO_SRR1582648.RSEM$
GSNO^IGSNO_SRR1582647.RSEM$
GSNO^IGSNO_SRR1582646.RSEM$
WT^Iwt_SRR1582651.RSEM$
WT^Iwt_SRR1582649.RSEM$
WT^Iwt_SRR1582650.RSEM$
```

```^I``` characters are tabs and ```$``` characters are newlines. Make sure that your text file looks like the example above when using ```cat -te```. If it doesn't, you'll need to edit it until it does.

####Detect differentially expressed transcripts in ```edgeR```

Now we are going to use the ```run_DE_analysis.pl``` script that is included with the ```Trinity``` package to detect differentially expressed transcripts in edgeR. Note that this will only work if the R packages from ```Environment setup.md``` were installed properly.

Create a new job file, and select the short queue and 2GB of RAM. Load the Trinity module. The command will look like this:

```
run_DE_analysis.pl \
      --matrix Trinity_trans.isoform.counts.matrix \
      --samples_file samples.txt \
      --method edgeR \
      --output edgeR_trans
```

Save the job file to your ```/pool/genomics/<username>/RNAseq_SMSC``` directory and submit it to the cluster. Once it is finished, there will be a new directory called ```edgeR_trans```. Take a look at its contents:

```
$ ls -lh edgeR_trans
```

There should be four files in the directory:

```
-rw-rw-r-- 1 gonzalezv gonzalezv  21K Oct 17 01:28 Trinity_trans.isoform.counts.matrix.GSNO_vs_WT.edgeR.count_matrix
-rw-rw-r-- 1 gonzalezv gonzalezv  61K Oct 17 01:28 Trinity_trans.isoform.counts.matrix.GSNO_vs_WT.edgeR.DE_results
-rw-rw-r-- 1 gonzalezv gonzalezv  13K Oct 17 01:28 Trinity_trans.isoform.counts.matrix.GSNO_vs_WT.edgeR.DE_results.MA_n_Volcano.pdf
-rw-rw-r-- 1 gonzalezv gonzalezv 1.4K Oct 17 01:28 Trinity_trans.isoform.counts.matrix.GSNO_vs_WT.GSNO.vs.WT.EdgeR.Rscript
```

The file ```Trinity_trans.isoform.counts.matrix.GSNO_vs_WT.edgeR.DE_results``` contains the results from comparing the GSNO condition to the wt condition. Take a look:

```
$ head edgeR_trans/Trinity_trans.isoform.counts.matrix.GSNO_vs_WT.edgeR.DE_results | column -t
```


```
sampleA              sampleB logFC  logCPM             PValue            FDR
TRINITY_DN587_c0_g1  GSNO    WT     -5.91504637307699  13.4818421437055  5.78409571193267e-55  3.74809402133237e-52
TRINITY_DN283_c0_g1  GSNO    WT     -8.96426024608977  13.1156162622867  1.55356586112296e-53  5.03355339003839e-51
TRINITY_DN586_c0_g1  GSNO    WT     -4.85424799859937  13.20852422536    8.59574378939845e-51  1.85668065851007e-48
TRINITY_DN300_c0_g1  GSNO    WT     -2.65263866887907  14.1251864925867  3.86504266443156e-50  6.26136911637913e-48
TRINITY_DN545_c0_g1  GSNO    WT     -2.66907551183065  14.0393039359086  1.33959062965647e-47  1.73610945603479e-45
TRINITY_DN425_c0_g1  GSNO    WT     -8.94862298901304  13.1002798609347  2.14783943672528e-47  2.3196665916633e-45
TRINITY_DN151_c0_g1  GSNO    WT     -5.63202406802949  12.8612509729179  2.93794317454296e-44  2.7196959672912e-42
TRINITY_DN283_c1_g1  GSNO    WT     -6.92986764039549  12.8711053889033  5.05433879863526e-34  4.09401442689456e-32
TRINITY_DN41_c0_g1   GSNO    WT     -5.22098506442703  12.8317674611305
```

As you can see, ```edgeR``` calculates log fold change (```logFC```), the log counts per million (```logCPM```), the p-value from the exact test (```PValue```), and the false discovery rate (```FDR```). 

_Note: Since there is no header for gene name, the headers are shifted one column to the right, i.e. ```logFC``` should be over the first column of floating point numbers._

```edgeR``` also generated MA and Volcano plots for these data. We will now download them to our computer. If you are using Mac or Linux, we will do this with the ```scp``` command. Open a new terminal window and ```cd``` to the directory that you wish to download the files to. On Mac, I often download to my ```Downloads``` directory. You can go there with:

```
$ cd ~/Downloads
```

Now download the plot:

```
$ scp <username>@hydra-login01.si.edu:/pool/genomics/<username>/RNAseq_SMSC/edgeR_trans/*.pdf .
```

Go ahead and open it to examine its contents.

![Volcano Plot](Trinity_trans.counts.matrix.GSNO_vs_WT.edgeR.DE_results.MA_n_Volcano.pdf)

The points that are in red are determined to be significant with an ```FDR``` <= 0.05. To read more about these tests, you can follow the citations on the [edgeR bioconductor page](https://bioconductor.org/packages/release/bioc/html/edgeR.html).

You might wonder what you can do with these data. Luckily, Trinity also includes scripts to extract differentially expressed transcripts and to create heatmaps.

Change directories into your ```edgeR_trans``` directory:

```
$ cd edgeR_trans
```

Now we will extract any transcript that is 4-fold differentially expressed between the two conditions at a significance of ```<= 0.001```.

Make another job file and choose the short queue and reserve the default RAM (1GB). Load the ```bioinformatics/trinity/2.6.6``` module. Your command will be:

```
analyze_diff_expr.pl \
      --matrix ../Trinity_trans.isoform.counts.matrix \
      --samples ../samples.txt \
      -P 1e-3 -C 2 
```

This command will filter transcripts based on pvalue of less than Several files will be written as a part of this job. One is called ```diffExpr.P1e-3_C2.matrix```. You can count the number of differentially expressed genes at this threshold by counting the number of lines:

```
$ wc -l diffExpr.P1e-3_C2.matrix
```

You should subtract 1 from the number since there is a header line.

This script also generates a heatmap that compares the differentially expressed transcripts. The file is called, ```diffExpr.P1e-3_C2.matrix.log2.centered.genes_vs_samples_heatmap.pdf```.

Download that file and examine it on your computer.

_Hint: you can use ```scp``` as above. Or you can use a GUI interface like Filezilla/Cyberduck._

Now examine the heatmap

![Differential Expression Heatmap](diffExpr.P1e-3_C2.matrix.log2.centered.genes_vs_samples_heatmap.pdf)

You can use the heatmap to compare the two conditions. The left columns with the turquoise line on top are those under wt and the right columns under the red line are under GSNO. Upregulated expression is in yellow and downregulated expression is in purple. This is a nice visual way to compare expression across conditions.

####View transcript clusters

You can also cut the dendrogram to view transcript clusters that share similar expression profiles. To do this, run the following command into a job file. Be sure to load the ```bioinformatics/trinity/2.6.6``` module and choose a serial job with 1GB of RAM:

```
define_clusters_by_cutting_tree.pl --Ptree 60 -R diffExpr.P1e-3_C2.matrix.RData
```

You should have a new output that looks like the following graph, which shows transcripts with similar expression profiles:

![My cluster plots](my_cluster_plots.pdf)

####Now run on genes

Now we will run differential expression analysis on the gene level. This will be very similar to the isoform analysis, but we will use the follow command:

```
run_DE_analysis.pl \
      --matrix Trinity_trans.genes.counts.matrix\
      --samples_file samples.txt \
      --method edgeR \
      --output edgeR_gene
```

You can also run the other downstream analyses, but you should replace ```Trinity_trans``` with ```Trinity_genes``` in the commands.