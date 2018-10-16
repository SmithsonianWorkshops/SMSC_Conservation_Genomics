## Environment setup
#### Login to Hydra and install R libraries

```$ ssh <username>@hydra-login01.si.edu```

Now load the R module and open R

```$ module load tools/R/3.4.1```

```$ R```

Now we are going to install bioconductor packages that we will use during the differential expression portion of the analysis.

```> source("http://bioconductor.org/biocLite.R")```

Warning message:
> \>source("http://bioconductor.org/biocLite.R")
Warning in install.packages("BiocInstaller", repos = a["BioCsoft", "URL"]) :
  'lib = "/share/apps/tools/R/gcc/4.9.2/3.4.1/lib64/R/library"' is not writable
Would you like to use a personal library instead?  (y/n) 

`Reply "y"`

Warning Message:
 
>Would you like to create a personal library
~/R/x86_64-pc-linux-gnu-library/3.4
to install packages into?  (y/n)

`Reply "y"`

Then Continue: 

```> biocLite()```

```> biocLite('edgeR')```

```> biocLite('ctc')```

```> biocLite('Biobase')```

```> biocLite('ape')```

```> install.packages('gplots')```

Exit R

```> quit()```

When prompted to save your workspace image, reply 'n' and press enter:

```Save workspace image? [y/n/c]: n```


