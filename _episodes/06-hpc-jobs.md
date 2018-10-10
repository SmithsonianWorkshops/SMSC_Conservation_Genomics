---
title: "Running jobs on HPC scheduler"
teaching: 20
exercises: 20
questions:
- What kinds of jobs are too big for running on the log-in node?
objectives:
- Learn how to submit a job to the HPC scheduler.
keypoints:
- The HPC is great.
---

## SI Computing Cluster: Hydra

PDF presentation about Hydra: 

## Modules

Modules are a way to set your system paths and environmental variables for use of a particular program or pipeline.

You can view available modules with the command:

~~~
$ module avail
~~~
{: .bash}

This will output all modules installed on Hydra. A complete list can also be found at [https://www.cfa.harvard.edu/~sylvain/hpc/module-avail.html](https://www.cfa.harvard.edu/~sylvain/hpc/module-avail.html).

There are additional module commands to find out more about a particular module. Let's use the bioinformatics program *RAxML* to test these out.

~~~
$ module whatis bioinformatics/raxml
~~~
{: .bash}

~~~
bioinformatics/raxml : System paths to run RAxML 8.2.11
~~~
{: .output}

~~~
$ module help bioinformatics/raxml
~~~
{: .bash}

~~~
----------- Module Specific Help for 'bioinformatics/raxml/8.2.11' ---------------------------


Purpose
-------
This module file defines the system paths for RAxML 8.2.11
The compiled binaries that you can now call are:
serial version: raxmlHPC-SSE3
pthreads version: raxmlHPC-PTHREADS-SSE3

Documentation
-------------
http://sco.h-its.org/exelixis/web/software/raxml/index.html

<- Last updated: Tue Oct  3 08:21:16 EDT 2017 ->
~~~
{: .output}

~~~
$ module show bioinformatics/raxml
~~~
{: .bash}

~~~
-------------------------------------------------------------------
/share/apps/modulefiles/bioinformatics/raxml/8.2.11:

module-whatis	 System paths to run RAxML 8.2.11
prepend-path	 PATH /share/apps/bioinformatics/raxml/gcc/4.9.2/8.2.11/bin
module		 load gcc/4.9.2
conflict	 bioinformatics/raxml
conflict	 bioinformatics/raxml/8.2.11-mpi
conflict	 bioinformatics/raxml/8.2
conflict	 bioinformatics/raxml/8.2-IB
conflict	 bioinformatics/raxml/8.2.7
conflict	 bioinformatics/raxml/8.2.9
-------------------------------------------------------------------
~~~
{: .output}

Ok, now let's actually load RAxML.
~~~
$ module load bioinformatics/raxml
~~~
{: .bash}

Nothing happens, but let's run a quick command to show that we have RAxML loaded properly.

~~~
$ raxmlHPC-SSE3 -help
~~~
{: .bash}

~~~
This is RAxML version 8.2.11 released by Alexandros Stamatakis on June 2017.

With greatly appreciated code contributions by:
Andre Aberer      (HITS)
Simon Berger      (HITS)
Alexey Kozlov     (HITS)
Kassian Kobert    (HITS)
David Dao         (KIT and HITS)
Sarah Lutteropp   (KIT and HITS)
Nick Pattengale   (Sandia)
Wayne Pfeiffer    (SDSC)
Akifumi S. Tanabe (NRIFS)
Charlie Taylor    (UF)


Please also consult the RAxML-manual

Please report bugs via the RAxML google group!
Please send us all input files, the exact invocation, details of the HW and operating system,
as well as all error messages printed to screen.

....

~~~
{: .output}


### Submitting a job

### Checking on your job

### Transferring files

### Interactive queue





