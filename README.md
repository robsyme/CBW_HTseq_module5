---
layout: tutorial_page
permalink: /htseq_2020_module5_lab
title: HTSeq Lab 5
header1: Workshop Pages for Students
header2: Informatics on High-Throughput Sequencing Data Module 5 Lab
image: /site_images/CBW_High-throughput_icon.jpg
home: https://bioinformaticsdotca.github.io/htseq_2020
---

-----------------------

**This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.**

-----------------------

# CBW HT-seq Module 5 - Structural Variant Calling

 Created by Mathieu Bourgey, _Ph.D_, then modified by Pascale Marquis and Rob Syme


## Table of contents
1. [Introduction](#introduction)
2. [Original Setup](#setup)
2. [Align DNA with BWA-MEM](#align)
3. [Characterize the fragment size distribution](#fragments)
4. [Run DELLY to detect SVs](#delly)
5. [Setting up IGV for SV visualization](#IGV)
6. [Explore the SVs](#explore)
7. [Acknowledgements](#ackno)


## Introduction
<a name="introduction"></a>

The goal of this practical session is to identify structural variants (SVs) in a human genome by identifying both discordant paired-end alignments and split-read alignments that.

Discordant paired-end alignments conflict with the alignment patterns that we expect (i.e., concordant alignments) for the DNA library and sequencing technology we have used.

For example, given a ~500bp paired-end Illumina library, we expect pairs to align in F/R orientation and we expect the ends of the pair to align roughly 500bp apart. Pairs that align too far apart suggest a potential deletion in the DNA sample's genome. As you may have guessed, the trick is how we define what "too far" is --- this depends on the fragment size distribution of the data.


Split-read alignments contain SV breakpoints and consequently, then DNA sequences up- and down-stream of the breakpoint align to disjoint locations in the reference genome.


We will be working on a 1000 genome sample, NA12878. You can find the whole raw data on the 1000 genome website:
<http://www.1000genomes.org/data>

The dataset comes from the [Illumina Platinum Genomes Project](http://www.illumina.com/platinumgenomes/), which is a 50X-coverage dataset of the NA12891/NA12892/NA12878 trio. The raw data can be downloaded from the following [URL](http://www.ebi.ac.uk/ena/data/view/ERP001960).

NA12878 is the child of the trio while NA12891 and NA12892 are her parents.

<img src="./img/Pedigree.png?raw=true" alt="Pedigree" width="100%" />

For practical reasons we subsampled the reads from the sample because running the whole dataset would take way too much time and resources.
We're going to focus on the reads extracted from the chromosome 20.

## Original Setup
<a name="setup"></a>


### Software requirements

These are all already installed, but here are the original links.

  * [BWA-MEM](http://bio-bwa.sourceforge.net/)
  * [samtools](http://sourceforge.net/projects/samtools/)
  * [delly](https://github.com/tobiasrausch/delly)
  * [bcftools](https://samtools.github.io/bcftools/bcftools.html)
  * [LUMPY](https://github.com/arq5x/lumpy-sv)
  * [R](https://cran.r-project.org/)
  * [python](https://www.python.org/)


In this session, we will particularly focus on DELLY, a SV detection tool. DELLY is an integrated structural variant prediction method that can discover, genotype and visualize deletions, tandem duplications, inversions and translocations at single-nucleotide resolution in short-read massively parallel sequencing data. It uses paired-ends and split-reads to sensitively and accurately delineate genomic rearrangements throughout the genome.

If you are interested in DELLY, you can read the full manuscript [here](http://bioinformatics.oxfordjournals.org/content/28/18/i333.abstract).

### Environment setup

### Accessing a working node

When you log into the server, you are assigned to a "login" node (sometimes called a "head node"), which is shared by other users who are also logged in. As these nodes are a shared resouces, running computationally heavy workloads here can make the system unstable for everybody. In order to run your analysis in a stable environment without affecting other user you need to access a work node (sometimes called a "compute node"). Usually each job shoule be launched through the scheduler to run in a working environment, but our jobs in this workshop as are small and fast, so we can instead launch an interactive session on one of the work nodes by running:

```bash
salloc --mem 0 -n 8
```

The salloc command will assign us to a compute node and give us permission to use up to 8 cpus at a time. The interactive session will last for 1h, after which our session will end and we will be returned to the login node.

### Setting up the workspace

We'll create a directory in which to do the module 5 lab:

```bash
export WORK_DIR_M3=$HOME/workspace/HTseq/Module5
mkdir -p $WORK_DIR_M5
cd $WORK_DIR_M5
ln -s $HOME/CourseData/HT_data/Module5/* .

module load \
  mugqic/GenomeAnalysisTK/4.1.0.0 \
  mugqic/java/openjdk-jdk1.8.0_72 \
  mugqic/R_Bioconductor/3.5.0_3.7 \
  mugqic/samtools/1.10 \
  mugqic/bvatools/1.6 \
  mugqic/bcftools/1.6 \
  mugqic/tabix/0.2.6 \
  mugqic/Delly/0.7.8 \
  mugqic/bwa/0.7.17
```

***Note:***
    The `ln -s` command adds symbolic links of all of the files contained in the (read-only) `~/CourseData/HT_data/Module5` directory.

### Data files

The initial structure of your folders should look like this:


```console
[user@node1 ref]$ cd $WORK_DIR_M3
[user@node1 ref]$ tree -d *
bam
├── NA12878
├── NA12891
└── NA12892
fastq
reference
saved_results
├── Moleculo_bam
└── SVvariants
scripts
```

### Cheat sheets

* [Unix comand line cheat sheet](http://sites.tufts.edu/cbi/files/2013/01/linux_cheat_sheet.pdf)
* [commands file of this module](https://github.com/mbourgey/CBW_HTseq_module5/blob/master/scripts/commands.sh)



## Align DNA with BWA-MEM
<a name="align"></a>

This step has been done for you in the interest of time, but the commands are shown so that you can reproduce the results later. The advantage of using BWA-MEM in the context of SV discovery is that it produces both paired-end and split-read alignments in a single BAM output file. In contrast, prior to BWA-MEM, one typically had to use two different aligners in order to produce both high quality paired-end and split-read alignments.

In the alignment commands, note the use of the -M parameter to mark shorter split hits as secondary.

```bash
### NOTE: These preparatory commands have already been run for you and the outputs supplied.
###       The code is given here as a reference only.
###       You don't need to run this block.

## Index reference genome ##
bwa index reference/hg19.fa


## Align NA12878 ##
bwa mem -M -t 2 \
  reference/hg19.fa \
  fastq/NA12878_S1.chr20.20X.1.fq \
  fastq/NA12878_S1.chr20.20X.2.fq \
| samtools view -S -b - > bam/NA12878/NA12878_S1.chr20.20X.pairs.readSorted.bam


## Align NA12891 ##
bwa mem -M -t 2 \
  reference/hg19.fa \
  fastq/NA12891_S1.chr20.20X.1.fq \
  fastq/NA12891_S1.chr20.20X.2.fq \
| samtools view -S -b - > bam/NA12891/NA12891_S1.chr20.20X.pairs.readSorted.bam


## Align NA12892 ##
bwa mem -M -t 2 \
  reference/hg19.fa \
  fastq/NA12892_S1.chr20.20X.1.fq \
  fastq/NA12892_S1.chr20.20X.2.fq \
| samtools view -S -b - > bam/NA1289/NA12892_S1.chr20.20X.pairs.readSorted.bam
```


**Why should we mark shorter split hits as secondary ?** [solution](https://github.com/robsyme/CBW_HTseq_module5/blob/master/solutions/_aln1.md)


## Characterize the fragment size distribution
<a name="fragments"></a>

Before we can attempt to identify structural variants via discordant alignments, we must first characterize the insert size distribution

**What do we mean by "the insert size distribution"?** [solution](https://github.com/robsyme/CBW_HTseq_module5/blob/master/solutions/_fragment1.md)

**How can we use the fragment size distribution in SV detection ?** [solution](https://github.com/robsyme/CBW_HTseq_module5/blob/master/solutions/_fragment2.md)


The following script, taken from LUMPY extracts F/R pairs from a BAM file and computes the mean and stdev of the F/R alignments. It also generates a density plot of the fragment size distribution.

Let's calculate the fragment distribution for the three dataset:

```bash
# NA12878
mkdir -p SVvariants
samtools view bam/NA12878/NA12878_S1.chr20.20X.pairs.readSorted.bam \
| python scripts/pairend_distro.py \
  -r 101 \
  -X 4 \
  -N 10000 \
  -o SVvariants/NA12878_S1.chr20.20X.pairs.histo \
> SVvariants/NA12878_S1.chr20.20X.pairs.params


# NA12891
samtools view bam/NA12891/NA12891_S1.chr20.20X.pairs.readSorted.bam \
| python scripts/pairend_distro.py \
  -r 101 \
  -X 4 \
  -N 10000 \
  -o SVvariants/NA12891_S1.chr20.20X.pairs.histo \
> SVvariants/NA12891_S1.chr20.20X.pairs.params


# NA12892
samtools view bam/NA12892/NA12892_S1.chr20.20X.pairs.readSorted.bam \
| python scripts/pairend_distro.py \
  -r 101 \
  -X 4 \
  -N 10000 \
  -o SVvariants/NA12892_S1.chr20.20X.pairs.histo \
> SVvariants/NA12892_S1.chr20.20X.pairs.params
```

`-r` specifies the read length

`-X` specifies the number of stdevs from mean to extend

`-N` specfies the number of read to sample"

`-o` specifies the output file


Let's take a peak at the first few lines of the histogram file that was produced:

```console
[user@node1 Module5]$ head -n 10 SVvariants/NA12878_S1.chr20.20X.pairs.histo
0       0.0
1       0.00020046106043900973
2       0.00030069159065851457
3       0.00030069159065851457
4       0.00020046106043900973
5       0.00020046106043900973
6       0.00030069159065851457
7       0.00010023053021950487
8       0.00010023053021950487
9       0.00010023053021950487
```

Let's use R to plot the fragment size distribution. First, launch R from the command line.

```bash
R
```

Now, within R, execute the following commands:

```R
size_dist <- read.table('SVvariants/NA12878_S1.chr20.20X.pairs.histo')
pdf(file = "SVvariants/fragment.hist.pdf")
layout(matrix(1:3))
plot(size_dist[,1], size_dist[,2], type='h', main="NA12878 insert size")
size_dist <- read.table('SVvariants/NA12891_S1.chr20.20X.pairs.histo')
plot(size_dist[,1], size_dist[,2], type='h', main="NA12891 insert size")
size_dist <- read.table('SVvariants/NA12892_S1.chr20.20X.pairs.histo')
plot(size_dist[,1], size_dist[,2], type='h', main="NA12892 insert size")
dev.off()
quit("no")
```

At this point, you should have the following files:

```console
[user@node1 Module5]$ tree SVvariants/
SVvariants/
├── fragment.hist.pdf
├── NA12878_S1.chr20.20X.pairs.histo
├── NA12878_S1.chr20.20X.pairs.params
├── NA12891_S1.chr20.20X.pairs.histo
├── NA12891_S1.chr20.20X.pairs.params
├── NA12892_S1.chr20.20X.pairs.histo
└── NA12892_S1.chr20.20X.pairs.params
```

To look at the fragment.hist.pdf file we generated, we'll need to copy it from the cluster to our local computer. If you used the same directory names as I have, the pdf should be located at `~/workspace/HTseq/Module5/SVvariants/fragment.hist.pdf`. If you have used different directory names, be sure to change the path in the command below.

On your *local* computer, run:

```bash
rsync YOURUSERNAMEHERE@login1.CBW.calculquebec.cloud:~/workspace/HTseq/Module5/SVvariants/fragment.hist.pdf .
```

This will copy the file onto your local computer where you can view the distributions. Spend some time thinking about what this plot means for identifying discordant alignments.

**What does the mean fragment size appear to be? Are all 3 graphs the same?** [solution](https://github.com/robsyme/CBW_HTseq_module5/blob/master/solutions/_insert1.md)


## SV detection
<a name="delly"></a>

For germline SVs, calling is done by sample or in small batches to increase SV sensitivity & breakpoint precision.

Here are the steps adapted from the DELLY [readme](https://github.com/tobiasrausch/delly/blob/master/README.md).


### First call

We've used Delly to call SVs.

```bash
# NA12878
delly call \
 --genome reference/hg19.fa \
 --exclude reference/hg19.excl \
 --outfile SVvariants/NA12878.bcf \
 bam/NA12878/NA12878_S1.chr20.20X.pairs.posSorted.bam

## NA12891
delly call \
 --genome reference/hg19.fa \
 --exclude reference/hg19.excl \
 --outfile SVvariants/NA12891.bcf \
 bam/NA12891/NA12891_S1.chr20.20X.pairs.posSorted.bam

# NA12892
delly call \
 --genome reference/hg19.fa \
 --exclude reference/hg19.excl \
 --outfile SVvariants/NA12892.bcf \
 bam/NA12892/NA12892_S1.chr20.20X.pairs.posSorted.bam
```


The `.bcf` output files are compressed (binary) vcf files. To inspect the output from the NA12878 sample, run:

```bash
bcftools view SVvariants/NA12878.bcf | less -S
```


***Cheat:*** These commands take a little under 5 minutes each to run. If these commands are taking too long, simply run the command `cp saved_results/SVvariants/NA128*[128].bc* SVvariants/`

**How many variants delly found in each sample ?** [solution](https://github.com/robsyme/CBW_HTseq_module5/blob/master/solutions/_vcf1.md)

**How many variants by SV type are found in each sample ?** [solution](https://github.com/robsyme/CBW_HTseq_module5/blob/master/solutions/_vcf4.md)

### Merge calls

We need to merge the SV sites into a unified site list:

```bash
delly merge \
 --minsize 500 \
 --maxsize 1000000 \
 --bp-offset 500 \
 --rec-overlap 0.5 \
 --outfile SVvariants/sv.bcf \
 SVvariants/NA12878.bcf SVvariants/NA12891.bcf SVvariants/NA12892.bcf
```

Look at the output:

```bash
bcftools view SVvariants/sv.bcf | less -S
```


**How many variants by SV type are found ?** [solution](https://github.com/robsyme/CBW_HTseq_module5/blob/master/solutions/_vcf5.md)

**What can you notice something different from the individual bcf ?** [solution](https://github.com/robsyme/CBW_HTseq_module5/blob/master/solutions/_vcf2.md)


### Re-genotype in all samples

We need to re-genotype the merged SV site list across all samples. This can be run in parallel for each sample.

```bash
# NA12878
delly call \
 --genome reference/hg19.fa \
 --vcffile SVvariants/sv.bcf \
 --exclude reference/hg19.excl \
 --outfile SVvariants/NA12878.geno.bcf \
 bam/NA12878/NA12878_S1.chr20.20X.pairs.posSorted.bam &

# NA12891
delly call \
 --genome reference/hg19.fa \
 --vcffile SVvariants/sv.bcf \
 --exclude reference/hg19.excl \
 --outfile SVvariants/NA12891.geno.bcf \
 bam/NA12891/NA12891_S1.chr20.20X.pairs.posSorted.bam &

# NA12892
delly call \
 --genome reference/hg19.fa \
 --vcffile SVvariants/sv.bcf \
 --exclude reference/hg19.excl \
 --outfile SVvariants/NA12892.geno.bcf \
 bam/NA12892/NA12892_S1.chr20.20X.pairs.posSorted.bam &

 wait
```

Look at the output:

```bash
bcftools view SVvariants/NA12878.geno.bcf | less -S
```

We now have our good candidate list genotype for each individual.


### Merge the new calls

Merge all re-genotyped samples to get a single VCF/BCF using bcftools merge. Also index the resulting file and create vcf file for visualization.

```bash
bcftools merge \
 --output-type b \
 --output SVvariants/merged.bcf \
 SVvariants/NA12878.geno.bcf SVvariants/NA12891.geno.bcf SVvariants/NA12892.geno.bcf

bcftools index SVvariants/merged.bcf

bcftools view SVvariants/merged.bcf > SVvariants/merged.vcf

bgzip -c SVvariants/merged.vcf > SVvariants/merged.vcf.gz

tabix -fp vcf SVvariants/merged.vcf.gz
```


**Do you know how to look at the resulting file?** [solution](https://github.com/robsyme/CBW_HTseq_module5/blob/master/solutions/_vcf3.md)


## Setting up IGV for SV visualization
<a name="IGV"></a>

Launch IGV and load the merged calls and the germline calls using `File -> Load from URL` using:


 * https://datahub-39-cm2.p.genap.ca/HTseq/Module5/SVvariants/merged.vcf.gz
 * https://datahub-39-cm2.p.genap.ca/HTseq/Module5/SVvariants/merged.vcf.gz.tbi


Now load the bam files in the same way using:

 * https://datahub-39-cm2.p.genap.ca/HTseq/Module5/bam/NA12878/NA12878_S1.chr20.20X.pairs.posSorted.bam
 * https://datahub-39-cm2.p.genap.ca/HTseq/Module5/bam/NA12891/NA12891_S1.chr20.20X.pairs.posSorted.bam
 * https://datahub-39-cm2.p.genap.ca/HTseq/Module5/bam/NA12892/NA12892_S1.chr20.20X.pairs.posSorted.bam


Navigate to the following location to see a deletion: `chr20:31,308,410-31,315,294`

You should see something like this:

<img src="./img/deletion.png?raw=true" alt="Deletion" width="100%" />

You can try to configure IGV such that we can more clearly see the alignments that support the SV prediction.

***Note*** - Do you remember how to make sure that the alignments are colored by insert size and orientation?*


## Explore the SVs
<a name="explore"></a>

**Is the variant at chr20:31,310,769-31,312,959 found in each member of the trio?** [solution](https://github.com/robsyme/CBW_HTseq_module5/blob/master/solutions/_igv1.md)

**What are the genotypes for each member of the trio at the locus chr20:61,721,523-61,728,495  (e.g., hemizygous, homozygous)?** [solution](https://github.com/robsyme/CBW_HTseq_module5/blob/master/solutions/_igv2.md)

**What about the variant at chr20:52,632,182-52,664,108?** [solution](https://github.com/robsyme/CBW_HTseq_module5/blob/master/solutions/_igv3.md)

Now load the bam files in

 * https://datahub-39-cm2.p.genap.ca/HTseq/Module5/saved_results/Moleculo_bam/NA12878.molelculo.chr20.bam

**Does the evidence in the Moleculo track mimic the evidence in the Illumina track for NA12878?** [solution](https://github.com/robsyme/CBW_HTseq_module5/blob/master/solutions/_igv4.md)

**What about chr20:18,953,476-18,957,998?** [solution](https://github.com/robsyme/CBW_HTseq_module5/blob/master/solutions/_igv5.md)

Continue exploring the data!


## Acknowledgements
<a name="ackno"></a>

This module is heavily based on a previous module prepared by Aaron Quinlan and Guillaume Bourque.

