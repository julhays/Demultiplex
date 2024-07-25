# BI62 Demultiplexing - Lab Notebook
Jules Hays
## Due: 8/2/24

Python version: 3.12

Environments used: bgmp_py312 (matplotlib 3.9.1)

---

## Assignment the First
Goal: look through a lane of sequencing generated from the 2017 BGMP cohort’s library preps and determine the level of index swapping and undetermined index-pairs, before and after quality filtering of index reads
### 7/25/24

I logged into Talapas and cloned the Demultiplex repo in new directory for this assignment, located at
```
/projects/bgmp/jkhay/bioinfo/Bi622/Demultiplex
```

I navigated to the shared Talapas directory to obtain the sequencing data at
```
/projects/bgmp/shared/2017_sequencing/
```

This directory contains 4 files with the fastq sequencing data:
```
1294_S1_L008_R1_001.fastq.gz
1294_S1_L008_R2_001.fastq.gz
1294_S1_L008_R3_001.fastq.gz
1294_S1_L008_R4_001.fastq.gz
```

It also contains a text file called ```indexes.txt``` that has a tab separated file with the information for each index. It is formatted in the following fashion:
```
sample  group   treatment       index   index sequence
1       2A      control B1      GTAGCGTA
...
```

It would be very bad to copy or unzip these files so don't do it.


### Part 1 - Quality Score Distribution per-nucleotide
Goal: Perform initial data exploration and plot the mean quality score at each nucletide position for each file.

For initial data exploration I will start an interactive session in Talapas:
```
srun -A bgmp -p bgmp --mem=100gb -c 8 --pty bash
```

Initial Data Exploration

Get the size of the files:
```
$ ls -lah

-rw-r-xr--+  1 coonrod  is.racs.pirg.bgmp  20G Jul 30  2018 1294_S1_L008_R1_001.fastq.gz
-rw-r-xr--+  1 coonrod  is.racs.pirg.bgmp 2.6G Jul 30  2018 1294_S1_L008_R2_001.fastq.gz
-rw-r-xr--+  1 coonrod  is.racs.pirg.bgmp 2.8G Jul 30  2018 1294_S1_L008_R3_001.fastq.gz
-rw-r-xr--+  1 coonrod  is.racs.pirg.bgmp  21G Jul 30  2018 1294_S1_L008_R4_001.fastq.gz
```

R1 and R4 files are about 20-21 GB ZIPPED :open_mouth:, and R2 and R3 are about 2.6-2.8 GB ZIPPED. Dang I best not unzip these.

Look at the file structure:
```
$ zcat 1294_S1_L008_R1_001.fastq.gz | head 
```
Alright, it looks like a normal fastq file.

Compare the first entry across the 4 files:
```
$ zcat 1294_S1_L008_R1_001.fastq.gz | head -4
$ zcat 1294_S1_L008_R2_001.fastq.gz | head -4
$ zcat 1294_S1_L008_R2_001.fastq.gz | head -4
$ zcat 1294_S1_L008_R4_001.fastq.gz | head -4
```
The first entry in each file have nearly the name header: 

@K00337:83:HJKJNBBXX:8:1101:1265:1191 X:N:0:1

But the X is replaced with the R number in the file name. So each entry number is linked across the 4 files. It also appears that matching R2 and R3s are the reverse compliment of each other.

Count the number of lines/records in the files:
```
$ zcat 1294_S1_L008_R1_001.fastq.gz | wc -l
$ zcat 1294_S1_L008_R2_001.fastq.gz | wc -l
$ zcat 1294_S1_L008_R3_001.fastq.gz | wc -l
$ zcat 1294_S1_L008_R4_001.fastq.gz | wc -l
```
Each of the files has the same number of lines too.

1452986940 lines / 4 = 363,246,735 records

There are 363,246,735 records per file.

Find the read length of each file:
```
$ zcat 1294_S1_L008_R1_001.fastq.gz | sed -n '2~4p' | awk '{print length($0)}' | uniq
$ zcat 1294_S1_L008_R2_001.fastq.gz | sed -n '2~4p' | awk '{print length($0)}' | uniq
$ zcat 1294_S1_L008_R3_001.fastq.gz | sed -n '2~4p' | awk '{print length($0)}' | uniq
$ zcat 1294_S1_L008_R4_001.fastq.gz | sed -n '2~4p' | awk '{print length($0)}' | uniq
```

R1 and R4 have a read length of 101 nts.
R2 and R3 have a read length of 8 nts.

Based on the lengths of the reads in the files, it appears that R1 and R4 are biological reads and R2 and R3 are index reads.

Determine the Phred score encoding:
```
$ zcat 1294_S1_L008_R1_001.fastq.gz | sed -n '4~4p' | head
$ zcat 1294_S1_L008_R2_001.fastq.gz | sed -n '4~4p' | head
$ zcat 1294_S1_L008_R3_001.fastq.gz | sed -n '4~4p' | head
$ zcat 1294_S1_L008_R4_001.fastq.gz | sed -n '4~4p' | head
```
The Phred score encoding appears to be Phred+33 because the quality scores contain symbols like '<' and '#' (ASCII values of 60 and 35, respectively). If the encoding was Phred+64 these characters would not be present because the lowest score would be '@' (64).

Below is a summary table of my intial data exploration:

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | forward read | 101 | Phred+33 |
| 1294_S1_L008_R2_001.fastq.gz | barcode 1 | 8 | Phred+33 |
| 1294_S1_L008_R3_001.fastq.gz | barcode 2 | 8 | Phred+33 |
| 1294_S1_L008_R4_001.fastq.gz | reverse read | 101 | Phred+33 |

Now I need to write a Python script to graph the distribution of quality scores for each file. This script is called ```




4 files
R2 and R3 should be rev comps, if not uk or index diff

if comp, add to header


I need to record my answers in ```Answers.md```.