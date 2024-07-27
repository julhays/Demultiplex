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
/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz
/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz
/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz
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

Now I need to write a Python script to graph the distribution of quality scores for each file. This script is called ```qual_dist.py``` and is in the ```Assignment_the_first``` directory. I added my bioinfo.py to this repo as well so I can utilize the functions in my script.

I added an argparse that takes in the following:
* -f - Specify the path of the input fastq file
* -n - Specify the read length of each read in the file

The file will output a png file that had the name ```RX.png``` where X is the R number of the input file.

To run my script, I made sure I was on an interactive Talapas session and activated my matplotlib environment, then ran my script:
```
srun -A bgmp -p bgmp --mem=100gb -c 8 --pty bash
conda activate bgmp_py312
./qual_dist.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -n 101
```

I made a test file that is just the first 12 lines of the R1 file and that ran through successfully.

It takes a LONG TIME to run the actual files... :sob:

Ok it took so long that I am writing an SBATCH script to run each file. There's a script for each file called ```R#_dist.sh``` where # is the R value, and I will run them all at once and hope it works :pray:

```
$ sbatch R1_dist.sh 
Submitted batch job 7638320
$ sbatch R2_dist.sh 
Submitted batch job 7638321
$ sbatch R3_dist.sh 
Submitted batch job 7638322
$ sbatch R4_dist.sh 
Submitted batch job 7638323
```

### 7/26/24
### Still Part 1

It worked! Here's a summary of the runtime output.

| File name | run time | % CPU |
|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | 3:09:23.00 | 99% |
| 1294_S1_L008_R2_001.fastq.gz | 23:02.98 | 99% |
| 1294_S1_L008_R3_001.fastq.gz | 22:22.43 | 99% |
| 1294_S1_L008_R4_001.fastq.gz | 3:44:49.00 | 99% |

My histograms are located in the ```Assignment-the-first``` directory. They have the name ```R#_dist.png``` where # is the R numbers.

It is now occuring to be that I made a scatterplot and not a histogram, but I think the scatterplot looks better. I need to ask if this is ok.

* Quality score cutoff?

Based on the graphs, it looks like a good cutoff for the biological read pairs is 30 because the average is above 30 for all positions, and slightly above 30 for the lowest quality positions. I think 30 is also good for the index reads.

26 for index?

* How many indicies have N base calls?

```
$ for file in /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz; do zcat $file | sed -n '2~4p' | grep -c 'N'; done

3976613
3328051
```
There are 3976613 index 1 barcodes with an N and 3328051 index 2 barcodes with an N.

Questions:
* It is now occuring to be that I made a scatterplot and not a histogram, but I think the scatterplot looks better. I need to ask if this is ok.
* 30 reasonable cutoff
* sum the N indicies or report both files (if both barcodes for a certain sample have N, would this count as 1 or 2)

I need to record my answers in ```Answers.md```.

### Part 2

1. Define the problem

Given 4 FASTQ files (2 with biological reads and 2 with index reads) and a list of acceptable matched indexes, write a script that demultiplexes the reads and reports index-hopping and low quality index reads.

2. Describe output

The script will output the following files:
* an R1 and R2 file per matching index pair
* an R1 and R2 file with hopped index pairs
* an R1 and R2 file with unknown or low quality index pairs (don't match the acceptable index list or contain an N)

All files should have the sequence of the index pairs in the header of both reads

The script will also output some statistics:
* the number of read pairs with matching indexes per index pair
* the number of read pairs with index hopping per index pair
* the number of read pairs with unknown indexes

3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).


4. Pseudocode
```
#make list of acceptable indexes called indexes
while read indexes.txt:
    for line in indexes.txt:
        extract index sequence and add to indexes

#create dictionaries to keep track of statistics as it goes
matched = {} # stores the acceptable indexes as keys and the number of correctly paired occurances as values
hopped = {} # stores the acceptable indexes as keys and the number of incorrectly paired occurances as values
unknown = 0 # integer count of the number of unknown pairs

#parse through input files and sort into respective output files
#while doing this keep count of the occurance of a match or mismatch for each index pair
while read R1, R2, R3, R4 file:
    while True:

        extract a fastq record (4 lines) for each file and store in a list

        if all the lines in each file are empty:
            break

        new header 1, new header 2 = barcodes_to_header(first line of R1 record, first line of R4 record, second line of R2 record, second line of R3 record) 
        
        if line 2 of R2 not in indexes or line 2 of R3 not in indexes:
            write R1 and R4 record with new headers to unknown_R1.fq and unknown_R2.fq
            increment unknown by 1

        low quality = low_qual(forth line of R2 record, forth line of R3 record)

        if low quality:
            write R1 and R4 record with new headers to unknown_R1.fq and unknown_R2.fq
            increment unknown by 1

        reverse compliment = rev_comp(second line of R2 record, second line of R3 record)

        if not reverse compliment:
            write R1 and R4 record with new headers to hopped_R1.fq and hopped_R2.fq
            if second line of R2 record in hopped dictionary:
                increment value by 1
            else:
                insert second line of R2 record into hopped dictionary with a value of 1

        else:
            write R1 and R4 record with new headers to {second line of R2 record}_R1.fq and {second line of R2 record}_R2.fq
             if second line of R2 record in matched dictionary:
                increment value by 1
            else:
                insert second line of R2 record into matched dictionary with a value of 1

#print statistics
print "Index    Matched    Hopped"
for barcode in indexes:
    print barcode, value in matched, and value in hopped if they exist (tab separated)

print "Number of read pairs with unknown indexes"
print unknown

```

5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
```
def barcodes_to_header(R1_head: str, R4_head: str, R2_seq: str, R3_seq: str) -> str, str:
    ```Takes in the R1 and R4 header and appends the barcode sequences 
    listed in R2 and the reverse complement of R3 to the end of the file header to make it headerBARCODE1-rcBARCODE2```
    return R1_head, R4_head
Input: ('@K00337:83:HJKJNBBXX:8:1101:1265:1191 1:N:0:1', '@K00337:83:HJKJNBBXX:8:1101:1265:1191 4:N:0:1', 'AAAAAAAA', 'TTTTTTTT')
Expected output: '@K00337:83:HJKJNBBXX:8:1101:1265:1191 1:N:0:1AAAAAAAA-AAAAAAAAA', '@K00337:83:HJKJNBBXX:8:1101:1265:1191 4:N:0:1AAAAAAAA-AAAAAAAAA'

def low_qual(R2_qual: str, R3_qual: str) -> bool:
    ```Takes in the barcode quality score lines in R2 and R3 and returns True if 
    the average quality score of either does not meet or exceed the cutoff```
    return True or return False
Input: ('!!!!!!!!', '!!!!!!!!')
Expected output: True

def rev_comp(R2_seq: str, R3_seq: str) -> bool:
    ```Takes in the barcode sequences listed in R2 and R3 and returns True if 
    they are reverse compliments and False if not```
    return True or return False
Input: ('GTAGCGTA', 'TACGCTAC')
Expected output: True

```

NOTES:

check other things first then check if complementary (not in list or have N)

R2 and R3 should be rev comps, if not uk (N in index or not in lis of indexes) or hopped
if complementary:
function to see if passes qual score thresh and if it doesnt pass then get added to uk
if it does pass it add R1 to file_R1.fq and R4 to file_R2.fq

if comp, add to header

R2 will match lists of indexes, R3 will be rev complement



