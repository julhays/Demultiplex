# BI62 Demultiplexing - Lab Notebook
Jules Hays
## Due: 8/2/24

Python version: 3.12

Environments used: bgmp_py312 (matplotlib 3.9.1, numpy 2.0.0)

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

@K00337:83:HJKJNBBXX:8:1101:1265:1191 #:N:0:1

But the # is replaced with the R number in the file name. So each entry number is linked across the 4 files. It also appears that matching R2 and R3s are the reverse compliment of each other.

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
---
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

Update: scatterplots are ok :tada:

* Quality score cutoff?

Based on the graphs, it looks like a good cutoff for the biological read pairs is 30 because the average is above 30 for all positions, and slightly above 30 for the lowest quality positions. Additionally, Ilumina claims that "Q30 is considered a benchmark for quality in next-generation sequencing". Finally, the quality score cutoff for downstream analysis of biological reads can be generous because aligner algorithms will throw out reads that can't be aligned.

source: https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/quality-scores.html

For the index reads, it is important to consider the hamming distance between indexes. If the quality cutoff is too low, it is possible that reads could be misidentified as belonging to another sample. Since the indexes are presumably well designed, it would take a few low quality reads/mistakes for a sample to be misidentified because the minimun hamming score is likely more than a couple base pairs. Therefore, I think rather than computing a quality threshold by base position, it would be better to calculate an average quality score for the index as a whole. Then if 1 or 2 base pairs are low quality, which wouldn't be enough to cause the sample to be misidentified, the index won't be thrown out from these low reads becasue the average will balance out the outliers. If there are enough low quality reads to make it possible for a sample to be misidentified, then the average would be below the threshold and the read would be kicked out. Based on this, I think a quality score somewhere below 20 and 30 (1 in 100 or 1 in 1000 chance or error) would be reasonable becasue the chance of 1 or more bases being misidentified in an 8 base index with a 1 in 100 error is very low. Therefore, I will select my quality score cutoff for index reads to be 26 (so that reads in the 27+ illumina quality score bin will be included.)

* How many indicies have N base calls?

```
$ for file in /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz; do zcat $file | sed -n '2~4p' | grep -c 'N'; done

3976613
3328051
```
There are 3976613 index 1 barcodes with an N and 3328051 index 2 barcodes with an N.

Questions:
* It is now occuring to be that I made a scatterplot and not a histogram, but I think the scatterplot looks better. I need to ask if this is ok.

yes

* 30 reasonable cutoff

idk... use hamming distance and probabilities

* sum the N indicies or report both files (if both barcodes for a certain sample have N, would this count as 1 or 2)

each file is good

I need to record my answers in ```Answers.md```.

### Part 2 - Write Pseudocode to Demultiplex the Files

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

7/29/24 update: Pseudocode feedback
* change the `while open R1....` to `with open ...`
* add a command line input for q score cutoff
* check if reverse complement of R3 in the index list, not just R3
* statistics are calculated wrong, I need to count the occurance of each possible pair of indexes
* only need to rev comp R3


NOTES:

check other things first then check if complementary (not in list or have N)

R2 and R3 should be rev comps, if not uk (N in index or not in lis of indexes) or hopped
if complementary:
function to see if passes qual score thresh and if it doesnt pass then get added to uk
if it does pass it add R1 to file_R1.fq and R4 to file_R2.fq

if comp, add to header

R2 will match lists of indexes, R3 will be rev complement


---
### 7/27/24
### Assignment the Second
Goal: write feedback to others about their pseudocode

I wrote pseudocode feedback to my assigned people. Here is the feedback that I wrote:

```
Bendycar:
Hi Ben,

Your flowchart is easy to follow and your algorithm makes sense. I like how you included a schematic of how the sample sequences were obtained and the desired output files in the top right corner. One minor thing your algorithm is missing is a way to report the statistics like the number of read-pairs with properly matched and hopped indexes and the number of read-pairs with unknown index(es), but that should be an easy addition. Your functions all seem reasonable and helpful for this algorithm. Adding some more detail to your proposed functions such as what the inputs and outputs would be and a description of what the function does would also make this easier for a reviewer to understand and easier for you to code up later. Overall, this seems like a solid plan and I hope coding this up goes smoothly!

-Jules

graceHach:
Hi Grace,

Your logic is easy to follow and your algorithm makes sense. I like the dictionary system you described to keep track of the statistics for the number of read-pairs for each index pair. I believe your algorithm meets all the requirements outlined in the assignment description. I'm a little confused about the purpose of having a dictionary that holds the status of if a file is open of not - I think you could append each new entry onto the file after each iteration of the loop but maybe I am missing something. One thing that you could do to improve and simplify your algorithm is separate some of the aspects of this code into functions, such as computing the reverse compliment or creating the new header, but this is not entirely necessary. It's also redundant to have a check for both if N is in the index and if the index is in the valid set of 24 because if the index does contain an N it won't be in the set of 24. Overall, this seems like a solid plan and I hope coding this up goes smoothly!

-Jules

claire-j-wells:
Hi Claire,

Your logic is easy to follow and your algorithm makes sense. I like how clear your problem and outputs descriptions are because they give a lot of context to better understand your pseudocode. One minor thing your algorithm is missing is a way to report the statistics like the number of read-pairs with properly matched and hopped indexes and the number of read-pairs with unknown index(es), but that should be an easy addition. Also, what happens to records that make it to step 8 but don't pass the quality threshold? Your functions all seem reasonable and helpful for this algorithm. For the index_seq_to_header(), you might consider just outputting each individual header and then writing that header as you make the output file, rather than outputting new files with new headers. For your write_record_to_file function, you mentioned you weren't sure what the output would be. I think what you currently have (an empty return statement) is appropriate because the body of the function is writing to your file, so there would be no need to return anything after that. Overall, this seems like a solid plan and I hope coding this up goes smoothly!

-Jules
```

7/29/24 update: Grace responded to my feedback and explained why she wrote a dictionary to open and keep track of file handles for writing and close them after. Lesley also mentioned that reopening and appending to files would be slow. I think I am going to implement Grace's strategy into my own code because I was just planning to append before.

---
### 7/30/24
### Still Part 1 :confused:

I need to write my unit test files :rage:

I will write 4 fastq files for each of the 4 illumina output read files that will be in the ```TEST-input_FASTQ``` directory:
```
test_r1.fastq
test_r2.fastq
test_r3.fastq
test_r4.fastq
```
Each file will have 3 reads that meet the following criteria:
* one dual matched read
* one index hopped read
* one unknown index read
* one low quality index read

I will also include the desired output for my files, located in the ```TEST-output_FASTQ``` directory. I have 6 output files for my unit tests:
```
hopped_R1.fastq
hopped_R2.fastq
TCGGATTC_TCGGATTC_R1.fastq
TCGGATTC_TCGGATTC_R2.fastq
unknown_R1.fastq
unknown_R2.fastq
```
When I get my assignment the third code working, I will compare my outputs to these files before running the main file.

--- 
### 8/1/24
### Assignment the Third
Goal: code up my pseudocode to make a working script

Ok time to code. :fearful:

I made a python script called ```demultiplex.py``` in the ```Assigment-the-third``` directory.

Ok so I got the test file to run through, and the output matches the expected unit test output.

I also had my script create a summary of the results called ```demux_stats.txt```.

I then tested the script on a file that was 100,000 records long and it took about .63 seconds to run. Lets calculate the total run time:

363,246,735 records/100000 records * 0.63s / 60s/min / 60min/hr = 0.635

It should take a little over a half hour to run. I am going to make a sbatch script to run this called ```run_demux.sh```. It contains the following:

```
conda activate bgmp_py312
/usr/bin/time -v ./demultiplex.py
```

Time to runnnnnnnnnnnnnnnnnnn

```
sbatch run_demux.sh
Submitted batch job 7827370
```

---
### 8/3/24
### Assignment the Third

It worked!!!!!

JobId: 7827370
Runtime: 1 hour 15 mins
% CPU: 86%

The output is a file called ```demux_stats.txt``` which has all the statistics. And all the output files are in the ```outputs``` directory inside of ```Assignment-the-third```

But I realized I didn't calculate all the statistics. I need to add the percentage of reads from each sample and a figure. I will go add that.

Now, the script will output a file called ```demux_stats.txt``` with the statistics, and a graph called ```demux_stats.png``` that shows the percentage of the total reads that each sample contributes.

I am going to run the big boy files again through my script now that I have refined the script and added to my statistics output.

```
sbatch run_demux.sh
Submitted batch job 7892950
```
JobId: 7892950
Runtime: 1 hour 25 mins
% CPU: 86%

IT WORKED!!!

My graph looks good. The stats are complete now. Time to submit.


Part 3 notes:
- include bar chart of reads/sample
- add a space between the first part of header and the barcodes
- open and close the file once
- use itertools permutations to make all permutations of file outputs maybe?
- only need named files for each matched pair



