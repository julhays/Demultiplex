#!/usr/bin/env python

# Author: Jules Hays

# Demultiplex The Third: Algorithm to Demultiplex the Samples

from itertools import islice
import gzip
import bioinfo

# ADD ARGPARSE

index_file: str = '/projects/bgmp/shared/2017_sequencing/indexes.txt'

r1_path: str = '/projects/bgmp/jkhay/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test_R1.fastq'
r2_path: str = '/projects/bgmp/jkhay/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test_R2.fastq'
r3_path: str = '/projects/bgmp/jkhay/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test_R3.fastq'
r4_path: str = '/projects/bgmp/jkhay/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test_R4.fastq'

#r1_path: str = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz'
#r2_path: str = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz'
#r3_path: str = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz'
#r4_path: str = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz'


# make list of acceptable indexes from the indexes file
index_list: list = []

with open(index_file, 'rt') as fh:
    next(fh)  #skip the header line in the index file
    for line in fh:
        line = line.strip()
        index = line.split('\t')[4]   #get the 5th column, which is the index value
        index_list.append(index)

#turn it into a set for faster lookups
indexes: set = set(index_list)

#necessary functions
def reverse_complement(seq: str) -> str:
    '''Takes in the a DNA/RNA sequence and returns the reverse complement of that sequence 
    written 5' -> 3' '''
    complements = {'A': 'T', 'T': 'A', 'G': 'C', 'C':'G', 'N':'N'}   #key is nucleotide, value is reverse complement
    seq_list= list(seq)  #cast sequence to a list of chars
    comp_seq = [complements[base] for base in seq_list]  #make a new list of the rev comp of each nt in the sequence
    rev_comp_seq = ''.join(comp_seq[::-1])  #reverse the list and cast it back to a string
    return rev_comp_seq

def good_qual(R2_qual: str, R3_qual: str, threshold) -> bool:
    '''Takes in the barcode quality score lines in R2 and R3 and returns False if 
    the average quality score of either does not meet or exceed the cutoff'''
    #caculate the average quality score for index 1 and 2
    r2_avg = bioinfo.qual_score(R2_qual)
    r3_avg = bioinfo.qual_score(R3_qual)

    #if either as less than the threshold, returns false, else true
    return (r2_avg < threshold) or (r3_avg < threshold)



def barcodes_to_header(R1_head: str, R4_head: str, R2_seq: str, R3_seq: str) -> tuple:
    '''Takes in the R1 and R4 header and appends the barcode sequences 
    listed in R2 and the reverse complement of R3 to the end of the file header to make it 'header BARCODE1-rcBARCODE2' '''
    new_header_1 = f'{R1_head} {R2_seq}-{R3_seq}'
    new_header_2 = f'{R4_head} {R2_seq}-{R3_seq}'
    return new_header_1, new_header_2

#make and open all the write files
output_files: dict = {}   #keys are barcode, values are list of 2 file handles for R1 and R2
#output_names = ['hopped_R1', 'hopped_R2', 'unknown_R1', 'unknown_r2']
for barcode in indexes:
    output_name1 = f'outputs/{barcode}_{barcode}_R1.fastq'
    output_name2 = f'outputs/{barcode}_{barcode}_R2.fastq'
    output_files[barcode] = [open(output_name1, 'wt'), open(output_name2, 'wt')]

output_files['hopped'] = [open('outputs/hopped_R1.fastq', 'wt'), open('outputs/hopped_R2.fastq', 'wt')]
output_files['unknown'] = [open('outputs/unknown_R1.fastq', 'wt'), open('outputs/unknown_R2.fastq', 'wt')]


#open all the read files
#with gzip.open
with open(r1_path, "rt") as r1, open(r2_path, "rt") as r2, open(r3_path, "rt") as r3, open(r4_path, "rt") as r4:
    while True:
        #isolate a record from each file
        r1_lines = []
        r2_lines = []
        r3_lines = []
        r4_lines = []
        for i in range(4):
            r1_lines.append(r1.readline().strip('\n'))
            r2_lines.append(r2.readline().strip('\n'))
            r3_lines.append(r3.readline().strip('\n'))
            r4_lines.append(r4.readline().strip('\n'))
        if r1_lines[0] == '':
            break

        #main body of the code
        #extract the indexes and reverse complement R3 index
        r2_index = r2_lines[1]
        r3_index = reverse_complement(r3_lines[1])

        #add indexes to headers
        r1_lines[0], r4_lines[0] = barcodes_to_header(r1_lines[0], r4_lines[0], r2_index, r3_index)

        #check if the indexes exist in the set of valid indexes
        if (r2_index not in indexes) or (r3_index not in indexes):
            r1_out = f'{r1_lines[0]}\n{r1_lines[1]}\n{r1_lines[2]}\n{r1_lines[3]}\n'
            r2_out = f'{r4_lines[0]}\n{r4_lines[1]}\n{r4_lines[2]}\n{r4_lines[3]}\n'
            output_files['unknown'][0].write(r1_out)
            output_files['unknown'][1].write(r2_out)

        

#CLOSE FILES
for key in output_files:
    output_files[key][0].close()
    output_files[key][1].close()


exit()

