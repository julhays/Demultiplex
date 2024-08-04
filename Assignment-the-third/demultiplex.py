#!/usr/bin/env python

# Author: Jules Hays

# Demultiplex The Third: Algorithm to Demultiplex the Samples

import argparse
import itertools
import gzip
import bioinfo
import matplotlib.pyplot as plt

# ADD ARGPARSE
#define arg inputs, defaults are for the actual files we are demultiplexing for this assignment
def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Establish criteria for k-merizing fastq data")
    parser.add_argument("-r1", "--R1", help="Specify the R1 file", type=str, default ='/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz')
    parser.add_argument("-r2", "--R2", help="Specify the R2 file", type=str, default ='/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz')
    parser.add_argument("-r3", "--R3", help="Specify the R3 file", type=str, default ='/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz')
    parser.add_argument("-r4", "--R4", help="Specify the R4 file", type=str, default ='/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz')
    parser.add_argument("-i", "--index", help="Specify the index file", type=str, default ='/projects/bgmp/shared/2017_sequencing/indexes.txt')
    return parser.parse_args()

#call get_args to create args object
args = get_args()

#set path variables and assign them to the user inputted values at the function call
r1_path: str = args.R1
r2_path: str = args.R2
r3_path: str = args.R3
r4_path: str = args.R4
index_file: str = args.index

#hard code paths for test files
# index_file: str = '/projects/bgmp/shared/2017_sequencing/indexes.txt'

# r1_path: str = '/projects/bgmp/jkhay/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test_R1.fastq'
# r2_path: str = '/projects/bgmp/jkhay/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test_R2.fastq'
# r3_path: str = '/projects/bgmp/jkhay/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test_R3.fastq'
# r4_path: str = '/projects/bgmp/jkhay/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test_R4.fastq'

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

#set up for statistics
# make a dictionary to store all permutations of 2 indexes from the index list
# start at 0, as reads are assigned to files it will update the number of occurances 
matched_pairs = {}  #key = tuple of matched index pair, value = number of occurances in the files
hopped_pairs = {}   #key = tuple of nonmatched index pair, value = number of occurances in the files

barcode_permutations = list(itertools.permutations(index_list, 2))  #make a tuple for every possible combination (unmatched) of the indexes
#add tuple to dictionary and initialize to 0
for pair in barcode_permutations:
    hopped_pairs[pair] = 0
#add tuple of each possible matched pair to matched dictionary
for index in indexes:
    matched_pairs[(index, index)] = 0

#integer variable to store the number of matched, hopped, and unknown index pairs
matched = 0
hopped = 0
unknown = 0

#set a counter to keep track of number of reads
read_count = 0

#necessary functions
def reverse_complement(seq: str) -> str:
    '''Takes in the a DNA/RNA sequence and returns the reverse complement of that sequence 
    written 5' -> 3' '''
    complements = {'A': 'T', 'T': 'A', 'G': 'C', 'C':'G', 'N':'N'}   #key is nucleotide, value is reverse complement
    rev_comp_seq = ""  #make empty string to store rev comp
    for nt in seq:
        rev_comp_seq = complements[nt] + rev_comp_seq   #find complement than add to front of string
    return rev_comp_seq

def good_qual(R2_qual: str, R3_qual: str, threshold: int) -> bool:
    '''Takes in the barcode quality score lines in R2 and R3 and returns False if 
    the average quality score of either does not meet or exceed the cutoff'''
    #caculate the average quality score for index 1 and 2
    r2_avg = bioinfo.qual_score(R2_qual)
    r3_avg = bioinfo.qual_score(R3_qual)

    #if either as less than the threshold, returns false, else true
    return (r2_avg >= threshold) and (r3_avg >= threshold)

def barcodes_to_header(R1_head: str, R4_head: str, R2_seq: str, R3_seq: str) -> tuple:
    '''Takes in the R1 and R4 header and appends the barcode sequences 
    listed in R2 and the reverse complement of R3 to the end of the file header to make it 'header BARCODE1-rcBARCODE2' '''
    new_header_1 = f'{R1_head} {R2_seq}-{R3_seq}' #appends sequence info to header
    new_header_2 = f'{R4_head} {R2_seq}-{R3_seq}'
    return new_header_1, new_header_2

#make and open all the write files
output_files: dict = {}   #keys are barcode, values are list of 2 file handles for R1 and R2
#output_names = ['hopped_R1', 'hopped_R2', 'unknown_R1', 'unknown_r2']
for barcode in indexes:
    #define output file names
    output_name1 = f'outputs/{barcode}_{barcode}_R1.fastq'
    output_name2 = f'outputs/{barcode}_{barcode}_R2.fastq'
    #add output file to dictionary
    output_files[barcode] = [open(output_name1, 'wt'), open(output_name2, 'wt')] #keys are barcodes, values are a list of the R1 and R2 file handles for each
#create an entry for hopped and unknown files
output_files['hopped'] = [open('outputs/hopped_R1.fastq', 'wt'), open('outputs/hopped_R2.fastq', 'wt')]
output_files['unknown'] = [open('outputs/unknown_R1.fastq', 'wt'), open('outputs/unknown_R2.fastq', 'wt')]

#open all the read files
#with open(r1_path, "rt") as r1, open(r2_path, "rt") as r2, open(r3_path, "rt") as r3, open(r4_path, "rt") as r4:
with gzip.open(r1_path, "rt") as r1, gzip.open(r2_path, "rt") as r2, gzip.open(r3_path, "rt") as r3, gzip.open(r4_path, "rt") as r4:
    while True:
        #create an empty list to store the record from each file in
        r1_lines = []
        r2_lines = []
        r3_lines = []
        r4_lines = []
        #add 4 lines to each list
        for i in range(4):
            r1_lines.append(r1.readline().strip('\n'))
            r2_lines.append(r2.readline().strip('\n'))
            r3_lines.append(r3.readline().strip('\n'))
            r4_lines.append(r4.readline().strip('\n'))
        #end the while true loop if there are no more lines to read in
        if r1_lines[0] == '':
            break
        
        #increment record count since there was a record added
        read_count += 1

        #main body of the code
        #extract the indexes and reverse complement R3 index
        r2_index = r2_lines[1]
        r3_index = reverse_complement(r3_lines[1])

        #add indexes to headers
        r1_lines[0], r4_lines[0] = barcodes_to_header(r1_lines[0], r4_lines[0], r2_index, r3_index)

        #define what will be written to a file for this record
        r1_out = f'{r1_lines[0]}\n{r1_lines[1]}\n{r1_lines[2]}\n{r1_lines[3]}\n'
        r2_out = f'{r4_lines[0]}\n{r4_lines[1]}\n{r4_lines[2]}\n{r4_lines[3]}\n'

        #check if the indexes exist in the set of valid indexes
        if (r2_index not in indexes) or (r3_index not in indexes):
            output_files['unknown'][0].write(r1_out) #write to the file
            output_files['unknown'][1].write(r2_out)
            unknown += 1
            continue #go back to the top of the while loop

        #check if indexes are good quality
        if not good_qual(r2_lines[3], r3_lines[3], 26):
            output_files['unknown'][0].write(r1_out)
            output_files['unknown'][1].write(r2_out)
            unknown += 1
            continue

        #check if indexes match
        if r2_index == r3_index:
            output_files[r2_index][0].write(r1_out)
            output_files[r2_index][1].write(r2_out)
            matched_pairs[(r2_index, r3_index)] += 1
            matched += 1

        elif r2_index != r3_index:
            output_files['hopped'][0].write(r1_out)
            output_files['hopped'][1].write(r2_out)
            hopped_pairs[(r2_index, r3_index)] += 1
            hopped += 1

        else:
            print('you messed up')  #should not incur this condition because the above checks should capture everything
        

#CLOSE FILES
for key in output_files: #loops through the dictionary to close each file
    output_files[key][0].close()
    output_files[key][1].close()

#output statistics into a stats file
with open('demux_stats.txt', "wt") as stats:
    stats.write(f'Number of total reads: {read_count}\n')
    stats.write(f'Number of matched reads: {matched}\n')
    stats.write(f'Percent of matched reads: {100*matched/read_count}%\n')
    stats.write(f'Number of hopped reads: {hopped}\n')
    stats.write(f'Percent of hopped reads (amount of index swapping): {100*hopped/read_count}%\n')
    stats.write(f'Number of unknown reads: {unknown}\n')
    stats.write(f'Percent of unknown reads: {100*unknown/read_count}%\n\n')
    stats.write(f'Breakdown of each read pair and number of times it occurs:\n\n')
    stats.write('Matched index pair\tNumber of occurances\tPercentage\n')

    for key in matched_pairs:
        percent = 100*matched_pairs[key]/read_count
        stats.write(f'{key}\t{matched_pairs[key]}\t{percent}\n')

    stats.write('\nHopped index pair\tNumber of occurances\tPercentage\n')

    for key in hopped_pairs:
        percent = 100*hopped_pairs[key]/read_count
        stats.write(f'{key}\t{hopped_pairs[key]}\t{percent}\n')
    
# make a graph of percentage of reads in each sample
sample = list(indexes)
percent = [100*matched_pairs[(index, index)]/read_count for index in sample]

plt.barh(sample, percent)
plt.xlabel("Percentage of Sample in the Reads")
plt.ylabel("Macthed Sample Index Sequence")
plt.title("Percentage of Each Sample in Demultiplexed Reads")
plt.tight_layout()  #the graph was cutting off the y axis label, so this makes it not do that
plt.savefig('demux_stats.png')
plt.close()



#UNIT TEST CHECK
#compare the outputs with my expected outputs for my unit tests
#unit tests passed!
# unit_test_pred = ['/projects/bgmp/jkhay/bioinfo/Bi622/Demultiplex/TEST-output_FASTQ/TCGGATTC_TCGGATTC_R1.fastq',
#                   '/projects/bgmp/jkhay/bioinfo/Bi622/Demultiplex/TEST-output_FASTQ/TCGGATTC_TCGGATTC_R2.fastq',
#                   '/projects/bgmp/jkhay/bioinfo/Bi622/Demultiplex/TEST-output_FASTQ/hopped_R1.fastq',
#                   '/projects/bgmp/jkhay/bioinfo/Bi622/Demultiplex/TEST-output_FASTQ/hopped_R2.fastq',
#                   '/projects/bgmp/jkhay/bioinfo/Bi622/Demultiplex/TEST-output_FASTQ/unknown_R1.fastq',
#                   '/projects/bgmp/jkhay/bioinfo/Bi622/Demultiplex/TEST-output_FASTQ/unknown_R2.fastq']

# unit_test_actual = ['/projects/bgmp/jkhay/bioinfo/Bi622/Demultiplex/Assignment-the-third/outputs/TCGGATTC_TCGGATTC_R1.fastq',
#                     '/projects/bgmp/jkhay/bioinfo/Bi622/Demultiplex/Assignment-the-third/outputs/TCGGATTC_TCGGATTC_R2.fastq',
#                     '/projects/bgmp/jkhay/bioinfo/Bi622/Demultiplex/Assignment-the-third/outputs/hopped_R1.fastq',
#                     '/projects/bgmp/jkhay/bioinfo/Bi622/Demultiplex/Assignment-the-third/outputs/hopped_R2.fastq',
#                     '/projects/bgmp/jkhay/bioinfo/Bi622/Demultiplex/Assignment-the-third/outputs/unknown_R1.fastq',
#                     '/projects/bgmp/jkhay/bioinfo/Bi622/Demultiplex/Assignment-the-third/outputs/unknown_R2.fastq']

# for i in range(len(unit_test_actual)):
#     with open(unit_test_pred[i], 'rt') as pred, open(unit_test_actual[i], 'rt') as actual:
#         while True:
#             line_a = pred.readline().strip('\n')
#             line_b = actual.readline().strip('\n')
#             print(line_a)
#             print(line_b)
#             #check 
#             if line_a != line_b:
#                 print('you messed up')
#                 break
#             if line_a == "" and line_b == "":
#                 print('you did it!')
#                 break


exit()

