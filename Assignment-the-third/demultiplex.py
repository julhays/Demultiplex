#!/usr/bin/env python

# Author: Jules Hays

# Demultiplex The Third: Algorithm to Demultiplex the Samples

from itertools import islice
import gzip

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



#make and open all the write files
output_files: dict = {}   #keys are barcode, values are list of 2 file handles for R1 and R2
#output_names = ['hopped_R1', 'hopped_R2', 'unknown_R1', 'unknown_r2']
for barcode in indexes:
    output_name1 = f'{barcode}_{barcode}_R1.fastq'
    output_name2 = f'{barcode}_{barcode}_R2.fastq'
    output_files[barcode] = [open(output_name1, 'wt'), open(output_name2, 'wt')]

output_files['hopped'] = [open('hopped_R1.fastq', 'wt'), open('hopped_R2.fastq', 'wt')]
output_files['unknown'] = [open('unknown_R1.fastq', 'wt'), open('unknown_R2.fastq', 'wt')]



#open all the read files
#with gzip.open
with open(r1_path, "rt") as r1, open(r2_path, "rt") as r2, open(r3_path, "rt") as r3, open(r4_path, "rt") as r4:
    while True:
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

#CLOSE FILES
for key in output_files:
    output_files[key][0].close()
    output_files[key][1].close()


exit()

