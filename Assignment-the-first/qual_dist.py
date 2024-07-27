#!/usr/bin/env python

# Author: Jules Hays

# Demultiplex The First Part 1: Quality Score Distribution Per-nucleotide

import argparse
import matplotlib.pyplot as plt
import bioinfo
import numpy as np
import re
import gzip

#define arg inputs, defaults are for the lane1_NoIndex_L001_R1_003.fastq file
def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Establish criteria for k-merizing fastq data")
    parser.add_argument("-f", "--file_path", help="Specify the path of the input fastq file", type=str, required=True)
    parser.add_argument("-n", "--read_len", help="Specify the read length of each read in the file", type=int, default=101)
    return parser.parse_args()

#call get_args to create args object
args = get_args()

#set global variables and assign them to the user inputted values at the function call
file: str = args.file_path
n: int = args.read_len

#initializes array of 0.0 the length of the read length (n) to store the mean quality score per position
means: np.array = np.zeros(n)

#loop through input fastq file
with gzip.open(file, "rt") as fh:
    counter = 0
    for line in fh:
        counter+=1
        line = line.strip('\n')
        if counter%4 == 0:              #pulls only the quality score line
            #calculate the q_score for each value and add it to respective spot in qual_sums
            for base_index, value in enumerate(line):
                #convert value to score
                score = bioinfo.convert_phred(value)
                means[base_index] += score
            if counter%1000000 == 0:
                print(counter, '/1452986940')

#calculate the number of reads from the line counter
num_records: float = counter/4

#convert the sum of q scores to he mean qscore
for i in range(len(means)):
    means[i] = means[i]/num_records

#grab the R number from the file name
file_r = re.findall(r'.+_R([1-4])', file)[0]

#plot the mean quality score distribution
plt.plot(means, '-')
plt.xlabel("Nucleotide Position")
plt.ylabel("Average Quality Score")
plt.title(f'Average Quality Scores for Nucleotide Position in R{file_r}')
plt.savefig(f'R{file_r}_dist.png')
plt.close()

#plt.xlim(0,150)

