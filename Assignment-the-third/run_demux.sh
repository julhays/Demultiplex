#!/bin/bash

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH -c 8                 #optional: number of cpus, default is 1
#SBATCH --mem=100GB                        #optional: amount of memory, default is 4GB per cpu
#SBATCH --mail-user=jkhay@uoregon.edu     #optional: if you'd like email
#SBATCH --mail-type=ALL                   #optional: must set email first, what type of email you want
#SBATCH --job-name=demux           #optional: job name
#SBATCH --output=demux_%j.out       #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=demux_%j.err        #optional: file to store stderr from job, %j adds the assigned jobID


conda activate bgmp_py312

/usr/bin/time -v ./demultiplex.py

exit