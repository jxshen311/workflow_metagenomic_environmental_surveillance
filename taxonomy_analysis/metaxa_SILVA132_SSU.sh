#!/bin/bash
#SBATCH --job-name="cdc_metaxa_SILVA132_SSU"
#SBATCH -A p30892             
#SBATCH -p normal                # Queue/partition
#SBATCH -t 8:00:00             # Walltime/duration of the job
#SBATCH -N 1                    # Number of Nodes
#SBATCH --mem=0    # --mem=0 means you take the whole node  
#SBATCH --ntasks-per-node=24     # Number of Cores (Processors/CPU)
#SBATCH --mail-user=jiaxianshen2022@u.northwestern.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/projects/p30892/cdc/log/metaxa_SILVA132_SSU.log  # need to be a file
#SBATCH --error=/projects/p30892/cdc/log/metaxa_SILVA132_SSU.err  # need to be a file


# unload any modules that carried over from your command line session
module purge

# load modules you need to use
module load hmmer/3.1b2 blast/2.7.1 mafft/7.407 python/anaconda3 metaxa2/2.2

# use loop to run
sample_list=`ls input | grep kneaddata_paired_1.fastq | cut -d "_" -f1,2`
for input in $sample_list
do
	metaxa2 -1 input/"${input}_L001_R1_001_kneaddata_paired_1.fastq" -2 input/"${input}_L001_R1_001_kneaddata_paired_2.fastq" --mode metagenome -g SILVA132_SSU  -p /home/jsz1618/bin/metaxa2_db/SILVA132_SSU/HMMs/ -o output_SILVA132_SSU/${input%_*} --cpu 24 --multi_thread T --plus T
done