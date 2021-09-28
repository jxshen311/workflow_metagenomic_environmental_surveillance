#!/bin/bash
#SBATCH --job-name="cdc_metaphlan" 
#SBATCH -A p30892             
#SBATCH -p normal                # Queue/partition
#SBATCH -t 8:00:00             # Walltime/duration of the job
#SBATCH -N 1                    # Number of Nodes
#SBATCH --mem=0    # --mem=0 means you take the whole node  
#SBATCH --ntasks-per-node=24     # Number of Cores (Processors/CPU)
#SBATCH --mail-user=jiaxianshen2022@u.northwestern.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/projects/p30892/cdc/log/metaphlan.log  # need to be a file
#SBATCH --error=/projects/p30892/cdc/log/metaphlan.err  # need to be a file 
     

# unload any modules that carried over from your command line session
module purge

# load modules you need to use
module load singularity


# use loop to run metaphlan on multiple samples
knead_seq=`ls /projects/p30892/cdc/metaphlan/input`
for input in $knead_seq
do 
	output=`echo $input | cut -d "_" -f1`.txt
	singularity exec -B /home/jsz1618/ -B /projects/p30892/ /home/jsz1618/biobakery-diamondv0822.simg metaphlan2.py /projects/p30892/cdc/metaphlan/input/$input --input_type fastq --nproc 24 > /projects/p30892/cdc/metaphlan/output/$output
done
