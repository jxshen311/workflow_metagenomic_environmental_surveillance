#!/bin/bash

#SBATCH --job-name="cdc_kneaddata_1" 
#SBATCH -A p30892             
#SBATCH -p short               # Queue/partition
#SBATCH -t 4:00:00             # Walltime/duration of the job
#SBATCH -N 1                    # Number of Nodes
#SBATCH --mem=0    # --mem=0 means you take the whole node  
#SBATCH --ntasks-per-node=24     # Number of Cores (Processors/CPU)
#SBATCH --mail-user=jiaxianshen2022@u.northwestern.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/projects/p30892/cdc/log/kneaddata_1.log  # need to be a file
#SBATCH --error=/projects/p30892/cdc/log/kneaddata_1.err  # need to be a file

module purge

module load singularity

while read P1
do
	singularity exec -B /home/jsz1618/ -B /projects/p30892/ /home/jsz1618/biobakery-diamondv0822.simg kneaddata -t 24 --input /projects/p30892/cdc/raw_data/"${P1}_L001_R1_001.fastq.gz" --input /projects/p30892/cdc/raw_data/"${P1}_L001_R2_001.fastq.gz" -db /home/jsz1618/kneaddata/refer_db/Homo_sapiens -db /projects/p30892/cdc/kneaddata/refer_db/negative_general_db --output /projects/p30892/cdc/kneaddata/kneadout/batch1/"${P1%S*}out"

done < /projects/p30892/cdc/kneaddata/list1.txt

# compress knead_out fastq files
gzip /projects/p30892/cdc/kneaddata/kneadout/batch1/*/*.fastq
