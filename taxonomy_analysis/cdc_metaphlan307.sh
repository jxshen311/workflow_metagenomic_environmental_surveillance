#!/bin/bash
#SBATCH --job-name="cdc_metaphlan307" 
#SBATCH -A p30892             
#SBATCH -p short                # Queue/partition
#SBATCH -t 4:00:00             # Walltime/duration of the job
#SBATCH -N 1                    # Number of Nodes
#SBATCH --mem=0    # --mem=0 means you take the whole node  
#SBATCH --ntasks-per-node=24     # Number of Cores (Processors/CPU)
#SBATCH --mail-user=jiaxianshen2022@u.northwestern.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/projects/p30892/cdc/log/"%x_%j.log"  
#SBATCH --error=/projects/p30892/cdc/log/"%x_%j.err"
     

module purge
module load singularity

knead_seq=`ls /projects/p30892/cdc/metaxa/input | grep _paired_1 | cut -d "_" -f1,2,3,4,5,6,7`
for input in $knead_seq
do 
	output=`echo $input | cut -d "_" -f1`
	singularity exec -B /projects/p30892/cdc workflows_3.0.0.a.6.metaphlanv3.0.7.sif metaphlan metaxa/input/"${input}_1.fastq.gz",metaxa/input/"${input}_2.fastq.gz" --bowtie2out metaphlan307_inter_bowtie2/"${output}.bowtie2.bz2" --input_type fastq --nproc 24 --unknown_estimation -o metaphlan307/"${output}.txt"
done

