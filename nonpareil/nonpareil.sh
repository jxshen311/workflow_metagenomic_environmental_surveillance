#!/bin/bash
#SBATCH --job-name="cdc_nonpareil"
#SBATCH -A p30892             
#SBATCH -p short                # Queue/partition
#SBATCH -t 4:00:00             # Walltime/duration of the job
#SBATCH -N 1                    # Number of Nodes
#SBATCH --mem-per-cpu=48G 
#SBATCH --ntasks-per-node=16     # Number of Cores (Processors/CPU)
#SBATCH --mail-user=jiaxianshen2022@u.northwestern.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/projects/p30892/cdc/log/nonpareil.log  # need to be a file
#SBATCH --error=/projects/p30892/cdc/log/nonpareil.err  # need to be a file

module purge
 
module load anaconda3
source activate nonpareil

sample_list=`ls /projects/p30892/cdc/metaxa/input | grep _paired_1.fastq`
for input in $sample_list
do
  output=`echo $input | cut -d "_" -f1`
  nonpareil -s /projects/p30892/cdc/metaxa/input/$input -T kmer -f fastq -b output/$output -t 16
done