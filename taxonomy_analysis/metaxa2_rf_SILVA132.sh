#!/bin/bash
#SBATCH --job-name="cdc_metaxa2_rf_SILVA132"
#SBATCH -A p30892             
#SBATCH -p short                # Queue/partition
#SBATCH -t 04:00:00             # Walltime/duration of the job
#SBATCH -N 1                    # Number of Nodes
#SBATCH --mem=0    # --mem=0 means you take the whole node  
#SBATCH --ntasks-per-node=24     # Number of Cores (Processors/CPU)
#SBATCH --mail-user=jiaxianshen2022@u.northwestern.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/projects/p30892/cdc/log/metaxa2_rf_SILVA132.log  # need to be a file
#SBATCH --error=/projects/p30892/cdc/log/metaxa2_rf_SILVA132.err  # need to be a file


module purge

module load hmmer/3.1b2 blast/2.7.1 mafft/7.407 python/anaconda3 metaxa2/2.2

sample_list=`ls output_SILVA132_SSU | grep .taxonomy.txt | cut -d "." -f1`
for input in $sample_list
do
	metaxa2_rf -i output_SILVA132_SSU/"${input}.taxonomy.txt" -o output_diversity_tools/output_rf_SILVA132/"${input}_rf" --unknown T --cpu 24 --multi_thread T
done