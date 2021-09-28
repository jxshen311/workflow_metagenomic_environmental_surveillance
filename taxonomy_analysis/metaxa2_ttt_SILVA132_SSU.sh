#!/bin/bash
module purge

module load hmmer/3.1b2 blast/2.7.1 mafft/7.407 python/anaconda3 metaxa2/2.2

# metaxa_ttt (to count the number of species, genera et)
sample_list=`ls output_SILVA132_SSU | grep .taxonomy.txt | cut -d "." -f1`

for input in $sample_list
do
  metaxa2_ttt -i output_SILVA132_SSU/"${input}.taxonomy.txt" -o output_SILVA132_SSU_ttt/${input} --cpu 24 --multi_thread T --unknown T
done
