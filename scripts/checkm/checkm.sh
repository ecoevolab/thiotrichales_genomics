#!/bin/bash
#$ -N checkmbeggiatoas1
#$ -o $JOB_NAME.log 
#$ -e $JOB_NAME.error 
#$ -cwd
#$ -S /bin/bash
#$ -l h_rss=60G
#$ -l h_rt=336:00:00  
#$ -pe openmp 10

date
echo "===== Beginning pipeline ====="

eval "$(conda shell.bash hook)"
conda activate checkm

checkm lineage_wf \
  --tab_table \
  -t 10 \
  -x fna \
  /mnt/atgc-d3/sur/shared_data/beggiatoas/fna_links \
  /mnt/atgc-d3/sur/users/mreyesr/exp/beggiatoas_outputs/outputs_checkm2

echo "===== Pipeline done ====="
date

