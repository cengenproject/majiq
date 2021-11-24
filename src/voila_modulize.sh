#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=voila_mod
#SBATCH -c 15
#SBATCH --mem=20G
#SBATCH --time=2-05:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


source /home/aw853/bin/majiq_v2-3/env/bin/activate
module load Python/3.8.6-GCCcore-10.2.0 HTSlib/1.12-GCCcore-10.2.0




outdir='/gpfs/ycga/project/ysm/hammarlund/aw853/counts/majiq/2021-11-10_outs'

mkdir -p $outdir/voila_modulizer

voila modulize \
  --overwrite \
  --nproc $SLURM_CPUS_PER_TASK \
  --logger $outdir/logs \
  --keep-constitutive --keep-no-lsvs-modules --keep-no-lsvs-junctions \
  -d $outdir/voila_modulizer
  $outdir/build/splicegraph.sql $outdir/deltapsi

