#!/bin/bash
#SBATCH --partition=week
#SBATCH --job-name=modulzr
#SBATCH -c 5
#SBATCH --mem=20G
#SBATCH --time=3-05:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


# MAJIQ workflow v3: in v2 grouped all the biological and technical replicates by neuron. In v3, group the technical replicates (if any) by biological replicate.
# MAJIQ v4: use the merged bams (no technical replicate), group the biological replicates by neuron
# v5: also perform tests for every neuron pair
# v6: add modulizer
# v7: prepare ref with gene names and add heterogen
# Extract only the Modulizer step


set -ue

echo "Starting $(date)"

# ----------------------------------------------------------------------------
# -----------------------    Parameters definitions   ------------------------
# ---------------------------------------------------------------------------- 

# Locations and parameters

majiqdir='/gpfs/ycga/project/ysm/hammarlund/aw853/majiq/outs/2024-03-04_outs'
outdir='/gpfs/ycga/project/ysm/hammarlund/aw853/majiq/outs/2024-03-22_outs'


WSversion='WS289'








# the output subdirectory
mkdir -p $outdir/voila_modulizer


source /home/aw853/bin/majiq_v2-3/env/bin/activate
module load Python/3.8.6-GCCcore-10.2.0 HTSlib/1.12-GCCcore-10.2.0









echo "##################################################################"
echo "##################           Modulize           ##################"
echo "##################################################################"
echo
echo


# modulize ----

echo 
echo "Running Voila Modulizer on $(date)"
echo


voila modulize \
  --overwrite \
  --nproc $SLURM_CPUS_PER_TASK \
  --logger $outdir/logs/modulize.log \
  --keep-constitutive --keep-no-lsvs-modules --keep-no-lsvs-junctions \
  -d $outdir/voila_modulizer \
  $majiqdir/psi $majiqdir/build/splicegraph.sql

echo "------------ End modulizer ------------"
echo







echo
echo "All done $(date)"
