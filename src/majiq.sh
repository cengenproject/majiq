#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=majiq2
#SBATCH -c 15
#SBATCH --mem=15G
#SBATCH --time=2-05:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


# MAJIQ workflow v3: in v2 grouped all the biological and technical replicates by neuron. In v3, group the technical replicates (if any) by biological replicate.

# ----------------------------------------------------------------------------
# -----------------------    Parameters definitions   ------------------------
# ---------------------------------------------------------------------------- 

# Locations and parameters
bam_dir='/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/bsn3'
WSversion='WS277'
references_dir='/gpfs/ycga/project/ysm/hammarlund/aw853/references'

outdir='/gpfs/ycga/project/ysm/hammarlund/aw853/counts/2021-05-06_majiq'


# ----------------------------------------------------------------------------
# ---------      No modification required below this line     ----------------
# ----------------------------------------------------------------------------


# Define references
ref_gff_name='c_elegans.PRJNA13758.'$WSversion'.canonical_geneset_gtf2gff3.gff3'

ref_gff=$references_dir'/'$WSversion'/'$ref_gff_name

cfg_file=$(date +%F)"_majiq.cfg"

# the output subdirectories
mkdir -p $outdir/build
mkdir -p $outdir/psi
mkdir -p $outdir/logs


# --------    End definitions, start computations    --------

set -e

# List the samples to process ----
samplelist=($(ls $bam_dir/*bam | xargs basename -a -s .bam))
declare -A sampByNeur #for the config file to the BUILD command
declare -A sampByNeur_path # for the PSI command

for smp in "${samplelist[@]}"
do
	if [[ $smp =~ ^([A-Z1-4ef]{2,4}r[0-9]{1,4})[t12]*$ ]]
	then
		if [[ -z ${sampByNeur[${BASH_REMATCH[1]}]} ]]
		then
			#first time that neuron encountered
			sampByNeur[${BASH_REMATCH[1]}]="${BASH_REMATCH[1]}=$smp"
			sampByNeur_path[${BASH_REMATCH[1]}]=$outdir/build/$smp.majiq
		else
			#neuron already initialized
			sampByNeur[${BASH_REMATCH[1]}]+=","$smp
			sampByNeur_path[${BASH_REMATCH[1]}]+=" $outdir/build/$smp.majiq"
		fi
	else
		#no match
		echo "Error matching regex for sample $smp"
		exit 1
	fi
done


# remove the Ref samples that we don't want here
unset sampByNeur[Ref]
unset sampByNeur_path[Ref]


# Create config file ----
echo "[info]
bamdirs=$bam_dir
genome=$WSversion
strandness=forward
[experiments]" > $cfg_file
printf "%s\n" ${sampByNeur[@]} >> $cfg_file

echo "Config file generated and saved as $cfg_file"

# build ----

echo
echo "Running MAJIQ build on $(date) with gffread and autocreated conf file"
echo

source /home/aw853/bin/majiq/env/bin/activate


majiq build --nproc $SLURM_CPUS_PER_TASK \
			-o $outdir/build \
			--conf $cfg_file \
			--min-experiments 2 \
			$ref_gff

echo "Finished the build step."


# psi ----

echo
echo "Running MAJIQ PSI on $(date) with conf $cfg_file"
echo


for neur in ${!sampByNeur_path[@]}
do
	echo "---- PSI for $neur ----"
	majiq psi --output $outdir/psi \
		  --name $neur \
		  --nproc $SLURM_CPUS_PER_TASK \
		  --output-type all \
		  --logger $outdir/logs \
		  ${sampByNeur_path[$neur]}
done

echo "Done"
echo


# export relevant data ----
mv $outdir/build/splicegraph.sql $outdir/splicegraph.sql
mv $outdir/build/majiq.log $outdir/logs/majiq_build.log
mv $outdir/psi/psi_majiq.log $outdir/logs/majiq_psi.log

# to download on local computer for Voila and R analyses, then can be deleted from cluster
tar -czf $outdir/$(date +%y%m%d)_mjq_exprt.tar.gz $outdir/psi/ $outdir/splicegraph.sql

