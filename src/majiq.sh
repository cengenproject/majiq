#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=majiq2
#SBATCH -c 15
#SBATCH --mem=15G
#SBATCH --time=2-05:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


# MAJIQ workflow v3: in v2 grouped all the biological and technical replicates by neuron. In v3, group the technical replicates (if any) by biological replicate.
# MAJIQ v4: use the merged bams (no technical replicate), group the biological replicates by neuron
# v5: also perform tests for every neuron pair

# ----------------------------------------------------------------------------
# -----------------------    Parameters definitions   ------------------------
# ---------------------------------------------------------------------------- 

# Locations and parameters
#bam_dir='/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/bsn3'
bam_dir='/home/aw853/scratch60/2021-11-08_alignments'

WSversion='WS281'
references_dir='/gpfs/ycga/project/ysm/hammarlund/aw853/references'

outdir='/gpfs/ycga/project/ysm/hammarlund/aw853/counts/majiq/2021-11-09_outs'

gff_conversion="gtf2gff3" # or gffread



# ----------------------------------------------------------------------------
# ---------      No modification required below this line     ----------------
# ----------------------------------------------------------------------------


# Define references
ref_gff_name='c_elegans.PRJNA13758.'$WSversion'.canonical_geneset_'$gff_conversion'.gff3'

ref_gff=$references_dir'/'$WSversion'/'$ref_gff_name

cfg_file=$outdir/$(date +%F)_majiq.cfg


# the output subdirectories
mkdir -p $outdir/build
mkdir -p $outdir/psi
mkdir -p $outdir/logs
mkdir -p $outdir/deltapsi





# --------    End definitions, start computations    --------

set -e

echo
echo
echo "Starting majiq (config, build, quantif, test v5), $(date)"
echo
echo


# Create gff3 from gtf if needed
if [[ ! -f $ref_gff ]]
then
  echo "GFF3 file doesn't exist. Create it with $gff_conversion"
  exit 1
fi


echo "--------   List samples to process   --------"

samplelist=($(ls $bam_dir/*bam | xargs basename -a -s .bam))
declare -A sampByNeur #for the config file to the BUILD command
declare -A sampByNeur_path # for the PSI command

for smp in "${samplelist[@]}"
do
	if [[ $smp =~ ^([A-Z1-9ef]{2,4})r[0-9]{1,4}[t12]*$ ]]
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



echo "--------     Create config file     --------"

echo "[info]
bamdirs=$bam_dir
genome=$WSversion
strandness=forward
[experiments]" > $cfg_file

printf "%s\n" ${sampByNeur[@]} >> $cfg_file

echo "Config file generated and saved as $cfg_file"




# build ----

echo
echo "Running MAJIQ build on $(date) with $gff_conversion and autocreated conf file"
echo

source /home/aw853/bin/majiq/env/bin/activate




# --min-experiments  fraction or number of experiments that must pass filters [Default: 0.5]
#
###     SJ filters
# --minreads    minimum nb of reads for SJ [Default: 3]
# --minpos      min nb of read positions [Default: 2]
# --min-denovo  min nb reads on de-novo SJ [Default: 5]
#
###     Intron retention filters
# --irnbins           fraction of intronic read positions that must pass min-intronic-cov [Default: 0.5]
# --min-intronic-cov  min per-position normalized intronic readrate [Default: 0.01]
# 
# Simplifier options:
# 
#     --simplify simplification ignores junctions and introns with low usage [Default: -1] 

# Bootstrap coverage sampling:
# 
#     --markstacks p-value threshold used for dremoving read stacks [Default: 1e-07]
#     --m number of bootstraps to save in output SJ and MAJIQ files [Default: 30]

majiq build --nproc $SLURM_CPUS_PER_TASK \
			-o $outdir/build \
			--conf $cfg_file \
			--min-experiments 2 \
			--logger $outdir/logs \
			$ref_gff

echo "Finished the build step."
echo





# psi ----

echo "--------    PSI quantifications    --------"
echo "Running MAJIQ PSI on $(date) with conf $cfg_file"
echo


# --minreads min nb of reads at any positions in a LSV [Default: 10]
# --minpos min nb of start positions with at least 1 read in a LSV to considered. [Default: 3]
# --min-experiments Use to alter the threshold for group filters.
# --output-type {voila,tsv,all} Specify the type(s) of output files to produce
#   
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



echo "--------      DeltaPSI tests      --------"
echo "Running MAJIQ DeltaPSI on $(date) with conf $cfg_file"
echo

# get the neuron names in a non-associative array
all_neurs=(${!sampByNeur[@]})

for ((  i = 0; i < ${#all_neurs[@]}; i++ ))
do
  for (( j = $i+1; j < ${#all_neurs[@]}; j++))
	do
	  echo "##################################################################"
		
		neurA=${all_neurs[$i]}
		neurB=${all_neurs[$j]}
		
		echo "i: $i , j: $j ; Testing $neurA vs $neurB"
		echo
		
		
    # Mandatory arguments:
    # 
    #   -grp1 FILES1 [FILES1 ...]: Set of .majiq file[s] for the first condition
    #   -grp2 FILES2 [FILES2 ...]: Set of .majiq file[s] for the second condition
    #   -n/--names NAMES [NAMES ...]: _cond_id1_ _cond_id2_: group identifiers for grp1 and grp2 respectively.
    #   -o/--output OUTDIR: PSI output directory. It will contain the deltapsi.voila file once the execution is finished.
    # 
    # Optional arguments:
    # 
    #   --minreads MINREADS: Minimum number of reads [Default: 10]
    #   --minpos MINPOS: Minimum number of start positions [Default: 3]
    #   --min-experiments MIN_EXP: Use to alter the threshold for group filters.
    #   --binsize BINSIZE: The bins for PSI values. With a BINSIZE of 0.025 (default), we have 40 bins
    #   --default-prior: Use a default prior instead of computing it using the empirical data
    #   --prior-minreads PRIORMINREADS: Minimum number of reads (for the 'best set' calculation). [Default: 20]
    #   --prior-minnonzero PRIORMINNONZERO: Minimum number of positions for the best set.
    #   --prior-iter ITER: Max number of iterations of the EM
    #   --output-type {voila,tsv,all} Specify the type(s) of output files to produce [Default: all]

		majiq deltapsi \
		  -grp1 ${sampByNeur_path[$neurA]} \
		  -grp2 ${sampByNeur_path[$neurB]} \
		  -n $neurA $neurB \
		  -o $outdir/deltapsi \
		  --nprocs $SLURM_CPUS_PER_TASK \
		  --output-type all \
		  --logger $outdir/logs
		
		echo
		
	done
done



echo "------------ End tests ------------"
echo


# export relevant data ----
mv $outdir/build/splicegraph.sql $outdir/splicegraph.sql
mv $outdir/build/majiq.log $outdir/logs/majiq_build.log
mv $outdir/psi/psi_majiq.log $outdir/logs/majiq_psi.log

# to download on local computer for Voila and R analyses, then can be deleted from cluster
tar -czf $outdir/$(date +%y%m%d)_mjq_exprt.tar.gz $outdir/psi/ $outdir/splicegraph.sql


echo "All done $(date)"
