MAJIQ analysis conducted using script `src/majiq.sh`. To produce correct GFF format, this script additionally calls `src/add_gene_names.awk`.

After MAJIQ has run, results are analyzed (and figures created) with `R/deltapsi.R`.

Script `R/format_literature.R` was used to compile the list of DAS genes used in the main script.


Other, more exploratory scripts are stored in separate branch `exploratory`, including VOILA modulizer separately ran with `src/modulize.sh` (not included in paper).


