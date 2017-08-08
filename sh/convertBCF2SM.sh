#!/bin/sh
INPUT_DIR='/home/ob219/scratch/pid/all_ai_exons_by_chr/bcf_pid_only/'
OUTPUT_DIR="${INPUT_DIR}/../snpMatrix_pid_only/"

for i in `\ls $INPUT_DIR/*.bcf` ;do
file=$(basename $i .bcf)
ofile="${OUTPUT_DIR}${file}.RData"
cmd="Rscript --vanilla /home/ob219/git/pid/R/enrich/convert_bcf2snpMatrix.R -f $i -o $ofile"
echo "$cmd"
done

INPUT_DIR='/home/ob219/scratch/pid/all_ht_exons_by_chr/bcf_pid_only/'
OUTPUT_DIR="${INPUT_DIR}/../snpMatrix_pid_only/"

for i in `\ls $INPUT_DIR/*.bcf` ;do
file=$(basename $i .bcf)
ofile="${OUTPUT_DIR}${file}.RData"
cmd="Rscript --vanilla /home/ob219/git/pid/R/enrich/convert_bcf2snpMatrix.R -f $i -o $ofile"
echo "$cmd"
done
