#!/bin/sh
INPUT_DIR='/home/ob219/scratch/pid/all_ai_exons_by_chr'
OUTPUT_DIR="${INPUT_DIR}/bcf_pid_only/"
BCF_TOOLS_BIN="/home/ob219/bin/bcftools-1.4/bcftools"
VCF_DIR="/scratch/WGS10K/data/release/latest/merged-vcf/no_hgmd/gnomad/"
SAMPLE_FILE="/home/ob219/scratch/ob219/pid/BRIDGE-PID_INDEX_CASES.txt"

for i in `\ls $INPUT_DIR/*.bed` ;do
file=$(basename $i .bed)
ofile="${OUTPUT_DIR}ai_${file}.bcf"
vcf="${VCF_DIR}chr${file}_agg3_dedup_vep_gnomad.bcf"
cmd="$BCF_TOOLS_BIN view -S $SAMPLE_FILE -R $i -O b -o $ofile -f PASS  $vcf"
echo "$cmd"
done
