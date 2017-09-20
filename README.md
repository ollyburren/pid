# pid
Code for PID analysis

## createRDS

These are R scripts that extract PID data without a molecular diagnosis

1. ./R/analysis/createRDS/get_deletions.R - Create and RDS file for SV > 0.03 (3%). bnd - break point, del - deletion, dup - duplication, ins - insertion, inv - inversion.
2. ./R/analysis/createRDS/get_exons.R - using gencode v26 lifted over to build 37 get a list of exonic regions - write a list of bcftools to extract the list of variants for pid samples that overlap these and save as a bcf file. We run this on HPC cluster.
3. ./R/analysis/createRDS/get_genes.R - using gencode v26 lifted over to build 37 get a list of genic regions.
4. ./R/analysis/createRDS/get_pchic.R - using Javierre et al data integrate Hnisz et al. data - create GRanges objects for other ends that overlap; se - superenhancer, enh - enhancer, se_enh - union of previous.
5. ./R/analysis/createRDS/get_promoters.R - using gencode v26 lifted over to build 37 get a list of promoter regions (500BP downstream of any transcriptional start site)

## annotate

1. shallow_overlap.R - Takes the RDS objects above and integrates them. It does this by finding genes that are interesting by promoter,exon or pchic,se overlap and then interrogating the protein coding space for individuals wotjh co-occurring rare variants in that gene (operates on the HPC using bcftools).
2. central_functions.R - Deprecated.
3. findPoss.R - Deprecated ?
4. parselofVCF.R - Rscript that Filters a gene for lof variants.
5.                                                                                                                                                                                                                                                   
