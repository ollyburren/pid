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
3. findPoss.R - This script generates the lists that we sent around.
4. parselofVCF.R - Rscript that Filters a gene for lof variants.                                                                                                                                                                                                             

## assignment of SE calls to PCHiC tissues

|    | Javierre Cell Type     | Super Enhancer Cell Type assignment                                                     |
|---:|------------------------|-----------------------------------------------------------------------------------------|
|  1 | Monocytes              | Peaks.Monocytes_merged_1e-5_broad_peaks.SuperEnhancers.bed                              |
|  2 | Macrophages_M0         | Peaks.Macrophage_merged_1e-5_broad_peaks.SuperEnhancers.bed                             |
|  3 | Macrophages_M1         | Peaks.InflammatoryMacrophage_merged_1e-5_broad_peaks.SuperEnhancers.bed                 |
|  4 | Macrophages_M2         | Peaks.AlternativelyActivatedMacrophage_merged_1e-5_broad_peaks.SuperEnhancers.bed       |
|  5 | Neutrophils            | Peaks.Neutrophils_Merged_1e-5_broad_peaks.SuperEnhancers.bed                            |
|  6 | Megakaryocytes         | Peaks.MK_merged_1e-5_broad_peaks.SuperEnhancers.bed                                     |
|  7 | Endothelial_precursors | Peaks.BOEC_progenitors_1e-5_broad_peaks.SuperEnhancers.bed                              |
|  8 | Erythroblasts          | Peaks.EB_Merged_1e-5_broad_peaks.SuperEnhancers.bed                                     |
|  9 | Foetal_thymus          |                                                                                         |
| 10 | Naive_CD4              | Peaks.CD4_CentralMemoryTcell_merged_1e-5_broad_peaks.SuperEnhancers.bed                 |
| 11 | Total_CD4_MF           |                                                                                         |
| 12 | Total_CD4_Activated    | Peaks.CD4Proliferating_Merged_1e-5_broad_peaks.SuperEnhancers.bed                       |
| 13 | Total_CD4_NonActivated | Peaks.CD4NonActive_Merged_1e-5_broad_peaks.SuperEnhancers.bed                           |
| 14 | Naive_CD8              | Peaks.CentralmemoryCD8positiveAlphaBetaTcell_merged_1e-5_broad_peaks.SuperEnhancers.bed |
| 15 | Total_CD8              | Peaks.CD8positiveAlphaBetaTcell_merged_1e-5_broad_peaks.SuperEnhancers.bed              |
| 16 | Naive_B                | Peaks.NaiveBcell_merged_1e-5_broad_peaks.SuperEnhancers.bed                             |
| 17 | Total_B                | Peaks.CentralmemoryCD8positiveAlphaBetaTcell_merged_1e-5_broad_peaks.SuperEnhancers.bed |
