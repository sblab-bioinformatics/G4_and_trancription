Differential binding analysis ATAC-seq at G4-positive active promoters - hypoxia, normoxia, pyPDS, DMSO
================

After collecting read coverage at the regions of interest, the differential binding is performed.

Different biological replicates have been sequenced together and 2 sequencing runs were performed to achieve good sequencing depth.

R customised script [DBA\_at\_G4\_pol2\_tss\_plus\_minus500](./sDBA_additional_analysis_pyPDS_hypo_normo_REVISED.r) performing the differential analysis is [here](./sDBA_additional_analysis_pyPDS_hypo_normo_REVISED.r)

``` bash

## ==================== ==================== ==================== ====================
## differential analysis ATAC-seq
## ==================== ==================== ==================== ====================
rm(list = ls())
path1="/Users/simeon01/Documents/Karen/20201104_karen_pyPDS_hypoxia_normoxia_pol2_atac"
pattern_file_extension="dba_consensus_G4_in_pol2_TSS_minus500.coverage.bed$"
selection_files_pattern="ATAC"
threshold_over_input=0
pattern_folder_prefix="DBA_ATAC_in_G4_pol2_TSS_minus_plus500_regions_thr0"
filter_pyPDSnormo=0
input_flag=0
flag_batch=0
case=2
flag_BG4_at_TSS=0
file_subset_TSS='/Users/simeon01/Documents/Karen/20201104_karen_hypoxia_normoxia_pol2_atac_bg4/overlap_of_merge_BG4_POL2_normo_hypo.at_TSS_plus_minus_500bp.bed'
source("~/sDBA_additional_analysis_pyPDS_hypo_normo_REVISED.r")
```
