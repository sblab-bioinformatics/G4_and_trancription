Differential binding analysis ATAC-seq at G4-positive active promoters
================

After collecting read coverage at the regions of interest, the differential binding can be performed.

Different biological replicates have been sequenced together and 2 sequencing runs were performed to achieve good sequencing depth.

R customised script [DBA\_at\_G4\_pol2\_tss\_plus\_minus500](./DBA_at_G4_pol2_tss_plus_minus500.R) performing the differential analysis is [here](./DBA_at_G4_pol2_tss_plus_minus500.R)

``` bash


## ==================== ==================== ==================== ====================
## differential analysis ATAC-seq
## ==================== ==================== ==================== ====================

## ATAC ##################################
rm(list = ls())
path1='/Users/simeon01/Documents/Karen/20201104_karen_hypoxia_normoxia_pol2_atac_bg4'
pattern_file_extension="dba_consensus_overlap_of_merge_BG4_POL2_normo_hypo.coverage.bed$"
selection_files_pattern="SLX-18056_SLX-18895"
pattern_folder_prefix="DBA_atac_in_G4_pol2_TSS_minus_plus500_regions_thr0_updated"
case=1
input_flag=0
threshold_over_input=2
flag_BG4_at_TSS=0
file_subset_TSS='/Users/simeon01/Documents/Karen/20201104_karen_hypoxia_normoxia_pol2_atac_bg4/overlap_of_merge_BG4_POL2_normo_hypo.at_TSS_plus_minus_500bp.bed'
source('/Users/simeon01/Documents/Karen/20201104_karen_hypoxia_normoxia_pol2_atac_bg4/DBA_at_G4_pol2_tss_plus_minus500.R')
```
