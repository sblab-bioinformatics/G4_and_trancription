Differential binding analysis at active G4-positive promoters - DRB, TPL vs controls
================

After collecting read coverage at the regions of interest, the differential binding can be performed.

For G4 regions, several sequencing runs have been performed and they correspond to different biological replicates.

Two customised R script, one for [DRB](./DBA_DRB_Revised_2conditions_Nov3.R) and a second one for [TPL](./DBA_TPL_Revised_2conditions_Nov3.R) perform the differential analysis.

``` bash

## ------------------------------
## differential analysis DRB
## ------------------------------
rm(list = ls())
threshold_over_input=0
flag_BG4_at_TSS=1
file_subset_TSS='/Users/simeon01/Documents/Karen/20201023_revise_TPL_DRB/DRB_coverages/consensus_noDRB_DRB_1hr.merged_TSS_plus_minus500.bed'
dir_selected='/Users/simeon01/Documents/Karen/20201023_revise_TPL_DRB/DRB_coverages/only_2Cond_DRBL_updated_filter_over_input_thr0_TSS_plus_minus500'
source('./DBA_DRB_Revised_2conditions_Nov3.R')

rm(list = ls())
threshold_over_input=2
flag_BG4_at_TSS=1
file_subset_TSS='/Users/simeon01/Documents/Karen/20201023_revise_TPL_DRB/DRB_coverages/consensus_noDRB_DRB_1hr.merged_TSS_plus_minus500.bed'
dir_selected='/Users/simeon01/Documents/Karen/20201023_revise_TPL_DRB/DRB_coverages/only_2Cond_DRBL_updated_filter_over_input_thr2_TSS_plus_minus500'
source('/Users/simeon01/Documents/Karen/20201023_revise_TPL_DRB/DBA_DRB_Revised_2conditions_Nov3.R')


## ------------------------------
## differential analysis TPL
## ------------------------------

rm(list = ls())
threshold_over_input=2
flag_BG4_at_TSS=1
file_subset_TSS='/Users/simeon01/Documents/Karen/20201023_revise_TPL_DRB/TPL_coverages/consensus_noTPL_10uM_TPL_2hr.merged_TSS_plus_minus500.bed'
dir_selected='/Users/simeon01/Documents/Karen/20201023_revise_TPL_DRB/TPL_coverages/only_2Cond_TPL_updated_filter_over_input_thr2_TSS_plus_minus500'
source('/Users/simeon01/Documents/Karen/20201023_revise_TPL_DRB/DBA_TPL_Revised_2conditions_Nov3.R')
```
