README
================


This repository collect details about data and computational methods developed for the paper "Promoter G-quadruplex folding precedes transcription and is controlled by chromatin".

DNA G-quadruplexes (G4) and transcription
=========================================

This study investigates the interplay between DNA G4 and transcription. Data used in this study have been deposited on GEO under the accession [GSE162299](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162299).

K562 cells have been profiled for DNA G4 structures (G4-ChIP-seq) and Pol2 (Pol2-ChIP-seq) occupancy in normoxia. Additionally accessibility have been profiled by ATAC-seq. Cells have been exposed initially to two chemicals (DRB, TPL) respectively and DNA G4s have been profiled. 
Next, after explosing cells to hypoxia, G4s, Pol2 and accessibility have been profiled. 
Finally, pyPDS (a G4-stabilising ligand) has been used in conjunction to hypoxic conditions: in this experimental case (pyPDS, hypoxia) Pol2 and accessibility have been profiled and compared to control cases (no pyPDS - DMSO, normoxia).

In this repository we report the general processing steps followed for the analysis of the sequencing data. In addition we report the custom code for the higher level analysis (differential binding anlysis, signal density analysis at TSS).



#### Basic data processing

-   [Processing of G4-ChIP-seq](./G4_seq_processing.md)

-   [Processing of Pol2-ChIP-seq](./Pol2_seq_processing.md)

-   [Processing of ATAC-ChIP-seq](./atac_seq_processing.md)

#### Higher level analysis

Regions of interest for differential analysis DRB and TPL:

-   [Define regions of interest for differential binding analysis and compute read coverage in those regions](./data_preparation_DRB_TPL.md)

    DRB

    -   [Differential binding analysis G4-ChIP-seq DRB](./wrapper_for_DBA_G4_DRB_TPL.md)

    TPL

    -   [Differential binding analysis G4-ChIP-seq TPL](./wrapper_for_DBA_G4_DRB_TPL.md)

Regions of interest for differential analysis Hypoxia - Normoxia and Hypoxia - Normoxia - pyPDS - DMSO:

-   [Define regions of interest for differential binding analysis and compute read coverage in those regions](./data_preparation.md)

    hypoxia vs normoxia

    -   [Differential binding analysis G4-ChIP-seq hypoxia vs normoxia](./wrapper_for_DBA_G4.md)

    -   [Differential binding analysis Pol2-ChIP-seq hypoxia vs normoxia](./wrapper_for_DBA_pol2_hypo_normo.md)

    -   [Differential binding analysis ATAC-seq hypoxia vs normoxia](./wrapper_for_DBA_atac_hypo_normo.md)

    hypoxia - normoxia - pyPDS - DMSO

    -   [Differential binding analysis Pol2-ChIP-seq hypoxia normoxia pyPDS DMSO](./wrapper_for_DBA_Pol2_hypo_normo_pyPDS_DMSO.md)

    -   [Differential binding analysis ATAC-seq hypoxia vs normoxia pyPDS DMSO](./wrapper_for_DBA_atac_hypo_normo_pyPDS_DMSO.md)

Genomic regions fold enrichments

-   [Fold-enrichments at genomic regions (exons, introns, promoters, intergenic, etc ...)](./assess_fold_enrichment_genomic_features.md)

Density signal and meta-genes signal at regions of interest

-   [Heatmap and TSS signal](./heatmaps_and_densities_maps.md)

Validtion of hypoxia effect on G4 signal in U2OS cells

-   [Basic processing and differential analysis](./processing_U2OS_dba.md)

#### Software, tools and environment used

<table style="width:100%;">
<colgroup>
<col width="50%" />
<col width="45%" />
<col width="4%" />
</colgroup>
<thead>
<tr class="header">
<th>Step</th>
<th>Software name and version</th>
<th></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>fastqc</td>
<td>FastQC v0.11.7</td>
<td></td>
</tr>
<tr class="even">
<td>adapter trimmin</td>
<td>cutadapt version 1.16</td>
<td></td>
</tr>
<tr class="odd">
<td>alignment</td>
<td>bwa mem 0.7.17-r1188</td>
<td></td>
</tr>
<tr class="even">
<td>duplicate marking</td>
<td>picard-2.20.3</td>
<td></td>
</tr>
<tr class="odd">
<td>bam file indexing, sorting and handeling</td>
<td>samtools Version: 1.8 (using htslib 1.8)</td>
<td></td>
</tr>
<tr class="even">
<td>bigWig track generation</td>
<td>bamCoverage 3.3.0 (deepTools 3.3.0)</td>
<td></td>
</tr>
<tr class="odd">
<td>peak calling</td>
<td>macs2 2.1.2</td>
<td></td>
</tr>
<tr class="even">
<td>bed files processing and manipulation</td>
<td>bedtools v2.27.1</td>
<td></td>
</tr>
<tr class="odd">
<td></td>
<td></td>
<td></td>
</tr>
<tr class="even">
<td>System info:</td>
<td></td>
<td></td>
</tr>
<tr class="odd">
<td>Linux kernel version</td>
<td>3.10.0-1127.18.2.el7.x86_64</td>
<td></td>
</tr>
<tr class="even">
<td>cluster management and job scheduling system</td>
<td>slurm 20.02.4</td>
<td></td>
</tr>
</tbody>
</table>

R session info:

``` bash
print(sessionInfo())
R version 3.6.1 (2019-07-05)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS High Sierra 10.13.6

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

locale:
[1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] ComplexHeatmap_2.0.0 dplyr_1.0.2          edgeR_3.26.8         limma_3.40.6        

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.5          pillar_1.4.6        compiler_3.6.1      RColorBrewer_1.1-2  tools_3.6.1         digest_0.6.25       evaluate_0.14      
 [8] lifecycle_0.2.0     tibble_3.0.3        lattice_0.20-41     clue_0.3-57         pkgconfig_2.0.3     png_0.1-7           rlang_0.4.7        
[15] rstudioapi_0.11     yaml_2.2.1          parallel_3.6.1      xfun_0.17           knitr_1.29          cluster_2.1.0       generics_0.0.2     
[22] GlobalOptions_0.1.2 vctrs_0.3.4         locfit_1.5-9.4      tidyselect_1.1.0    glue_1.4.2          R6_2.4.1            GetoptLong_1.0.2   
[29] rmarkdown_2.3       purrr_0.3.4         magrittr_1.5        ellipsis_0.3.1      htmltools_0.5.0     splines_3.6.1       shape_1.4.5        
[36] circlize_0.4.10     colorspace_1.4-1    crayon_1.3.4        rjson_0.2.20   
```
