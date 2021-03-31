Data preparation for differential analysis
================
Angela Simeone

Two main experiments with hypoxia-normoxia have been conducted in the [study](hppt:link). In the first experiment, K562 cells have been exposed to (A) hypoxic and (B) normoxic conditions profiled (G4-ChIP-seq, Pol2-ChIP-seq, ATAC-seq). Regions with enrichments (confirmed across multiple replicates) have been identified for each treatment case. Upon normoxia, G4 regions overlapping Pol2 enrichments and also placed in promoters (TSS+-500bp) have been selected to perform differential binding analysis. Those regions represent active promoter loci with a formed G4 structure.

Briefly, the steps perfomed here are:

1.  definintion of ROIs;
2.  count reads at ROIs (read coverage);
3.  fetching total library size (for each individual library used in the analysis
4.  differential binding analysis.

``` bash
### ---- promoters (TSS+-500) ---- ### 
bedtools slop -i gencode.v28.TSS.bed -g hg38.sorted.genome -s -b 500 > gencode.v28.TSS_plus_minu500.bed
wc -l gencode.v28.TSS_plus_minu500.bed
# 9217 G4_in_pol2_TSS_minus500.bed

### ---- define ROIs ---- ### 
TSS_500=/reference_genomes/hg38/gencode.v28.TSS_plus_minus_500bp.bed
g4_normo=normoxia.bio2.bed
g4_hypo=hypoxia.bio2.bed
pol2_normo=normoxia.multi3.bed
pol2_hypo=hypoxia.multi3.bed

# create reference of G4_in pol2 TSS plus minus 500
intersectBed -a $g4_normo -b $pol2_normo -wa | sort | uniq | intersectBed -a - -b $TSS_500 -wa | sort | uniq > G4_in_pol2_TSS_minus500.bed


### ---- coverage at G4_in_pol2_TSS_plus_minus500bp ROIs  for experiments hypoxia normoxia ---- ### 
## coverage at G4_in_pol2_TSS_plus_minus500bp
bam_atac=/scratchb/sblab/simeon01/Atac_karen_hypo_normo_SLX18056_18895
bam_BG4=/scratcha/sblab/simeon01/Data/20200108_Karen_hypoxia_k562_BG4/SLX-18891/trimmed/aligned/merged_runs
bam_pol2=/scratcha/sblab/simeon01/Data/20200216_karen_pol2chip_hypoxia/SLX-18974/trimmed/aligned
path_out=/scratcha/sblab/simeon01/Data/20201104_karen_hypoxia_normoxia_pol2_atac_bg4
mkdir $path_out
region1_BG4_pol2=/scratcha/sblab/simeon01/Data/20201023_revise_atac_pol2_hypox_normox/G4_in_pol2_TSS_minus500.bed

bams_hypo_normo=($bam_atac $bam_BG4 $bam_pol2)
for path in ${bams_hypo_normo[@]}
  do
  cd $path
  for bam in *markduplicates.bam
    do
    basename_region1=`basename $region1_BG4_pol2`
    cmd_coverage_G4_1="bedtools coverage -a $region1_BG4_pol2 -b ${bam} -counts > $path_out/${bam%%.bam}.dba_consensus_${basename_region1%%.bed}.coverage.bed"
    echo $cmd_coverage_G4_1
    sbatch --time 01:00:00 --mem 4G --wrap "$cmd_coverage_G4_1"
    
    done
  echo " ========= "
done

# recover files with total number of reads and move them into the same folder as the coverage files
for path in ${bams_hypo_normo[@]}
  do
  cd $path
  cp *stat5 $path_out
  
done


### ---- coverage at G4_in_pol2_TSS_plus_minus500bp ROIs  for experiments hypoxia normoxia in presence of pyPDS or DMSO ---- ### 
# do the same as above using the same genomics regions but referring to the experiments pyPDS, hypoxia and normoxia
bam_atac_pyPDS=/scratcha/sblab/simeon01/Data/20200708_karen_atac_pyPDS_hypoxia_rep/merged_atac_SLX-19339_SLX-19340
bam_pol2_pyPDS=/scratcha/sblab/simeon01/Data/20200727_karen_pol2_hypoxia_pds_normoxia_dmso/merged_pol2chip_SLX-19338_SLX-19411
#bam_BG4_pyPDS=
bams_hypo_normo_pds=($bam_atac_pyPDS $bam_pol2_pyPDS)
path_out=/scratcha/sblab/simeon01/Data/20201104_karen_pyPDS_hypoxia_normoxia_pol2_atac
mkdir $path_out
region1_BG4_pol2=/scratcha/sblab/simeon01/Data/20201023_revise_atac_pol2_hypox_normox/G4_in_pol2_TSS_minus500.bed
region1_BG4_pol2_second=/scratcha/sblab/simeon01/Data/20201019_BG4hypoxia_normoxia_SLX-198891_SLX-198892/overlap_of_merge_BG4_POL2_normo_hypo.bed
cd $path_out
ls


for path in ${bams_hypo_normo_pds[@]}
  do
  cd $path
  cp *stat5 $path_out
  for bam in *markduplicates.bam
    do
    basename_region1=`basename $region1_BG4_pol2`
    cmd_coverage_G4_1="bedtools coverage -a $region1_BG4_pol2 -b ${bam} -counts > $path_out/${bam%%.bam}.dba_consensus_${basename_region1%%.bed}.coverage.bed"
    echo $cmd_coverage_G4_1
    sbatch --time 01:00:00 --mem 4G --wrap "$cmd_coverage_G4_1"
    
    basename_region2=`basename $region1_BG4_pol2_second`
    cmd_coverage_G4_2="bedtools coverage -a $region1_BG4_pol2_second -b ${bam} -counts > $path_out/${bam%%.bam}.dba_consensus_${basename_region2%%.bed}.coverage.bed"
    echo $cmd_coverage_G4_1
    sbatch --time 01:00:00 --mem 4G --wrap "$cmd_coverage_G4_2"
    done
  echo " ========= "
done
```

For normoxia, chech overlap of G4 regions with ATAC-regions.

``` bash

## overlap between BG4 regions and ATAC regions in normoxia
```

``` bash

g4_normo=/scratcha/sblab/simeon01/Data/20201019_BG4hypoxia_normoxia_SLX-198891_SLX-198892/normoxia.bio2.bed
g4_hypo=/scratcha/sblab/simeon01/Data/20201019_BG4hypoxia_normoxia_SLX-198891_SLX-198892/hypoxia.bio2.bed
pol2_normo=/scratcha/sblab/simeon01/Data/20201019_BG4hypoxia_normoxia_SLX-198891_SLX-198892/normoxia.multi3.bed
pol2_hypo=/scratcha/sblab/simeon01/Data/20201019_BG4hypoxia_normoxia_SLX-198891_SLX-198892/hypoxia.multi3.bed
#atac_normo=/scratchb/sblab/simeon01/Atac_karen_hypo_normo_SLX18056_18895/normoxia_1hr_wholeK562_fragment.multi2.bed
#atac_hypo=/scratchb/sblab/simeon01/Atac_karen_hypo_normo_SLX18056_18895/hypoxia_1hr_wholeK562_fragment_1hr.multi2.bed

atac_normo=/scratcha/sblab/simeon01/Data/20201019_BG4hypoxia_normoxia_SLX-198891_SLX-198892/normoxia_1hr_wholeK562_fragment.multi2.bed
atac_hypo=/scratcha/sblab/simeon01/Data/20201019_BG4hypoxia_normoxia_SLX-198891_SLX-198892/hypoxia_1hr_wholeK562_fragment_1hr.multi2.bed

cd /scratcha/sblab/simeon01/Data/20201019_BG4hypoxia_normoxia_SLX-198891_SLX-198892/

mkdir intersection_G4_atac_normoxia
intervene venn -i $g4_normo $atac_normo -o ./intersection_G4_atac_normoxia
mkdir intersection_G4_atac_hypoxia
intervene venn -i $g4_hypo $atac_hypo -o ./intersection_G4_atac_hypoxia
mkdir intersection_G4_atac_normoxia_hypoxia
intervene venn -i $g4_normo $atac_normo $g4_hypo $atac_hypo -o ./intersection_G4_atac_normoxia_hypoxia

mkdir intersection_G4_atac_pol2_normoxia
intervene venn -i $g4_normo $atac_normo $pol2_normo -o ./intersection_G4_atac_pol2_normoxia
mkdir intersection_G4_atac_pol2_hypoxia
intervene venn -i $g4_hypo $atac_hypo $pol2_hypo -o ./intersection_G4_atac_pol2_hypoxia
mkdir intersection_G4_atac_pol2_normoxia_hypoxia
intervene venn -i $g4_normo $atac_normo $g4_hypo $atac_hypo $pol2_normo $pol2_hypo -o ./intersection_G4_atac_pol2_normoxia_hypoxia

bed_to_use=($g4_normo $g4_hypo $pol2_normo $pol2_hypo $atac_normo $atac_hypo)
oqs=/scratcha/sblab/simeon01/reference_genomes/OQs/OQ_hits.lifted_hg19_to_hg38_no_.bed

for file in ${bed_to_use[@]}
do
  a=`wc -l $file |awk '{print $1}'`
  b=`intersectBed -a $file -b $oqs -wa | sort | uniq | wc -l |awk '{print $1}'`
  c=$file
  
  echo $c "\t" $a "\t" $b
done
```
