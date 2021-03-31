G4-ChIP-seq processing
================

Sequencing data processing included the following steps:

-   adaptor trimming
-   alignment
-   merge resequenced libraries (when present)
-   duplicate removal
-   peak calling

Software, tools and environment used:

<table style="width:100%;">
<colgroup>
<col width="50%" />
<col width="45%" />
<col width="4%" />
</colgroup>
<thead>
<tr class="header">
<th>Software name and version</th>
<th>step</th>
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

Quality control with fastqc

``` bash
mkdir fastqc
for file in *.fq.gz
do
sbatch --mem 4G --wrap "fastqc --noextract --nogroup -q -o fastqc/ $file"
done
```

Adaptor trimming

``` bah
mkdir trimmed
for fq in *.fq.gz
do
bname=`basename $fq`
sbatch -o %j.$bname.tmp.out -e %j.$bname.tmp.err --mem 16000 --wrap "cutadapt -q 20 -O 3 -a CTGTCTCTTATACACATCT -o trimmed/${bname%%.fq.gz}.trimmed.fq.gz $fq"
done
```

Aligment of trimmed reads

``` bash

path=~/trimmed
g='~/reference_genomes/hg38/hg38_selected.fa'
w='~/reference_genomes/hg38/hg38.whitelist.sorted.bed'
mkdir aligned
for f in *gz
do
sbatch -o %j.$f.tmp.out -e %j.$f.tmp.err --mem 16000 --wrap "bwa mem -M $g $f | samtools view -S -u -F2304 -q 10 -L $w - | samtools sort -@ 8 - > $path/aligned/${f%%.trimmed.fq.gz}.hg38.sort.bam"
done
```

Duplicates removal

``` bash
cd ~/trimmed/aligned
for bam in *hg38.sort.bam
do
sbatch -o %j.$bam.tmp.out -e %j.$bam.tmp.err --mem 32000  --wrap "java -Xmx7g -jar ~/picard-2.20.3.jar MarkDuplicates INPUT=$bam OUTPUT=${bam%%.bam}.markduplicates.bam REMOVE_DUPLICATES=true  AS=true METRICS_FILE=${bam%%.bam}.markduplicates_metrics.txt"
done
```

Indexing markduplicates bam

``` bash
cd ~/trimmed/aligned
for bam in *.sort.markduplicates.bam
do
sbatch -o %j.out.log --mem 16G --wrap "samtools index $bam"
done
```

Generate stats

``` bash
cd ~/trimmed/aligned
for file in *hg38.sort.bam
do
bam_hg38=$file
bam_hg38_nodup=${file%%.bam}.markduplicates.bam
ls -lth $bam # check bam exist
echo "====="
ls -lth $bam_hg38_nodup # check bam exist

sbatch -o %j.out.stat --mem 8G --wrap "samtools view -c -F 260 $bam_hg38 >  ${file%%.bam}.stat2"
sbatch -o %j.out.stat --mem 8G --wrap "samtools view -c -F 260 $bam_hg38_nodup >  ${file%%.bam}.stat5"
done

cd ~/trimmed/aligned

touch pyPDS_hypoxia_BG4_SLX-19336.out.stat
for bam in *.sort.markduplicates.bam
do
base_name=${bam%%.markduplicates.bam}
#base_name=${bam/.sort.markduplicates.bam/.sort.bam}
files_stat=`ls *.stat* | grep $base_name`
echo $bam > temp
cat $files_stat >> temp
paste -d '\t' pyPDS_hypoxia_BG4_SLX-19336.out.stat temp >> temp2
mv temp2 pyPDS_hypoxia_BG4_SLX-19336.out.stat
done
```

Generate tracks for human

``` bash

cd ~/trimmed/aligned
for bam in *.sort.markduplicates.bam
do
tot_r_hg38=`cat ${bam%%.markduplicates.bam}.stat5`
scal_factor_hg38=$(awk -v m=$tot_r_hg38 'BEGIN { print 1000000/m }')
echo $scal_factor_hg38
echo sbatch -o %j.out.log --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"
sbatch -o %j.out.log --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"
#sbatch -o %j.out.log --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 50 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w50.rpm.bw"
#sbatch -o %j.out.log --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 100 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w100.rpm.bw"
done
```

Call macs2 peaks with default options and also broad option

``` bash
cd ~/trimmed/aligned

rep=(1 3)
condition=(1uM_pyPDS_2hr_hypoxia_biorep1 1uM_pyPDS_2hr_hypoxia_biorep2 1uM_pyPDS_2hr_hypoxia_biorep3 DMSO_control_2hr_hypoxia_biorep1 DMSO_control_2hr_hypoxia_biorep2 DMSO_control_2hr_hypoxia_biorep3)
counter_l=1
out_dir_narrow='~/trimmed/aligned/macs2_individual_rep'
out_dir_broad='~/trimmed/aligned/macs2_broad_individual_rep'
mkdir $out_dir_narrow
mkdir $out_dir_broad
for currcon in ${condition[@]};
do
  for i in ${rep[@]}; 
  do
  echo "${c}"
  echo "$i"
  #echo ${c}*ChIP_biorep$i*.sort.markduplicates.bam
  t=`ls ${currcon}*ChIP$i*.sort.markduplicates.bam`
  c=`ls ~/Data/20200619_karen_bg4_pyPDS_hypoxia_rep/SLX-19269/trimmed/aligned/${currcon}*Input_*.sort.markduplicates.bam`
  
  ls -lht $t
  ls -lht $c
  #sbatch --mem 8000 -J macs --wrap "macs2 callpeak --keep-dup all --broad -t $t -c $c -n $out_dir_broad/${currcon}_ChIP_biorep${i}_broad" # human is default
  #sbatch --mem 8000 -J macs --wrap "macs2 callpeak --keep-dup all -t $t -c $c -n $out_dir_narrow/${currcon}_ChIP_biorep${i}" # human is default
  
  #sbatch --mem 8000 -J macs --wrap "macs2 callpeak --keep-dup all -t $t -c $c -p 1e-3 -n $out_dir_narrow/${currcon}_ChIP_biorep${i}_1e-3" # human is default
  #sbatch --mem 8000 -J macs --wrap "macs2 callpeak --keep-dup all --broad -t $t -c $c -p 1e-3 -n $out_dir_broad/${currcon}_ChIP_biorep${i}_broad_1e-3" # human is default
  echo "====="
  done    
done
```

Narrow Peaks consensus

``` bash
# create a new macs2 folder with also previous peaks
mkdir ~/trimmed/aligned/macs2_individual_rep_including_previous_run
cd ~/trimmed/aligned/macs2_individual_rep_including_previous_run

cp ~/trimmed/aligned/macs2_individual_rep/*narrowPeak .
cp ~/Data/20200619_karen_bg4_pyPDS_hypoxia_rep/SLX-19269/trimmed/aligned/macs2_individual_rep/*narrowPeak .


out_consensum=~/trimmed/aligned/macs2_individual_rep_including_previous_run/consensus_narrow_peaks

condition=(1uM_pyPDS_2hr_hypoxia DMSO_control_2hr_hypoxia)

mkdir $out_consensum

for case in ${condition[@]}
do
echo $case
bed_files=`ls *$case*.narrowPeak`

echo $bed_files

multiIntersectBed -i $bed_files | awk '{if($4>=6) print $0}' | sortBed -i | mergeBed -i - > $out_consensum/$case.multi6.bed
echo " ================ . ================"
done
```

Overlap with OQs

``` bash
cd ~/trimmed/aligned/macs2_individual_rep_including_previous_run

oqs=~/reference_genomes/OQs/OQ_hits.lifted_hg19_to_hg38_no_.bed
jochen=~/Data/Jochen_K562_peaks/20180108_K562_async_rep1-3.mult.5of8.hg19_lifted_hg38.bed

for file in *narrowPeak
do
n_peaks=`wc -l $file | awk '{print $1}'`
n_peaks_in_oqs=`intersectBed -a $file -b $oqs -wa | sort | uniq | wc -l`
n_peaks_in_jochen=`intersectBed -a $file -b $jochen -wa | sort | uniq | wc -l`
echo  $file$'\t'$n_peaks$'\t'$n_peaks_in_oqs$'\t'$n_peaks_in_jochen
echo  $file$'\t'$n_peaks$'\t'$n_peaks_in_oqs$'\t'$n_peaks_in_jochen >> N_peaks_pyPDS_hypoxia_BG4_SLX-19336.individualLib.txt
done

out_consensum=~/trimmed/aligned/macs2_individual_rep_including_previous_run/consensus_narrow_peaks
cd $out_consensum

for file in *multi*.bed
do
n_peaks=`wc -l $file | awk '{print $1}'`
n_peaks_in_oqs=`intersectBed -a $file -b $oqs -wa | sort | uniq | wc -l`
n_peaks_in_jochen=`intersectBed -a $file -b $jochen -wa | sort | uniq | wc -l`
echo  $file$'\t'$n_peaks$'\t'$n_peaks_in_oqs$'\t'$n_peaks_in_jochen
echo  $file$'\t'$n_peaks$'\t'$n_peaks_in_oqs$'\t'$n_peaks_in_jochen >> N_peaks_pyPDS_hypoxia_BG4_SLX-19336.multi.txt
done
```

Coverage for differential binding analysis
------------------------------------------

In order to perfom the differential analysis computed the coverate at the merge of all consensus peaks (observed in each experimental conditions).

For this purpose we follower those steps:

1.  merge of all confirmed peaks
2.  compute coverage at the final consensus
3.  perform differential analysis in EdgeR importing also the individual library sizes

``` bash
# select consensus peaks from before lockdown - set0
# select consensus peaks 1st june sequenncing - set1
# select consensus peaks 2nd june sequencing - set2

set0=~trimmed/aligned/macs2_individual_rep/consensus_narrow_peaks/*bio*multi2.bed
#set1=/trimmed2/aligned/macs2_individual_rep/consensus_narrow_peaks/*bio*multi2.bed # this might contain the additional sequencing run performed for the additional biological replicate
#set2=~/trimmed3/aligned/macs2_individual_rep_including_previous_run/consensus_narrow_peaks/*.multi6.bed

# merge peaks for coverage -

mkdir ~/Data/final_bg4_coverage
cd ~/Data/final_bg4_coverage
cat $set0 $set1 $set2 | sortBed -i - | mergeBed -i - > dba_consensus.20201001.experimental_conditions.bed


# coverag for DBA

bed_consensus=~/Data/final_bg4_coverage/dba_consensus.20201001.experimental_conditions.bed
path_out=~/Data/final_bg4_coverage
path_bam0=~/trimmed/aligned
path_bam1=~/trimmed2/aligned
path_bam2=~/trimmed3/aligned

all_bam_paths=($path_bam0 $path_bam1 $path_bam2)

cd $path_out

for path in ${all_bam_paths[@]}
do
  cd $path
  for bam in *markduplicates.bam
  do
    cmd_coverage="bedtools coverage -a $bed_consensus -b ${bam} -counts > $path_out/${bam%%.bam}.dba_consensus.20201001.experimental_conditions.bed"
    echo $cmd_coverage
    echo "*****"
    sbatch --mem 4G --wrap "$cmd_coverage"
  done
cd $path_bed
echo " ========= "
cd $path_out
done
```

After extracting the coverages, proceed with differential analysis. For this we used R.

-   [Differential analysis TPL treatment (4 experimental conditions)](./Rscripts/DBA_TPL_all_samples_revised_Oct2020.R)

-   [Differential analysis TPL treatment (2 experimental conditions)](./Rscripts/DBA_TPL_all_samples_revised_ONLY2conditions_Oct2020.R)

<!-- -->

    > print(.libPaths())
    [1] "/Library/Frameworks/R.framework/Versions/3.6/Resources/library"
    > print(sessionInfo())
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
    > print(version)
                   _                           
    platform       x86_64-apple-darwin15.6.0   
    arch           x86_64                      
    os             darwin15.6.0                
    system         x86_64, darwin15.6.0        
    status                                     
    major          3                           
    minor          6.1                         
    year           2019                        
    month          07                          
    day            05                          
    svn rev        76782                       
    language       R                           
    version.string R version 3.6.1 (2019-07-05)
    nickname       Action of the Toes
