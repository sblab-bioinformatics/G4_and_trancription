Processing hypoxia in U2OS cells
================
Angela Simeone
15/02/2021

Processing of data files
------------------------

Basic processing involved:

-   check base quality (fastqc)
-   adaptor trimming
-   aligment
-   merge bams of the same library
-   duplicate removal
-   stats generation
-   track files (\*bw) generation
-   peak calling
-   bio rep consensus region extraction
-   consensus set for DBA
-   coverages across consensus set for DBA
-   differential binding analysis (DBA)

Quality control fastqc and adaptor trimming trimming

``` bash

for file in *.fq.gz
  do
  
  # fastqc
  sbatch --time 02:00:00 --mem 4G --wrap "fastqc --noextract --nogroup -q -o fastqc/ $file"
  
  #adaptor trimming
  bname=`basename $file`
  sbatch --time 02:00:00 --mem 12000 --wrap "cutadapt -q 20 -O 3 -a CTGTCTCTTATACACATCT -o trimmed/${bname%%.fq.gz}.trimmed.fq.gz $file"
done
```

Alignment (bwa-mem)

``` bash
g='~/reference_genomes/hg38/hg38_selected.fa'
w='~/reference_genomes/hg38/hg38.whitelist.sorted.bed'
seq_runs=(SLX-20346 SLX-20347)

for run in ${seq_runs[@]}
do
  path=~/G4_chipseq_u2os/$run/trimmed
  cd $path

  mkdir aligned
  for f in *gz
  do
  sbatch --time 12:00:00 --mem 16G --wrap "bwa mem -M $g $f | samtools view -S -u -F2304 -q 10 -L $w - | samtools sort -@ 8 - > $path/${f%%.trimmed.fq.gz}.hg38.sort.bam && java -Xmx7g -jar ~/applications/picard-2.20.3.jar MarkDuplicates INPUT=$path/${f%%.trimmed.fq.gz}.hg38.sort.bam OUTPUT=$path/${f%%.trimmed.fq.gz}.hg38.sort.markduplicates.bam REMOVE_DUPLICATES=true  AS=true METRICS_FILE=$path/${f%%.trimmed.fq.gz}.markduplicates_metrics.txt"
  done
done
```

Merge separate runs, identify duplicates and remove them

``` bash
 # === merge the 2 sequencing runs ===
path_run1=~/G4_chipseq_u2os/SLX-20346
path_run2=~/G4_chipseq_u2os/SLX-20347
merged_runs=~/G4_chipseq_u2os/merged_BG4_SLX-20346_SLX-20347

mkdir $merged_runs
cd $path_run2/trimmed

list_experiments=`ls *hg38.sort.bam`
  
for exp in ${list_experiments[@]}
do
  base_name=${exp%%_SLX*}
  echo $base_name
  bam1=`ls $path_run1/trimmed/$base_name*.hg38.sort.bam`
  bam2=`ls $path_run2/trimmed/$base_name*.hg38.sort.bam`
  ls -lth $bam1
  ls -lth $bam2
  out_bam=$merged_runs/$base_name.SLX-20346_SLX-20347.merged.bam
  cmd="samtools merge -@ 4 $out_bam $bam1 $bam2"
  sbatch --time 04:00:00 --mem 16G --wrap "$cmd"
  echo $out_bam
  echo $cmd
  echo " ============== "
done
  
cd ~/G4_chipseq_u2os/merged_BG4_SLX-20346_SLX-20347
seq_runs=(SLX-20346 SLX-20347)
for file in *bam
do
  sbatch --mem 16G --wrap "samtools sort -@ 8 $file -o ${file/.bam/.sorted.bam}"
done

# === remove duplicates ====
cd ~/G4_chipseq_u2os/merged_BG4_SLX-20346_SLX-20347
for file in *merged.sorted.bam
  do
  sbatch --mem 16G --time 02:00:00 --wrap "java -Xmx7g -jar ~/applications/picard-2.20.3.jar MarkDuplicates INPUT=$file OUTPUT=${file%%.bam}.markduplicates.bam REMOVE_DUPLICATES=true  AS=true METRICS_FILE=${file%%.bam}.markduplicates_metrics.txt"
  echo "==========="
  echo "          "
done

# indexing bams
cd ~/G4_chipseq_u2os/merged_BG4_SLX-20346_SLX-20347
for bam in *.markduplicates.bam
  do
  sbatch --time 00:10:00 --mem 16G --wrap "samtools index $bam"
done
  
## generate stats  (stat2: tot reads; stat5: tot reads no dup)
cd ~/G4_chipseq_u2os/merged_BG4_SLX-20346_SLX-20347
for file in *sorted.bam
  do
  bam_hg38=$file
  #bam_f_nodup=${bam%%.hg38.sort.bam}.hg38.sort.markduplicates.bam
  bam_hg38_nodup=${file%%.bam}.markduplicates.bam
  
  ls -lth $bam
  echo "====="
  ls -lth $bam_hg38_nodup
  
  sbatch --time 00:05:00 --mem 8G --wrap "samtools view -c -F 260 $bam_hg38 >  ${file%%.bam}.stat2"
  sbatch --time 00:05:00 --mem 8G --wrap "samtools view -c -F 260 $bam_hg38_nodup >  ${file%%.bam}.stat5"
done

cd ~/G4_chipseq_u2os/merged_BG4_SLX-20346_SLX-20347

touch U2OS_hypoxia_BG4_SLX-20346_SLX-20347.out.stat
for bam in *.markduplicates.bam
do
base_name=${bam%%.markduplicates.bam}
#base_name=${bam/.sort.markduplicates.bam/.sort.bam}
files_stat=`ls *.stat* | grep $base_name`
echo $bam > temp
cat $files_stat >> temp
paste -d '\t' U2OS_hypoxia_BG4_SLX-20346_SLX-20347.out.stat temp >> temp2
mv temp2 U2OS_hypoxia_BG4_SLX-20346_SLX-20347.out.stat
done
```

### generate tracks for human

``` bash

cd ~/G4_chipseq_u2os/merged_BG4_SLX-20346_SLX-20347
for bam in *.sorted.markduplicates.bam
do
tot_r_hg38=`cat ${bam%%.markduplicates.bam}.stat5`
scal_factor_hg38=$(awk -v m=$tot_r_hg38 'BEGIN { print 1000000/m }')
echo $scal_factor_hg38
echo sbatch --time 06:00:00 --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"
sbatch --time 06:00:00 --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"
```

### call peaks

Peak calling: macs2 with default options

``` bash
cd ~/G4_chipseq_u2os/merged_BG4_SLX-20346_SLX-20347

rep=(1 2 3)
condition=(hypoxia_1hr_U2OS_biorep1 hypoxia_1hr_U2OS_biorep3 hypoxia_1hr_U2OS_biorep4 normoxia_1hr_U2OS_biorep1 normoxia_1hr_U2OS_biorep3 normoxia_1hr_U2OS_biorep4)
counter_l=1
out_dir_narrow='~/G4_chipseq_u2os/merged_BG4_SLX-20346_SLX-20347/macs2_individual_rep'
out_dir_broad='~/G4_chipseq_u2os/merged_BG4_SLX-20346_SLX-20347/macs2_broad_individual_rep'
mkdir $out_dir_narrow
mkdir $out_dir_broad
cd ~/G4_chipseq_u2os/merged_BG4_SLX-20346_SLX-20347
for currcon in ${condition[@]};
do
  for i in ${rep[@]}; 
  do
      echo "${c}"
      echo "$i"
      
      t=`ls ${currcon}*ChIP$i*.markduplicates.bam`
      c=`ls ${currcon}*Input*.markduplicates.bam`
      
      ls -lht $t
      ls -lht $c
      sbatch --time 02:00:00 --mem 8000 -J macs --wrap "macs2 callpeak --keep-dup all -t $t -c $c -n $out_dir_narrow/${currcon}_ChIP_biorep${i}" # human is default
      echo "====="
  done    
done
```

Narrow Peaks consensus

``` bash
out_consensum=~/G4_chipseq_u2os/merged_BG4_SLX-20346_SLX-20347/macs2_individual_rep/consensus_narrow_peaks
condition=(hypoxia_1hr_U2OS_biorep1 hypoxia_1hr_U2OS_biorep3 hypoxia_1hr_U2OS_biorep4 normoxia_1hr_U2OS_biorep1 normoxia_1hr_U2OS_biorep3 normoxia_1hr_U2OS_biorep4)
#n_thr=(2 3)
mkdir $out_consensum

for case in ${condition[@]}
do
echo $case
bed_files=`ls *$case*.narrowPeak`

echo $bed_files

multiIntersectBed -i $bed_files | awk '{if($4>=2) print $0}' | sortBed -i | mergeBed -i - > $out_consensum/$case.multi2.bed
echo " ================ . ================"
done

wc -l *bed | awk '{print $2"\t"$1}' 
#hypoxia_1hr_U2OS_biorep1.multi2.bed    9163
#hypoxia_1hr_U2OS_biorep3.multi2.bed    10721
#hypoxia_1hr_U2OS_biorep4.multi2.bed    12048
#normoxia_1hr_U2OS_biorep1.multi2.bed   17567
#normoxia_1hr_U2OS_biorep4.multi2.bed   21098
#normoxia_1hr_U2OS_biorep3.multi2.bed   11782
```

Overlap with OQs

``` bash
oqs=~/reference_genomes/OQs/OQ_hits.lifted_hg19_to_hg38_no_.bed


cd ~/G4_chipseq_u2os/merged_BG4_SLX-20346_SLX-20347/macs2_individual_rep
for file in *narrowPeak
do
n_peaks=`wc -l $file | awk '{print $1}'`
n_peaks_in_oqs=`intersectBed -a $file -b $oqs -wa | sort | uniq | wc -l`

echo  $file$'\t'$n_peaks$'\t'$n_peaks_in_oqs
echo  $file$'\t'$n_peaks$'\t'$n_peaks_in_oqs >> N_peaks_U2OS_hypoxia_BG4_SLX-20346_SLX-20347.individualLib.txt
done

cd ~/G4_chipseq_u2os/merged_BG4_SLX-20346_SLX-20347/macs2_individual_rep/consensus_narrow_peaks
for file in *2.bed
do
n_peaks=`wc -l $file | awk '{print $1}'`
n_peaks_in_oqs=`intersectBed -a $file -b $oqs -wa | sort | uniq | wc -l`

echo  $file$'\t'$n_peaks$'\t'$n_peaks_in_oqs
echo  $file$'\t'$n_peaks$'\t'$n_peaks_in_oqs >> N_peaks_U2OS_hypoxia_BG4_SLX-20346_SLX-20347.multi2.bio2.txt
done
```

Overlap between biological rep

``` bash
out_consensum=~/G4_chipseq_u2os/merged_BG4_SLX-20346_SLX-20347/macs2_individual_rep/consensus_narrow_peaks
cond=(normoxia hypoxia) 
cd ~/G4_chipseq_u2os/merged_BG4_SLX-20346_SLX-20347/macs2_individual_rep/consensus_narrow_peaks
for case in ${cond[@]}
do
echo $case
bed_files=`ls *$case*.multi2.bed`

echo $bed_files

multiIntersectBed -i $bed_files | awk '{if($4>=2) print $0}' | sortBed -i | mergeBed -i - > $out_consensum/$case.bio2.bed
echo " ================ . ================"
done
```

Coverages for differential binding analysis analysis
----------------------------------------------------

We considered for the analysis all regions obtained as the merge of all `*bio2.bed`. For the analysis, libraries have been normalised to the total sequencing depth (sequencing detphs have been imported).

``` bash
cd ~/G4_chipseq_u2os/merged_BG4_SLX-20346_SLX-20347/macs2_individual_rep/consensus_narrow_peaks

ref2_u2os=u2osBG4.SLX-20346_SLX-20347.merge_bio2.bed

cat *bio2.bed | sortBed -i - | mergeBed -i - > $ref2_u2os


#create folders for coverage
bed_consensus2=~/G4_chipseq_u2os/merged_BG4_SLX-20346_SLX-20347/macs2_individual_rep/consensus_narrow_peaks/u2osBG4.SLX-20346_SLX-20347.merge_bio2.bed

path_out=~/G4_chipseq_u2os/merged_BG4_SLX-20346_SLX-20347/macs2_individual_rep/consensus_narrow_peaks/u2os_hypoxia_normoxia_cov
mkdir $path_out

path_bam0=~/G4_chipseq_u2os/merged_BG4_SLX-20346_SLX-20347
all_bam_paths=($path_bam0)

cd $path_out

for bam in *markduplicates.bam
  do
  cmd_coverage2="bedtools coverage -a $bed_consensus2 -b ${bam} -counts > $path_out/${bam%%.bam}.dba_consensus.merge_bio2.coverage.bed"
  echo $cmd_coverage2
  sbatch --time 01:00:00 --mem 4G --wrap "$cmd_coverage2"
done
```

copy data locally to run differential analysis

``` bash

rm(list = ls())
path1='/Users/simeon01/Documents/Karen/G4_chipseq_u2os/u2os_hypoxia_normoxia_cov'
pattern_file_extension="dba_consensus.merge_bio2.coverage.bed$"
selection_files_pattern="SLX-20346_SLX-20347"
pattern_folder_prefix="DBA_G4_in_G4bio2_filterOverInput"
case=1
flag_BG4_at_TSS=0
input_flag=1
threshold_over_input=2
file_subset_TSS='/Users/simeon01/Documents/Karen/20201104_karen_hypoxia_normoxia_pol2_atac_bg4/overlap_of_merge_BG4_POL2_normo_hypo.at_TSS_plus_minus_500bp.bed'
source('/Users/simeon01/Documents/Karen/20201104_karen_hypoxia_normoxia_pol2_atac_bg4/DBA_at_G4_pol2_tss_plus_minus500.R')

rm(list = ls())
path1='/Users/simeon01/Documents/Karen/G4_chipseq_u2os/u2os_hypoxia_normoxia_cov'
pattern_file_extension="dba_consensus.merge_bio2.coverage.bed$"
selection_files_pattern="SLX-20346_SLX-20347"
pattern_folder_prefix="DBA_G4_in_G4bio2"
case=1
flag_BG4_at_TSS=0
input_flag=0
threshold_over_input=2
file_subset_TSS='/Users/simeon01/Documents/Karen/20201104_karen_hypoxia_normoxia_pol2_atac_bg4/overlap_of_merge_BG4_POL2_normo_hypo.at_TSS_plus_minus_500bp.bed'
source('/Users/simeon01/Documents/Karen/20201104_karen_hypoxia_normoxia_pol2_atac_bg4/DBA_at_G4_pol2_tss_plus_minus500.R')
```
