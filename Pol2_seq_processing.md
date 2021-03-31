Pol2-chip-seq
================
Angela Simeone

The data processing pipeline involved different steps: fastQC, adaptor trimmin, alignment, deduplication and identification of regions with local enrichments (peaks). In the following, those steps are delineated. Pol2-ChIP-seq have been conducted over 5 independent biological replicates therefore they have been processed independently. To identify the regions representing the consensus, peak files have been compared and only regions observed in 3 over 5 replicates have been selected and included in the consensus (i.e. regions consistently observed in 60% of the cases).

### adaptor trimming

``` bah
cd /local/Pol2ChIP/SLX-18974
mkdir trimmed
for fq in *.fq.gz
do
  bname=`basename $fq`
  sbatch -o %j.$bname.tmp.out -e %j.$bname.tmp.err --mem 16000 --wrap "cutadapt -q 20 -O 3 -a CTGTCTCTTATACACATCT -o trimmed/${bname%%.fq.gz}.trimmed.fq.gz $fq"
done
```

### aligment

``` bash

path=/local/Pol2ChIP/SLX-18974/trimmed
g='~/reference_genomes/hg38/hg38_selected.fa'
w='~/reference_genomes/hg38/hg38.whitelist.sorted.bed'
mkdir aligned
for f in *gz
do
    sbatch -o %j.$f.tmp.out -e %j.$f.tmp.err --mem 16000 --wrap "bwa mem -M $g $f | samtools view -S -u -F2304 -q 10 -L $w - | samtools sort -@ 8 - > $path/aligned/${f%%.trimmed.fq.gz}.hg38.sort.bam"
done
  
```

remove duplicates
=================

``` bash
cd /local/Pol2ChIP/SLX-18974/trimmed/aligned
for bam in *hg38.sort.bam
do
  sbatch -o %j.$bam.tmp.out -e %j.$bam.tmp.err --mem 32000  --wrap "java -Xmx7g -jar ~/applications/picard-2.20.3.jar MarkDuplicates INPUT=$bam OUTPUT=${bam%%.bam}.markduplicates.bam REMOVE_DUPLICATES=true  AS=true METRICS_FILE=${bam%%.bam}.markduplicates_metrics.txt"
done
```

### generate stats

In turn 2 bam files are checked: bam obtained after alignment and bam obtained after deduplications. The total number of reads of each one are counted and stored in `stat2` and `stat5` respectively.

``` bash
for file in *hg38.sort.bam
  do
  bam_hg38=$file
  #bam_f_nodup=${bam%%.hg38.sort.bam}.hg38.sort.markduplicates.bam
  bam_hg38_nodup=${file%%.bam}.markduplicates.bam
  
  ls -lth $bam
  echo "====="
  ls -lth $bam_hg38_nodup
  
  sbatch -o %j.out.stat --mem 8G --wrap "samtools view -c -F 260 $bam_hg38 >  ${file%%.bam}.stat2"
  sbatch -o %j.out.stat --mem 8G --wrap "samtools view -c -F 260 $bam_hg38_nodup >  ${file%%.bam}.stat5"
done
```

Collect stats and print them to file.

``` batch
cd /local/Pol2ChIP/SLX-18974/trimmed/aligned

touch Pol2_hypoxia_SLX-18974.out.stat
for bam in *.sort.markduplicates.bam
  do
  base_name=${bam%%.markduplicates.bam}
  #base_name=${bam/.sort.markduplicates.bam/.sort.bam}
  files_stat=`ls *.stat* | grep $base_name`
  echo $bam > temp
  cat $files_stat >> temp
  paste -d '\t' Pol2_hypoxia_SLX-18974.out.stat temp >> temp2
  mv temp2 Pol2_hypoxia_SLX-18974.out.stat
done
```

### generate tracks for human

``` bash
#index bam file
cd /local/Pol2ChIP/SLX-18974/trimmed/aligned
for bam in *.sort.markduplicates.bam
    do
  sbatch -o %j.out.log --mem 16G --wrap "samtools index $bam"
done

# generate track files 
for bam in *.sort.markduplicates.bam
    do
        tot_r_hg38=`cat ${bam%%.markduplicates.bam}.stat5`
        scal_factor_hg38=$(awk -v m=$tot_r_hg38 'BEGIN { print 1000000/m }')
        echo $scal_factor_hg38
      echo sbatch -o %j.out.log --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"
        sbatch -o %j.out.log --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"
done
```

### call macs2 peaks with default options

``` bash
cd /local/Pol2ChIP/SLX-18974/trimmed/aligned

rep=(1 2 3 4 5)
condition=(hypoxia_1hr_Pol2_ normoxia_1hr_Pol2_)
counter_l=1
out_dir_narrow='/local/Pol2ChIP/SLX-18974/trimmed/aligned/macs2_individual_rep'
out_dir_broad='/local/Pol2ChIP/SLX-18974/trimmed/aligned/macs2_broad_individual_rep'
mkdir $out_dir_narrow
mkdir $out_dir_broad
for currcon in ${condition[@]};
do
  for i in ${rep[@]}; 
  do
      echo "${c}"
      echo "$i"
      #echo ${c}*ChIP_biorep$i*.sort.markduplicates.bam
      t=`ls ${currcon}*ChIP_biorep$i*.sort.markduplicates.bam`
      c=`ls ${currcon}*Input_biorep$i*.sort.markduplicates.bam`
      
      ls -lht $t
      ls -lht $c
      sbatch --mem 8000 -J macs --wrap "macs2 callpeak --keep-dup all -t $t -c $c -n $out_dir_narrow/${currcon}_ChIP_biorep${i}" # human is default
     
      echo "====="
  done    
done
```

### create consesnsum beteween the 2 different conditions

``` bash
mkdir 1e-3
mv *1e-3_peaks.narrowPeak ./1e-3/


# Narrow Peaks
out_consensum=/local/Pol2ChIP/SLX-18974/trimmed/aligned/macs2_individual_rep/consensus_narrow_peaks
n_thr=(2 3 4 5)
mkdir $out_consensum
for n in ${n_thr[@]}
do
  multiIntersectBed -i hypoxia*_peaks.narrowPeak | awk -v var="$n" '{if($4>=var) print $0}' | sortBed -i - | mergeBed -i - > $out_consensum/hypoxia.multi$n.bed
  multiIntersectBed -i normoxia*_peaks.narrowPeak | awk -v var="$n" '{if($4>=var) print $0}' | sortBed -i - | mergeBed -i - > $out_consensum/normoxia.multi$n.bed
done

for n in ${n_thr[@]}
do
  A=`ls hypoxia*multi${n}.bed`
  B=`ls normoxia*multi${n}.bed`
  common=`intersectBed -a $A -b $B |sortBed -i - | mergeBed -i - | wc -l`
  A_specific=`intersectBed -a $A -b $B -wa -v |sortBed -i - | mergeBed -i - | wc -l`
  B_specific=`intersectBed -b $A -a $B -wa -v|sortBed -i - | mergeBed -i -  | wc -l`
  tot_A=`wc -l $A | awk '{print $1}'`
  tot_B=`wc -l $B | awk '{print $1}'`
  
  echo  multi$n$'\t'$tot_A$'\t'$tot_B$'\t'$common$'\t'$A_specific$'\t'$B_specific
  #echo  multi$n$'\t'$tot_A$'\t'$tot_B$'\t'$common$'\t'$A_specific$'\t'$B_specific >> stats_consensus_peaks.txt
done
```
