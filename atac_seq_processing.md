ATAC-seq processing
================
Angela Simeone

### Adaptor trimming

``` bah
cd /atac_seq_folder/SLX-run1
for fq1 in *r_1.fq.gz
do
  fq2=${fq1%%r_1.fq.gz}r_2.fq.gz
  ls $fq1
  ls $fq2
  sbatch -o %j.out -e %j.err --mem 16000 --wrap "cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAG -o ${fq1%%.fq.gz}.trimmed.fq.gz -p ${fq2%%.fq.gz}.trimmed.fq.gz $fq1 $fq2 "
done
```

### Aligment

``` bash
# === aligning ===
# do this in both folders - 2 individual sequencing runs of the same set of libraries
cd /atac_seq_folder/SLX-run1

cd /atac_seq_folder/SLX-run1
g='/scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.fa'
w='/scratcha/sblab/simeon01/reference_genomes/hg38/hg38.whitelist.sorted.bed'

mkdir aligned
for fq1_trimmed in *r_1.trimmed.fq.gz
do
  fq2_trimmed=${fq1_trimmed%%r_1.trimmed.fq.gz}r_2.trimmed.fq.gz
  echo sbatch -o %j.out -e %j.err --mem 16G --wrap="bwa mem -M -t 12 $g $fq1_trimmed $fq2_trimmed | samtools view -Sb -F780 -q 10 -L $w - | samtools sort -@  > aligned/${fq1_trimmed%%_R1_001.trimmed.fq.gz}.hg38.sort.bam"
  sbatch -o %j.out -e %j.err --mem 16G --wrap="bwa mem -M -t 12 $g $fq1_trimmed $fq2_trimmed | samtools view -Sb -F780 -q 10 -L $w - | samtools sort -@ 8 - > aligned/${fq1_trimmed%%_R1_001.trimmed.fq.gz}.hg38.sort.bam"
done

cd aligned
for file in *bam
  do
  sbatch --mem 16G --wrap "samtools sort -@ 8 $file -o ${file/.bam/.sorted.bam}"
done
```

### Remove duplicates

``` bash

cd /atac_seq_folder/SLX-run1/aligned
for file in *.hg38.sort.bam
do
  sbatch --mem 16G  --wrap "java -Xmx7g -jar /home/simeon01/applications/picard-2.20.3.jar MarkDuplicatesWithMateCigar INPUT=$file OUTPUT=${file%%.bam}.markduplicates.bam REMOVE_DUPLICATES=true  AS=true METRICS_FILE=${file%%.bam}.markduplicates_metrics.txt"
  echo "==========="
  echo "          "
done
```

    ## bash: line 1: cd: /atac_seq_folder/SLX-run1/aligned: No such file or directory
    ## bash: line 4: sbatch: command not found
    ## ===========
    ## 

### generate stats

``` bash
cd /atac_seq_folder/SLX-run1/trimmed/aligned
for file in *hg38.sort.bam
do
bam_hg38=$file
#bam_f_nodup=${bam%%.hg38.sort.bam}.hg38.sort.markduplicates.bam
bam_hg38_nodup=${file%%.bam}.markduplicates.bam

ls -lth $bam
echo "====="
ls -lth $bam_hg38_nodup

#sbatch -o %j.out.stat --mem 8G --wrap "samtools view -c -F 260 $bam_hg38 >  ${file%%.bam}.stat2"
#sbatch -o %j.out.stat --mem 8G --wrap "samtools view -c -F 260 $bam_hg38_nodup >  ${file%%.bam}.stat5"

sbatch -o %j.out.stat --mem 1G --wrap "samtools idxstats $bam_hg38_nodup | grep chrM | cut -f 3 > ${file%%.bam}.statchrMtot"

done

cd /atac_seq_folder/SLX-run1/trimmed/aligned

touch pyPDS_hypoxia_atac_SLX-run1.out.stat
for bam in *.sort.markduplicates.bam
do
base_name=${bam%%.markduplicates.bam}
#base_name=${bam/.sort.markduplicates.bam/.sort.bam}
files_stat=`ls *.stat* | grep $base_name`
echo $bam > temp
cat $files_stat >> temp
paste -d '\t' pyPDS_hypoxia_atac_SLX-run1.out.stat temp >> temp2
mv temp2 pyPDS_hypoxia_atac_SLX-run1.out.stat
done
```

### Index bam

``` bash
cd /atac_seq_folder/SLX-run1/trimmed/aligned
for bam in *.sort.markduplicates.bam
do
sbatch -o %j.out.log --mem 16G --wrap "samtools index $bam"
done
```

### Estimate insert size distribution on the deduplicated bams (insert\_size\_estimation)

Picard CollectInsertSizeMetrics has been used to estimate insert size

``` bash

for file in *.markduplicates.bam
do
sbatch --mem 12G --wrap "samtools view $file | cut -f 9 > ${file%%.markduplicates.bam}.insert_size.txt"
done



for file in *markdupl*bam
do
  #/Users/marsic01/applications/picard-2.8.2.jar MarkDuplicates ==> this was the one used before
  # now picard in home
  #picard MarkDuplicates --version
  #2.18.12-SNAPSHOT
  sbatch --mem 16G  --wrap "java -Xmx7g -jar /home/simeon01/applications/picard-2.20.3.jar CollectInsertSizeMetrics INPUT=$file OUTPUT=${file%%.markduplicates.bam}.insert_size_metrics.txt H=${file%%.markduplicates.bam}.insert_size_metrics.pdf"
  echo "==========="
  echo "          "
done
```

### Generate tracks for human

``` bash

cd /atac_seq_folder/SLX-run1/trimmed/aligned
for bam in *.sort.markduplicates.bam
do
tot_r_hg38=`cat ${bam%%.markduplicates.bam}.stat5`
scal_factor_hg38=$(awk -v m=$tot_r_hg38 'BEGIN { print 1000000/m }')
echo $scal_factor_hg38
echo sbatch -o %j.out.log --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"
sbatch -o %j.out.log --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"

# in the following, step size is larger: useful to get files smaller in size but there is a loss in resolution. 
#sbatch -o %j.out.log --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 50 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w50.rpm.bw"
#sbatch -o %j.out.log --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 100 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w100.rpm.bw"
done

for bam in *markduplicates.bam
    do
        tot_r_hg38=`cat ${bam%%.markduplicates.bam}.stat5`
        scal_factor_hg38=$(awk -v m=$tot_r_hg38 'BEGIN { print 1000000/m }')
        echo $scal_factor_hg38
      echo sbatch --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"
        sbatch --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"
        echo "=============="
done
```

Call macs2 peaks with default options
-------------------------------------

``` bash
cd /atac_seq_folder/SLX-run1/aligned

out_dir_narrow='/atac_seq_folder/SLX-run1/aligned/macs2_individual_rep'
# === call peaks on atac-seq libraries where Chromosome M is excluded
mkdir $out_dir_narrow
for file in *.markduplicates.bam
  do 
  tmpBam=${file%%.markduplicates.bam}.markduplicates.temp
  bname=${file%%.markduplicates.bam}.macs2_peaks
  sbatch --mem 16G --wrap="samtools view -h -F1024 $file | grep -v -P '\tchrM\t' | samtools view -b - > $tmpBam && macs2 callpeak --keep-dup all -t $tmpBam -n $out_dir_narrow/${bname}"
done
```

### Generate consesus for each experimental condition

``` bash

cd /atac_seq_folder/SLX-run1/aligned/macs2_individual_rep
mkdir consensus_peaks

conditions=(1uM_pyPDS_2hr_hypoxia_whole_cell_ATACseq 1uM_pyPDS_2hr_normoxia_whole_cell_ATACseq DMSO_control_2hr_hypoxia_whole_cell_ATACseq DMSO_control_2hr_normoxia_whole_cell_ATACseq)

for con in ${conditions[@]}
do
multiIntersectBed -i ${con}*narrowPeak | awk '{if ($4>=2) print $0}' | sortBed -i - | mergeBed -i - >./consensus_peaks/${con}.multi2.bed
done

multiIntersectBed -i hypoxia_1hr_wholeK562_fragment_1hr*narrowPeak | awk '{if ($4>=2) print $0}' | sortBed -i - | mergeBed -i - > hypoxia_1hr_wholeK562_fragment_1hr.multi2.bed

multiIntersectBed -i no_garcinol_4hr_wholeK562_fragment*narrowPeak | awk '{if ($4>=2) print $0}' | sortBed -i - | mergeBed -i - > no_garcinol_4hr_wholeK562_fragment.multi2.bed

multiIntersectBed -i normoxia_1hr_wholeK562_fragment*narrowPeak | awk '{if ($4>=2) print $0}' | sortBed -i - | mergeBed -i - > normoxia_1hr_wholeK562_fragment.multi2.bed
```

Narrow Peaks consensus
======================

``` bash
out_consensum=/atac_seq_folder/SLX-run1/trimmed/aligned/macs2_individual_rep/consensus_narrow_peaks
condition=(1uM_pyPDS_2hr_hypoxia DMSO_control_2hr_hypoxia)
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
```

### merge two sequencing runs

``` bash
  # === merge the 2 sequencing runs ===
  
  path_run1=/atac_seq_folder/SLX-run1/aligned
  path_run2=/atac_seq_folder/SLX-run2/aligned
  
  merged_runs=/atac_seq_folder/merged_atac_SLX-run1_SLX-run2
  mkdir $merged_runs
  
  cd $path_run2
  
list_experiments=`ls *hg38.sort.bam`
  
for exp in ${list_experiments[@]}
do
  base_name=${exp%%_SLX*}
  echo $base_name
  bam1=`ls $path_run1/$base_name*.hg38.sort.bam`
  bam2=`ls $path_run2/$base_name*.hg38.sort.bam`
  ls -lth $bam1
  ls -lth $bam2
  out_bam=$merged_runs/$base_name.SLX-run1_SLX-run2.merged.bam
  cmd="samtools merge -@ 4 $out_bam $bam1 $bam2"
  sbatch --mem 16G --wrap "$cmd"
  echo $out_bam
  echo $cmd
  echo " ============== "
done
  
cd /atac_seq_folder/merged_atac_SLX-run1_SLX-run2
for file in *bam
do
  sbatch --mem 16G --wrap "samtools sort -@ 8 $file -o ${file/.bam/.sorted.bam}"
done
  # === remove duplicates ====
cd /atac_seq_folder/merged_atac_SLX-run1_SLX-run2
for file in *.merged.sorted.bam 
  do
  sbatch --mem 16G  --wrap "java -Xmx7g -jar /home/simeon01/applications/picard-2.20.3.jar MarkDuplicatesWithMateCigar INPUT=$file OUTPUT=${file%%.bam}.markduplicates.bam REMOVE_DUPLICATES=true  AS=true METRICS_FILE=${file%%.bam}.markduplicates_metrics.txt"
  echo "==========="
  echo "          "
done
  
  
  ## index markduplicates bam
  
```

\`\`\`bash cd /atac\_seq\_folder/SLX-run2/trimmed/aligned for bam in \*.sorted.markduplicates.bam do sbatch -o %j.out.log --mem 16G --wrap "samtools index $bam" done

\`\`\`

\#\# generate stats

\`\`\`bash cd /atac\_seq\_folder/merged\_atac\_SLX-run1\_SLX-run2 for file in \*.merged.sorted.bam do bam\_hg38=$file  \#bam\_f\_nodup=${bam%%.hg38.sort.bam}.hg38.sort.markduplicates.bam bam\_hg38\_nodup=${file%%.bam}.markduplicates.bam

\#ls -lth $bam\_hg38 echo "=====" ls -lth $bam\_hg38\_nodup

sbatch -o %j.out.stat --mem 8G --wrap "samtools view -c -F 260 $bam\_hg38 &gt; ${file%%.bam}.stat2" sbatch -o %j.out.stat --mem 8G --wrap "samtools view -c -F 260 $bam\_hg38\_nodup &gt; ${file%%.bam}.stat5"

sbatch -o %j.out.stat --mem 1G --wrap "samtools idxstats $bam\_hg38\_nodup | grep chrM | cut -f 3 &gt; ${file%%.bam}.statchrMtot"

done

cd /atac\_seq\_folder/merged\_atac\_SLX-run1\_SLX-run2

touch pyPDS\_hypoxia\_merged\_atac\_SLX-run1\_SLX-run2.out.stat for bam in \*.merged.sorted.bam do base\_name=${bam%%.bam}  \#base\_name=${bam/.sort.markduplicates.bam/.sort.bam} files\_stat=`ls *.stat* | grep $base_name` echo $bam &gt; temp cat $files\_stat &gt;&gt; temp paste -d '' pyPDS\_hypoxia\_merged\_atac\_SLX-run1\_SLX-run2.out.stat temp &gt;&gt; temp2 mv temp2 pyPDS\_hypoxia\_merged\_atac\_SLX-run1\_SLX-run2.out.stat done \`\`\`

### insert\_size\_estimation after combining 2 sequencing runs

I used picard CollectInsertSizeMetrics

``` bash

for file in *.markduplicates.bam
do
sbatch --mem 12G --wrap "samtools view $file | cut -f 9 > ${file%%.markduplicates.bam}.insert_size.txt"
done

for file in *markdupl*bam
do
  #/Users/marsic01/applications/picard-2.8.2.jar MarkDuplicates ==> this was the one used before
  # now picard in home
  #picard MarkDuplicates --version
  #2.18.12-SNAPSHOT
  sbatch --mem 16G  --wrap "java -Xmx7g -jar /home/simeon01/applications/picard-2.20.3.jar CollectInsertSizeMetrics INPUT=$file OUTPUT=${file%%.markduplicates.bam}.insert_size_metrics.txt H=${file%%.markduplicates.bam}.insert_size_metrics.pdf"
  echo "==========="
  echo "          "
done
```

### generate tracks for human after combining sequencing runs

``` bash

cd /atac_seq_folder/SLX-run2/trimmed/aligned
for bam in *.markduplicates.bam
do
tot_r_hg38=`cat ${bam%%.markduplicates.bam}.stat5`
scal_factor_hg38=$(awk -v m=$tot_r_hg38 'BEGIN { print 1000000/m }')
echo $scal_factor_hg38
echo sbatch -o %j.out.log --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"
#sbatch -o %j.out.log --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"
#sbatch -o %j.out.log --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 50 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w50.rpm.bw"
sbatch -o %j.out.log --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 100 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w100.rpm.bw"
done
```

### Call macs2 peaks with default options and also broad option

``` bash
cd /atac_seq_folder/merged_atac_SLX-run1_SLX-run2

out_dir_narrow='/atac_seq_folder/merged_atac_SLX-run1_SLX-run2/macs2_individual_rep'
# === call peaks on atac-seq libraries where Chromosome M is excluded
mkdir $out_dir_narrow
for file in *.markduplicates.bam
do 
tmpBam=${file%%.markduplicates.bam}.markduplicates.temp
bname=${file%%.markduplicates.bam}.macs2_peaks
sbatch --mem 16G --wrap="samtools view -h -F1024 $file | grep -v -P '\tchrM\t' | samtools view -b - > $tmpBam && macs2 callpeak --keep-dup all -t $tmpBam -n $out_dir_narrow/${bname}"
done



for bam in *markduplicates.bam
do
tot_r_hg38=`cat ${bam%%.markduplicates.bam}.stat5`
scal_factor_hg38=$(awk -v m=$tot_r_hg38 'BEGIN { print 1000000/m }')
echo $scal_factor_hg38
echo sbatch --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"
sbatch --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"
echo "=============="
done
```

### Generate consesus for each experimental condition

``` bash

cd /atac_seq_folder/merged_atac_SLX-run1_SLX-run2/macs2_individual_rep
mkdir consensus_peaks

conditions=(1uM_pyPDS_2hr_hypoxia_whole_cell_ATACseq 1uM_pyPDS_2hr_normoxia_whole_cell_ATACseq DMSO_control_2hr_hypoxia_whole_cell_ATACseq DMSO_control_2hr_normoxia_whole_cell_ATACseq)

for con in ${conditions[@]}
do
multiIntersectBed -i ${con}*narrowPeak | awk '{if ($4>=2) print $0}' | sortBed -i - | mergeBed -i - >./consensus_peaks/${con}.multi2.bed
done
 
```
