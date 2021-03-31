Fold enrichment at genomic regions
================
Angela Simeone

Fold enrichments have been assessed using the [Genomic Association Tester (GAT)](https://gat.readthedocs.io/en/latest/contents.html).

``` bash

g4_normo=~/normoxia.bio2.bed

g4_hypo=~/hypoxia.bio2.bed

pol2_normo=~/normoxia.multi3.bed

pol2_hypo=~/hypoxia.multi3.bed

out_dir=/G4_pol2_dir

mkdir $out_dir

# copy files into directory
cp $g4_normo $out_dir
cp $g4_hypo $out_dir
cp $pol2_normo $out_dir
cp $pol2_hypo $out_dir
cd $out_dir

# check annotation at promoters for BG4 normo and hypo
mkdir FE_normo_hypo
workspace_gen=~/reference_genomes/hg38/hg38.whitelist.bed
genomic_annotation=~/reference_genomes/hg38/annotations_genomic_regions_gencode_v28.anno2.bed
#genomic_annotation_promoters=~/reference_genomes/hg38/genecode.v28.various_definitions_promoters.bed
genomic_annotation_promoters=~/reference_genomes/hg38/binned_prom/binned_prom.bed

cat ~/reference_genomes/hg38/annotations_genomic_regions_gencode_v28.anno2.bed  ~/reference_genomes/hg38/genecode.v28.promoters_plus_minus500.bed > ~/reference_genomes/hg38/annotations_genomic_regions_gencode_v28.anno2.and_prom_tss_plus_minus500.bed

genomic_annotation_promoters_updated_prom=~/reference_genomes/hg38/annotations_genomic_regions_gencode_v28.anno2.and_prom_tss_plus_minus500.bed

for file in *bio[23].bed
#for file in *multi3.bed
do
  echo $file
  
  base_annotation=`basename $genomic_annotation`
  cmd_gat="gat-run.py --ignore-segment-tracks --workspace=${workspace_gen}  --annotations=$genomic_annotation --segments=${file} -S ./FE_normo_hypo/${file%%.bed}_intervals_${base_annotation%%.bed}.dat -t 4"
  echo ${file%%.bed}_intervals_${base_annotation%%.bed}.dat
  echo $cmd_gat
  #sbatch --mem 6G --wrap "$cmd_gat"
  
  base_annotation_promoters=`basename $genomic_annotation_promoters`
  cmd_gat2="gat-run.py --ignore-segment-tracks --workspace=${workspace_gen}  --annotations=$genomic_annotation_promoters --segments=${file} -S ./FE_normo_hypo/${file%%.bed}_intervals_${base_annotation_promoters%%.bed}.dat -t 4"
  echo ${file%%.bed}_intervals_${base_annotation_promoters%%.bed}.dat
  echo $cmd_gat2
  #sbatch --mem 6G --wrap "$cmd_gat2"
  
  base_annotation_promoters_updated_prom=`basename $genomic_annotation_promoters_updated_prom`
  cmd_gat3="gat-run.py --ignore-segment-tracks --workspace=${workspace_gen}  --annotations=$genomic_annotation_promoters_updated_prom --segments=${file} -S ./FE_normo_hypo/${file%%.bed}_intervals_${base_annotation_promoters_updated_prom%%.bed}.dat -t 4"
  echo ${file%%.bed}_intervals_${base_annotation_promoters_updated_prom%%.bed}.dat
  echo $cmd_gat2
  sbatch --mem 6G --wrap "$cmd_gat3"
  
done

cd FE_normo_hypo
for file in *dat; do cut -f 1,2,3,8 $file > ${file%%.dat}.reduced.txt; done

#############################
## enrichments also for DRB and TPL
#############################
cd ~/Data/20201023_revise_TPL_DRB/DRB_bed
mkdir FE_DRB
workspace_gen=~/reference_genomes/hg38/hg38.whitelist.bed

genomic_annotation_promoters_updated_prom=~/reference_genomes/hg38/annotations_genomic_regions_gencode_v28.anno2.and_prom_tss_plus_minus500.bed

for file in *confirmed*bed
do

  base_annotation_promoters_updated_prom=`basename $genomic_annotation_promoters_updated_prom`
  cmd_gat3="gat-run.py --ignore-segment-tracks --workspace=${workspace_gen}  --annotations=$genomic_annotation_promoters_updated_prom --segments=${file} -S ./FE_DRB/${file%%.bed}_intervals_${base_annotation_promoters_updated_prom%%.bed}.dat -t 4"
  echo ${file%%.bed}_intervals_${base_annotation_promoters_updated_prom%%.bed}.dat
  echo $cmd_gat2
  sbatch --mem 6G --wrap "$cmd_gat3"
  
done
cd FE_DRB
for file in *dat; do cut -f 1,2,3,8 $file > ${file%%.dat}.reduced.txt; done


cd ~/Data/20201023_revise_TPL_DRB/TPL_bed
mkdir FE_TPL
workspace_gen=~/reference_genomes/hg38/hg38.whitelist.bed

genomic_annotation_promoters_updated_prom=~/reference_genomes/hg38/annotations_genomic_regions_gencode_v28.anno2.and_prom_tss_plus_minus500.bed

for file in *confirmed*bed
do

  base_annotation_promoters_updated_prom=`basename $genomic_annotation_promoters_updated_prom`
  cmd_gat3="gat-run.py --ignore-segment-tracks --workspace=${workspace_gen}  --annotations=$genomic_annotation_promoters_updated_prom --segments=${file} -S ./FE_TPL/${file%%.bed}_intervals_${base_annotation_promoters_updated_prom%%.bed}.dat -t 4"
  echo ${file%%.bed}_intervals_${base_annotation_promoters_updated_prom%%.bed}.dat
  echo $cmd_gat2
  sbatch --mem 6G --wrap "$cmd_gat3"
  
done
cd FE_TPL
for file in *dat; do cut -f 1,2,3,8 $file > ${file%%.dat}.reduced.txt; done
```
