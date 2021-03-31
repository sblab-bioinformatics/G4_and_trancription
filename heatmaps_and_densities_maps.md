Heatmaps and density plots
================
Angela Simeone

#### Heatmap and density plots at TSS

The following steps are necessary to generate the appropriate TSS (reference) files for the TSS plots.

``` bash
TSS_500=~/reference_genomes/hg38/gencode.v28.TSS_plus_minus_500bp.bed
TSS_plus_minus1000=~/reference_genomes/hg38/gencode.v28.TSS.slop1000.bed
g4_normo=~/Data/20201019_BG4hypoxia_normoxia_SLX-198891_SLX-198892/normoxia.bio2.bed
g4_hypo=~/Data/20201019_BG4hypoxia_normoxia_SLX-198891_SLX-198892/hypoxia.bio2.bed
pol2_normo=~/Data/20201019_BG4hypoxia_normoxia_SLX-198891_SLX-198892/normoxia.multi3.bed
pol2_hypo=~/Data/20201019_BG4hypoxia_normoxia_SLX-198891_SLX-198892/hypoxia.multi3.bed

# G4 overlapping TSS plus minus 500
intersectBed -a $g4_normo -b $TSS_500 -wa | sort | uniq > ${g4_normo/.bed/.at_TSS_plus_minus500.bed}

output_folder=~/Data/20201019_BG4hypoxia_normoxia_SLX-198891_SLX-198892/reference_bed_for_figure2_promoters_plus_minus500
mkdir $output_folder
cd ~/Data/20201019_BG4hypoxia_normoxia_SLX-198891_SLX-198892
G4_normoxia_in_TSS_plus_minus500=$output_folder/normoxia.bio2.at_TSS_plus_minus500.bed

# use TSS500
intersectBed -a $g4_normo -b $pol2_normo -wa| sortBed -i - |intersectBed -a $TSS_500 -b - -wa |sort | uniq | wc -l
#9893 

cd ~/Data/20201019_BG4hypoxia_normoxia_SLX-198891_SLX-198892
intersectBed -a $g4_normo -b $pol2_normo -wa|sort | uniq |intersectBed -a $TSS_500 -b - -wa |sort | uniq > $output_folder/Prom_with_G4_overlap_Pol2_Normo_in_TSSplus_minus500.bed
intersectBed -a $TSS_500 -b $pol2_normo -wa|sort | uniq | intersectBed -a - -b $output_folder/Prom_with_G4_overlap_Pol2_Normo_in_TSSplus_minus500.bed -wa -v > $output_folder/Prom_NoG4_Pol2_Normo_in_TSSplus_minus500.bed
intersectBed -a $TSS_500 -b $pol2_normo -wa -v | intersectBed -a - -b $output_folder/Prom_with_G4_overlap_Pol2_Normo_in_TSSplus_minus500.bed -v -wa| sort | uniq | intersectBed -a - -b  $output_folder/Prom_NoG4_Pol2_Normo_in_TSSplus_minus500.bed -wa -v | sort | uniq > $output_folder/Prom_NoG4_NoPol2_Normo_in_TSSplus_minus500.bed
cd $output_folder

## intersectBed -a $g4_normo -b $pol2_normo -wa| sortBed -i - |intersectBed -b $TSS_500 -a -  -wa|sort | uniq | wc -l ## this was used to select the G4 that overlaps the regions


for file in Prom*.bed
do
echo $file
cut -f 4 $file > ${file/.bed/.tmp}
LC_ALL=C grep -f ${file/.bed/.tmp} ~/reference_genomes/hg38/gencode.v28.genebody.bed > ${file/.bed/.genebody_ref.bed}
echo "=============="
done
```

Use deeptools to compute coverage signal from library normalised BW files.

``` bash

### perform computation of density signal (computeMatrix reference-point --> TSS)


#atac case
cd $bw_dir_atac
for file in  *ox*bw
do
sbatch --time 03:00:00 --mem 8G --wrap="computeMatrix reference-point -S ${file} -R $folder_ref_bed/*.genebody_ref.bed -a 1000 -b 1000 --outFileName $out_path1/${file}.matrix.gz --missingDataAsZero && 
  plotProfile --matrixFile $out_path1/${file}.matrix.gz -o $out_path1/${file}.matrix.pdf --averageType median --plotWidth 18 --plotHeight 10 &&
  plotHeatmap --matrixFile $out_path1/${file}.matrix.gz -o $out_path1/${file}.matrix_heatmap.pdf"
echo "========================"

#sbatch --time 03:00:00 --mem 8G --wrap="computeMatrix scale-regions -S ${file} -R $folder_ref_bed/*.genebody_ref.bed -a 1000 -b 1000 --outFileName $out_path1/TSS_TES_${file}.matrix.gz --missingDataAsZero && 
#  plotProfile --matrixFile $out_path1/TSS_TES_${file}.matrix.gz -o $out_path1/TSS_TES_${file}.matrix.pdf --averageType median --plotWidth 18 --plotHeight 10 &&
#  plotHeatmap --matrixFile $out_path1/TSS_TES_${file}.matrix.gz -o $out_path1/TSS_TES_${file}.matrix_heatmap.pdf"
#echo "========================"

done

#pol2 chip
cd $bw_dir_pol2
for file in  *ox*ChIP_*bw
do
sbatch --time 03:00:00 --mem 8G --wrap="computeMatrix reference-point -S ${file} -R $folder_ref_bed/*.genebody_ref.bed -a 1000 -b 1000 --outFileName $out_path2/${file}.matrix.gz --missingDataAsZero && 
  plotProfile --matrixFile $out_path2/${file}.matrix.gz -o $out_path2/${file}.matrix.pdf --averageType median --plotWidth 18 --plotHeight 10 &&
  plotHeatmap --matrixFile $out_path2/${file}.matrix.gz -o $out_path2/${file}.matrix_heatmap.pdf"
echo "========================"

sbatch --time 03:00:00 --mem 8G --wrap="computeMatrix scale-regions -S ${file} -R $folder_ref_bed/*.genebody_ref.bed -a 1000 -b 1000 --outFileName $out_path2/TSS_TES_${file}.matrix.gz --missingDataAsZero && 
  plotProfile --matrixFile $out_path2/TSS_TES_${file}.matrix.gz -o $out_path2/TSS_TES_${file}.matrix.pdf --averageType median --plotWidth 18 --plotHeight 10 &&
  plotHeatmap --matrixFile $out_path2/TSS_TES_${file}.matrix.gz -o $out_path2/TSS_TES_${file}.matrix_heatmap.pdf"
echo "========================"

done

#G4 chip
cd $bw_dir_BG4
for file in  *ox*ChIP*bw
do
sbatch --time 03:00:00 --mem 8G --wrap="computeMatrix reference-point -S ${file} -R $folder_ref_bed/*.genebody_ref.bed -a 1000 -b 1000 --outFileName $out_path3/${file}.matrix.gz --missingDataAsZero && 
  plotProfile --matrixFile $out_path3/${file}.matrix.gz -o $out_path3/${file}.matrix.pdf --averageType median --plotWidth 18 --plotHeight 10 &&
  plotHeatmap --matrixFile $out_path3/${file}.matrix.gz -o $out_path3/${file}.matrix_heatmap.pdf"
echo "========================"

sbatch --time 03:00:00 --mem 8G --wrap="computeMatrix scale-regions -S ${file} -R $folder_ref_bed/*.genebody_ref.bed -a 1000 -b 1000 --outFileName $out_path3/TSS_TES_${file}.matrix.gz --missingDataAsZero && 
  plotProfile --matrixFile $out_path3/TSS_TES_${file}.matrix.gz -o $out_path3/TSS_TES_${file}.matrix.pdf --averageType median --plotWidth 18 --plotHeight 10 &&
  plotHeatmap --matrixFile $out_path3/TSS_TES_${file}.matrix.gz -o $out_path3/TSS_TES_${file}.matrix_heatmap.pdf"
echo "========================"

done
```

In R combine libraries

``` r
ending_label=".bw.matrix.gz"

case_label_hypo='^hypoxia_1hr'
case_label_normo='^normoxia_1hr'

matrix_files_hypo <- list.files(".",paste0(case_label_hypo,".*",ending_label,"$"))

matrix_files_normo <- list.files(".",paste0(case_label_normo,".*",ending_label,"$"))


source('/Users/simeon01/Documents/Karen/Plot_combined_intensities_given_deeptoolsMatrix.R')

atac_hypo <- average_densities_single_case(matrix_files_hypo,'atac_hypo_signal')
atac_normo <- average_densities_single_case(matrix_files_normo,'atac_normo_signal')
atac_hypo_minus_normo <- average_densities_2(matrix_files_hypo,matrix_files_normo,'atac_hypo_signal','atac_normo_signal')

# plot average +-sem
cd /Users/simeon01/Documents/Karen/20201019_BG4hypoxia_normoxia_SLX-198891_SLX-198892/reference_bed_for_figure2_promoters_plus_minus500/atac_seq_signal_at_G4_regions
Rscript /Users/simeon01/Documents/Karen/TSS_generate_plot_with_confidence.R atac_hypo_signal_minus_atac_normo_signal_average.matrix.gz atac_hypo_signal_minus_atac_normo_signal_average_sem
```

After summarising various rep into a single case, use deepTools to generate heatmap. This is done in turn for all individual experimental cases.

``` bash

for file in  **avg*gz
do
  sbatch --time 00:05:00 --mem 1G --wrap="plotHeatmap --matrixFile $file --averageTypeSummaryPlot median -o  ${file%%.matrix.gz}.matrix_heatmap.pdf --zMax 1 --yMax 0.45 --yMin 0"
  echo "========================"
done

for file in  **minus*av*gz
do
  sbatch --time 00:05:00 --mem 1G --wrap="plotHeatmap --matrixFile $file --averageTypeSummaryPlot median -o  ${file%%.matrix.gz}.matrix_heatmap.pdf --zMax 0 --zMin -0.3"
  echo "========================"
done
```

    ## bash: line 3: sbatch: command not found
    ## ========================
    ## bash: line 9: sbatch: command not found
    ## ========================
