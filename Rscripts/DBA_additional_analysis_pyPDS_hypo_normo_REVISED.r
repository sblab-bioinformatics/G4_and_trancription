#' @title Script for perform differential bidngin analysis
#' @description Script that perform differential analysis of samples esposed to pyPDS, DMSO, normoxia, hypoxia. Multiple comparisons are performed (cross-comparisons)
#' @param coverage files where file names encode biological replicate and condition information (path of files); files prefixex and suffixes if multiple cases are in the same path
#' @return Plots (quality controls, non-supervised analysis plots), EdgeR table and plots with differential analysis outcome
#' @author angela simeone  \email{angela.simeone@gmail.com}
# 

## == fuctions  and packages ==
library(dendextend)
library(dplyr)
library(ggdendro)
library(ggplot2)
tss_file <- '/Users/simeon01/Documents/genomes/hg38/gencode.v28.TSS.bed'
source('/Users/simeon01/Documents/Karen/20190613_Karen_continue_pol2_G4_maps/G4_signal_DRB/function_to_export_tables.R')

norm_sign_to_tot_map_reads <- function(matrix,list_stat5_mapped_reads,pattern){
  Library_sizes <- c()
  matrix_norm <- matrix
  for (i in 1:dim(matrix)[2]){
    print(i)
    tmp_cond <- colnames(matrix)[i]
    tmp_cond <- gsub(pattern,'',tmp_cond)
    print(tmp_cond)
    tmp_stat <- grep(tmp_cond,list_stat5_mapped_reads,value = T)
    tot_reads <- as.numeric(system(paste0('cat ',tmp_stat),intern = TRUE))
    print(tot_reads)
    Library_sizes[i] <- tot_reads
    matrix_norm[,i] <- matrix[,i]/tot_reads*1e6
  }
  names(Library_sizes) <- colnames(matrix)
  NewList <- list("norm_cpm_M" = matrix_norm, "lib_size" = Library_sizes)
  return(NewList)
}


################################################
# load data
#################################################

#path1='/Users/simeon01/Documents/Karen/20201014_karen_FE_bg4chip_pyPDS_hypo_normo/coverages_at_various_regions_20201014'
#path1='/Users/simeon01/Documents/Karen/coverages_at_various_regions_20201020'
setwd(path1)
# mv SLX-19270.HNTFYBGXF.s_1.r_1.CAGAGAGG+AGAGGATA.hg38.sort.markduplicates.dba_consensus.20200821.multi2.coverage.bed 1uM_pyPDS_2hr_biorep2_ChIP3_SLX-19270.i708_i505.r_1.hg38.sort.markduplicates.dba_consensus.20200821.multi2.coverage.bed

# ## rename the replicated experiment (used for estimating the batch effect)
# for file in SLX-19270.HNTFYBGXF.s_1.r_1.CAGAGAGG+AGAGGATA.hg38.sort*
#   do
#   echo $file
#   echo ${file/SLX-19270.HNTFYBGXF.s_1.r_1.CAGAGAGG+AGAGGATA.hg38.sort/1uM_pyPDS_2hr_biorep2_ChIP3_SLX-19270.i708_i505.r_1.hg38.sort}
#   echo "===="
#   mv $file ${file/SLX-19270.HNTFYBGXF.s_1.r_1.CAGAGAGG+AGAGGATA.hg38.sort/1uM_pyPDS_2hr_biorep2_ChIP3_SLX-19270.i708_i505.r_1.hg38.sort}
# done
# ## rename the replicated experiment (the first one) - in order to exclude it
# for file in 1uM_pyPDS_2hr_hypoxia_biorep2_ChIP2_SLX-19269.i710_i505.r_1.hg38.sort.markduplicates*bed
# do
#   mv $file ${file}z
#   done
# ## remove lostreads
# # rm *lostreads*bed



#pattern_file_extension="markduplicates.dba_consensus_pyPDS_DMSO_hypo_normo.G4.20201014.coverage.bed$"
#pattern_folder_prefix="dba_consensus_pyPDS_DMSO_hypo_normo.G4.20201014.coverage"

list_files <- list.files(path = ".",pattern = pattern_file_extension)
list_files <- grep(selection_files_pattern,list_files,value = T)
list_files
temp_size1 <- read.table(list_files[1],stringsAsFactors=F)

Matrix_cov <- matrix(ncol = length(list_files), nrow = dim(temp_size1)[1])
for (i in 1:length(list_files)) {
  temp <- read.table(list_files[i],stringsAsFactors = F)
  Matrix_cov[,i] <- temp$V4
  peaks_coordinates <- temp[,c(1,2,3)]
}
colnames(Matrix_cov) <-list_files
colnames(Matrix_cov) <- gsub(pattern_file_extension,'',list_files)
#
# grep groups and define them
#colnames(Matrix_cov) <- gsub('.SLX*.*','',colnames(Matrix_cov))
#colnames(Matrix_cov) <- gsub('.SLX*.*','',colnames(Matrix_cov))

rownames(Matrix_cov) <- paste(peaks_coordinates$V1,peaks_coordinates$V2,peaks_coordinates$V3,sep="_")

# specific samples need to be removed
colnames(Matrix_cov)
grep('Input',colnames(Matrix_cov))
Matrix_cov_selected <- Matrix_cov[,-grep('Input',colnames(Matrix_cov))]


if(input_flag==1){
Matrix_cov_selected <- Matrix_cov[,-grep('Input',colnames(Matrix_cov))]
} else {Matrix_cov_selected=Matrix_cov}


stat5_files <- list.files(path = path1, pattern = "stat5")

norm_Matrix_cov_selected <- norm_sign_to_tot_map_reads(Matrix_cov_selected,stat5_files,".markduplicates.*.")
norm_Matrix_cov <- norm_sign_to_tot_map_reads(Matrix_cov,stat5_files,".markduplicates.*.")

setwd(path1)

# ## == clustering raw data == 
# 
# # # clustering - raw data
# pdf('hclust_raw_counts.pdf')
# h_Matrix_cov_selected <- hclust(dist(t(Matrix_cov_selected)),method="ward.D2")
# p <- ggdendrogram(h_Matrix_cov_selected,rotate=T) +  ggtitle('clustering - raw counts')
# print(p)
# dev.off()
# 
# ## == clustering - data norm to total number of reads == 
# pdf('hclust_RPM.pdf')
# h_norm_Matrix_cov <- hclust(dist(t(norm_Matrix_cov$norm_cpm_M)),method="ward.D2")
# p <- ggdendrogram(h_norm_Matrix_cov,rotate=T) +  ggtitle('clustering - RPM')
# print(p)  
# dev.off()
# 
# # ## == clustering - data norm to total number of reads -- NO INPUT == 
# pdf('hclust_RPM_no_input.pdf')
h_norm_Matrix_cov_selected_noInput <- hclust(dist(t(norm_Matrix_cov_selected$norm_cpm_M)),method="ward.D2")
p <- ggdendrogram(h_norm_Matrix_cov_selected_noInput,rotate=T) +  ggtitle('clustering - RPM - no input')
print(p)
# dev.off()
# print(p)
# 
# ## == clustering - data norm to total number of reads -- NO INPUT == 
# pdf('hclust_RPM_no_input_excluding_low10percent.pdf')
# A <- norm_Matrix_cov_selected$norm_cpm_M
# A_sum <- rowSums(A)
# hist(A_sum)
# summary(A_sum)
# quantile(A_sum, 0.1)
# A_selected <- A[A_sum>quantile(A_sum, 0.1),]
# h_A_selected <- hclust(dist(t(A_selected)),method="ward.D2")
# p <- ggdendrogram(h_A_selected,rotate=T) +  ggtitle('clustering - RPM - no input')
# plot(p)
# dev.off()
# plot(p)
# 
# #use spearman corr
# cols.cor <- cor(A_selected,method = 'spearman')
# h_A_selected_spearman <- hclust(as.dist(1-cols.cor)) 
# p <- ggdendrogram(h_A_selected_spearman,rotate=T) +  ggtitle('clustering - RPM - no input')
# plot(p)
# 
# 
# pdf('hclust_RPM_no_input_excluding_NOlow10percent_NOTop10percent.pdf')
# A <- norm_Matrix_cov_selected$norm_cpm_M
# A_sum <- rowSums(A)
# hist(A_sum)
# summary(A_sum)
# quantile(A_sum, 0.1)
# A_selected <- A[A_sum>quantile(A_sum, 0.1) & A_sum<quantile(A_sum, 0.9) ,]
# #A_selected <- A[A_sum>quantile(A_sum, 0.5) ,]
# h_A_selected <- hclust(dist(t(A_selected)),method="ward.D2")
# p <- ggdendrogram(h_A_selected,rotate=T) +  ggtitle('clustering - RPM - no input')
# plot(p)
# dev.off()
# plot(p)
# 
# pdf('hclust_RPM_no_input_DMSOonly.pdf')
# B <- norm_Matrix_cov_selected$norm_cpm_M[,grep('DMSO',colnames(norm_Matrix_cov_selected$norm_cpm_M))]
# C <- B[,-grep('SLX-19412_SLX-19413',colnames(B))]
# h_B <- hclust(dist(t(B)),method="ward.D2")
# p <- ggdendrogram(h_B,rotate=T) +  ggtitle('clustering - RPM - no input')
# plot(p)
# 
# h_C <- hclust(dist(t(C)),method="ward.D2")
# p_C <- ggdendrogram(h_C,rotate=T) +  ggtitle('clustering - RPM - no input')
# plot(p_C)
# dev.off()
# pdf('hclust_RPM_no_input_pyPDSonly.pdf')
# C <- norm_Matrix_cov_selected$norm_cpm_M[,grep('pyPDS',colnames(norm_Matrix_cov_selected$norm_cpm_M))]
# h_C <- hclust(dist(t(C)),method="ward.D2")
# p <- ggdendrogram(h_C,rotate=T) +  ggtitle('pyPDS lib. clustering - RPM - no input')
# plot(p)
# dev.off()
# 
# 


#################################
# Differential analysis
#################################

# in the differential analysis we will account for:
# bio rep [4 different ones]

# condition [2 experimental conditions DMSO  and Garcinol]

Counts_to_use_all <- Matrix_cov_selected
curr_lib_size <- norm_Matrix_cov_selected$lib_size
if(flag_BG4_at_TSS==1){
# here filter those regions that have an overlap with TSS
#  file_subset_TSS : file containing the set of peaks overlapping TSS
BG4_overlapping_TSS <- read.table(file_subset_TSS,stringsAsFactors = F)

BG4_overlapping_TSS_id <- paste(BG4_overlapping_TSS$V1,BG4_overlapping_TSS$V2,BG4_overlapping_TSS$V3,sep = "_")
Counts_to_use_all <- Counts_to_use_all[rownames(Counts_to_use_all)%in%BG4_overlapping_TSS_id ,]
}

Counts_input <- Matrix_cov[,grep('nput',colnames(Matrix_cov))]
curr_lib_size_input <- norm_Matrix_cov$lib_size[grep('nput',names(norm_Matrix_cov$lib_size))]
temp_names <- colnames(Counts_to_use_all)
conditions <- gsub('_biore.*.','',temp_names)

conditions[grep('DMSO_control_2hr_normoxia',conditions)] <- "DMSO_control_2hr_normoxia"
conditions[grep('1uM_pyPDS_2hr_normoxia',conditions)] <- "1uM_pyPDS_2hr_normoxia"
conditions[grep('DMSO_control_2hr_hypoxia',conditions)] <- "DMSO_control_2hr_hypoxia"
conditions[grep('1uM_pyPDS_2hr_hypoxia',conditions)] <- "1uM_pyPDS_2hr_hypoxia"
conditions <- factor(conditions) # there are 4 factors
levels(conditions)
level_to_use <- grep('DMSO_control_2hr_normoxia',levels(conditions),value = T)
conditions <- relevel(conditions,ref=level_to_use) # relevel as reference base the condition no_TPL_1hr_no_dm6
conditions

#biorep <- sub('_ChIP*.*','',gsub('.*._bio','',temp_names))
biorep <- sub('.SLX*.*','',gsub('.*._bio','',temp_names))
biorep

if (flag_batch==1){
batch <- gsub('.i70*.*|.merged.','',gsub('.*._SLX','SLX',temp_names)) 
batch <- gsub('.merged','',batch) 
batch
}
### complex heatmap here ####################################################
# 
# # complexHeatmap on all BG4 signal
# library(ComplexHeatmap)
# library(circlize)
# ha = HeatmapAnnotation(exp_condition = conditions , annotation_name_side = "left")
# ht_list = Heatmap(norm_Matrix_cov_selected$norm_cpm_M, name = "BG4signal", column_km  = 4,row_km = 4,
#                   col = colorRamp2(c(0, 0.6, 10), c("white", "pink", "red")),
#                   #top_annotation = ha, 
#                   show_column_names = T,show_row_names = FALSE, row_title = NULL, show_row_dend = FALSE,show_column_dend = T) 
# 
# pdf('clust_BG4signal.pdf', width = 8,height = 7)
# draw(ht_list)
# dev.off()


# conditions_D <- gsub('_bio.*.','',colnames(norm_Matrix_cov_selected$norm_cpm_M)[-grep('4h|SLX-19267',colnames(norm_Matrix_cov_selected$norm_cpm_M))])
# 
# ha_D = HeatmapAnnotation(exp_condition = conditions_D , annotation_name_side = "left")
# 
# ht_list = Heatmap(D, name = "BG4signal", column_km  = 4,row_km = 5,
#                   col = colorRamp2(c(0, 0.6, 10), c("white", "pink", "red")),
#                   top_annotation = ha_D, 
#                   show_column_names = T,show_row_names = FALSE, row_title = NULL, show_row_dend = FALSE,show_column_dend = T) 


#############################################################################
# filter out pyPDS_normoxia case
if (filter_pyPDSnormo==1){
  
  ind_pyPDSnormo_to_remove <- grep("^1uM_pyPDS_2hr$",conditions)
  conditions <- conditions[-ind_pyPDSnormo_to_remove]
  biorep <- biorep[-ind_pyPDSnormo_to_remove] 
  batch <- batch[-ind_pyPDSnormo_to_remove]
  norm_Matrix_cov_selected$lib_size <-norm_Matrix_cov_selected$lib_size[-ind_pyPDSnormo_to_remove]
  norm_Matrix_cov_selected$norm_cpm_M <-norm_Matrix_cov_selected$norm_cpm_M[,-ind_pyPDSnormo_to_remove]
  
  Counts_to_use_all <- Counts_to_use_all[,-ind_pyPDSnormo_to_remove]
}


#############################################################################


curr_lib_size <- norm_Matrix_cov_selected$lib_size

if (flag_batch==1){
design <- model.matrix(~0 + conditions+biorep + batch)
} else {
design <- model.matrix(~0 + conditions)
}

print(colnames(design))


# create EdgeR object
library(edgeR)
y <- DGEList(counts = Counts_to_use_all, group= conditions)
stopifnot(rownames(y$samples) == names(curr_lib_size))
y$samples$lib.size<- curr_lib_size
cpm_y <- cpm(Counts_to_use_all,lib.size = curr_lib_size)

dir_to_use <- paste0(path1,'/',pattern_folder_prefix,'_','exclude_nothing_include_all_pval0.05')
dir.create(dir_to_use)
setwd(dir_to_use)

if (input_flag==1){
print("flag_input is one")
# do the same for input libraries
y_input <- DGEList(counts = Counts_input)
stopifnot(rownames(y_input$samples) == names(curr_lib_size_input))
y_input$samples$lib.size<- curr_lib_size_input
cpm_y_input <- cpm(Counts_input,lib.size = curr_lib_size_input)
#threshold_over_input =2 # ive used this for most cases
treshold_to_use <- threshold_over_input * quantile(rowMeans(cpm_y_input),0.99)
message("======================================")
message(paste0('|| threshold_to_use:',treshold_to_use))
message(paste0('|| 99quantile:',quantile(rowMeans(cpm_y_input),0.99)))
message("======================================")
 

Counts_to_use_all2 <- Counts_to_use_all
curr_lib_size2 <- norm_Matrix_cov_selected$lib_size
temp_names2 <- colnames(Counts_to_use_all2)
conditions2 <- gsub('_biore.*.','',temp_names2)

conditions2[grep('DMSO_control_2hr_normoxia',conditions2)] <- "DMSO_control_2hr_normoxia"
conditions2[grep('1uM_pyPDS_2hr_normoxia',conditions2)] <- "1uM_pyPDS_2hr_normoxia"
conditions2[grep('DMSO_control_2hr_hypoxia',conditions2)] <- "DMSO_control_2hr_hypoxia"
conditions2[grep('1uM_pyPDS_2hr_hypoxia',conditions2)] <- "1uM_pyPDS_2hr_hypoxia"
conditions2 <- factor(conditions2) # there are 4 factors
levels(conditions2)
level_to_use <- grep('DMSO_control_2hr_normoxia',levels(conditions2),value = T)
conditions2 <- relevel(conditions2,ref=level_to_use) # relevel as reference base the condition no_TPL_1hr_no_dm6
conditions2
biorep2 <- biorep

design2 <- model.matrix(~ 0+ conditions2+ biorep2 )
#design2 <- model.matrix(~ 0 + conditions2)
colnames(design2)
y2 <- DGEList(counts = Counts_to_use_all2, group= conditions2)
stopifnot(rownames(y2$samples) == names(curr_lib_size2))
y2$samples$lib.size<- curr_lib_size2
cpm_y2 <- cpm(Counts_to_use_all2,lib.size = curr_lib_size2)

intermediate_cpm_y2 <- data.frame(cpm_y2)
unique_set_conditions <- unique(conditions2)
updated_intermediate_cpm_y2 <- matrix(nrow = dim(intermediate_cpm_y2)[1],ncol=length(unique_set_conditions))
for (i in 1:length(unique_set_conditions)){
  print(grep(unique_set_conditions[i],colnames(intermediate_cpm_y2)))
  temp_avg_count_condition <- rowMeans(intermediate_cpm_y2[,grep(unique_set_conditions[i],colnames(intermediate_cpm_y2))])
  updated_intermediate_cpm_y2[,i] <- temp_avg_count_condition
  rm(temp_avg_count_condition)
}

ind_selected <- rowSums(updated_intermediate_cpm_y2 >= treshold_to_use) >= 1
print(head(updated_intermediate_cpm_y2))
message("=================================================================")
message(c("=====> Matrix dimentsions before filtering: ",dim(y2)[1],' rows ',dim(y2)[2],' cols <====='))
tmp_dim1 <- dim(y2)
y2 = y2[ind_selected,] # at least 1 cell type has a average cpm > 1
tmp_dim1 <- rbind(tmp_dim1,dim(y2))
rownames(tmp_dim1) <- c('initial_set','filtered_set')
write.table(tmp_dim1,file=paste0(dir_to_use,'/summary_all_cases_DBA.hypo_vs_normo',case,'_size_dataset.csv'),quote = F, sep = ",", col.names = NA)
message(c("=====> Matrix dimentsions after filtering: ",dim(y2)[1],' rows ',dim(y2)[2],' cols <====='))
message("=================================================================")
}

if(input_flag==0) {
  Counts_to_use_all2 <- Counts_to_use_all
  curr_lib_size2 <- norm_Matrix_cov_selected$lib_size
  conditions2 <- gsub('_biore.*.','',temp_names)
  conditions2 <- factor(conditions2) # there are 4 factors
  level_to_use <- grep('DMSO_control_2hr_normoxia',levels(conditions2),value = T)
  if (input_flag==1){
    level_to_use <- grep('ChIP',levels(conditions2),value = T)
  }
  conditions2 <- relevel(conditions2,ref=level_to_use) # relevel as reference base the condition no_TPL_1hr_no_dm6
  conditions2
  biorep2 <- biorep
  
  design2 <- model.matrix(~ 0+ conditions2+ biorep2 )
  #design2 <- model.matrix(~ 0 + conditions2)
  y2 <- DGEList(counts = Counts_to_use_all2, group= conditions2)
  stopifnot(rownames(y2$samples) == names(curr_lib_size2))
  y2$samples$lib.size<- curr_lib_size2
  cpm_y2 <- cpm(Counts_to_use_all2,lib.size = curr_lib_size2) 
}

y2 <- calcNormFactors(y2, method= 'none')
y2 <- estimateDisp(y2,design2)
y2 <- estimateGLMCommonDisp(y2,design2)
y2 <- estimateGLMTagwiseDisp(y2,design2)


# check MDS plot
pdf('MDS_edger_mdsplot_g4.pdf', 7, 7)
par(las= 1, cex=1.2)
col_map <- conditions2
levels(col_map) <- 1:length(levels(col_map))
col_map <- as.numeric(factor(conditions2))
temp_rainbow <- rainbow(7)
col_map_rainbow <- temp_rainbow[col_map]
plotMDS(y2,
        #pch=c(0,3,0,3,0,3),
        #xlim=c(-5,8), ylim= c(-1.5,0.5),
        #xlim=c(-1.5,1.5), ylim= c(-1,0.7),
        #xaxp= c(-40, 40, 8),  yaxp= c(-10,10, 8),
        labels= paste0(conditions2,'-',biorep2),
        col=col_map_rainbow,
        main= 'MDS plot ', cex= 0.6)
dev.off()

pdf('MDS_edger_mdsplot_g4_symb.pdf', 7, 7)
par(las= 1, cex=1.2)
col_map <- conditions2
levels(col_map) <- 1:length(levels(col_map))
col_map <- as.numeric(factor(conditions2))
temp_rainbow <- rainbow(6)
col_map_rainbow <- temp_rainbow[col_map]
plotMDS(y2,
        pch=col_map,
        #xlim=c(-1.1,1.3), ylim= c(-1.3,0.4),
        #xlim=c(-1.5,1.5), ylim= c(-1,0.7),
        #xaxp= c(-40, 40, 8),  yaxp= c(-10,10, 8),
        #labels= paste(y$samples$group,bio_rep, sep= '--'),
        col=col_map_rainbow,
        main= 'MDS plot ', cex= 0.7)
grid()
dev.off()

pdf('PCA2_edger_mdsplot_g4.pdf', 9, 7)
pcaResult<-prcomp(t(y2$counts))
plot(pcaResult$x,
     main= 'Principal components, counts',
     xlab="PC1",
     ylab="PC2",
     #xlab= sprintf('PC1 (sd: %s%%)', round(100 * (pcaResult$sdev[1] / sum(pcaResult$sdev)))),
     #ylab= sprintf('PC2 (sd: %s%%)', round(100 * (pcaResult$sdev[2] / sum(pcaResult$sdev)))),
     type= 'n')#, xlim = c(-55000,5000),ylim =c(-1000,1000))
text(x= pcaResult$x[,1], y= pcaResult$x[,2], labels= paste0(conditions2,'-',biorep2), cex= 0.5,
     #col= c(rep('#CCFF00FF',5),rep('#00FF66FF',5),rep('#FF0000FF',5)))
     col = col_map_rainbow)
dev.off()


# fit the generalized linear model (quasi-likelihood negative binomial generalized log-linear model)
fit <- glmQLFit(y2, design2)
colnames(fit)
# [1] "conditions2DMSO_control_2hr"         "conditions21uM_pyPDS_2hr"            "conditions21uM_pyPDS_2hr_hypoxia"   
# [4] "conditions2DMSO_control_2hr_hypoxia" "biorep2rep2"                         "biorep2rep3"                        
# [7] "biorep2rep4"  

function_run_glmQLFTest <- function(fit_to_use,contrast_to_use,label_case,flag_exclude_print_bed=T){
  qlf <- glmQLFTest(fit_to_use, contrast = contrast_to_use) 
  #qlf <- glmLRT(fit_to_use, contrast = contrast_to_use) 
  de_pairing <- topTags(qlf, n= nrow(y$counts))$table
  res <- decideTestsDGE(qlf,adjust.method="none", p.value = 0.05,lfc=0)
  summary(res)
  summary_res <- summary(res)
  res
  
  if(flag_exclude_print_bed==T){
    tss_res <- export_bed_and_TSS_in_bed_regions(tss_file,'./',label_case,res)
    list_res <- list("summary_res"=summary_res,"res"=res,"tss_res"=tss_res,'edgeRtable'=de_pairing)
  }
  
  pdf(paste0('MD_MA_plot_differential_binding_analysis_',label_case,'.pdf'),8,7)
  # plotMD(qlf,status = res,values =c(1,0,-1),col=c("red","black","blue"), main = label_case)
  # title(paste0("\n\nDEup:",length(which(res==1)),
  #              " - DEdown: ",length(which(res==-1)),
  #              " - notDE: ",length(which(res==0))))
  plotMD(qlf, status = res,
         values =c(1,0,-1),
         col=c("#ff000080","#00000030","#0000ff80"),
         pch = 19,cex=c(1,0.8,1))
  title(main = list( paste0("\n\n\npval0.05 DEup:",
                            length(which(res==1)),
                            " - DEdown: ",length(which(res==-1)),
                            " - Tot: ",length(res==0)),
                     cex = .7,font = 3))
  dev.off()
  # write EdgeR table
  write.table(de_pairing,file= paste0('EdgeR_dba_',label_case,'.csv'),quote = F, sep = ",", col.names = NA)
  
  list_res <- list("summary_res"=summary_res,"res"=res,'edgeRtable'=de_pairing)
  return(list_res)
}


##  =========== DMSO2hrhypoxia_vs_DMSO2h ================
contrast_DMSO2hrhypoxia_vs_DMSO2h <- rep(0, length(colnames(fit)))
A <- grep("DMSO_control_2hr_hypoxia",colnames(fit))
B <- grep("DMSO_control_2hr_normoxia",colnames(fit))
A
B
# A vs B
contrast_DMSO2hrhypoxia_vs_DMSO2h[A] <- 1
contrast_DMSO2hrhypoxia_vs_DMSO2h[B] <- -1
colnames(fit)[which(contrast_DMSO2hrhypoxia_vs_DMSO2h!=0)]
label_case <-  "DMSO2hrhypoxia_vs_DMSO2h"



A_res_DMSO2hrhypoxia_vs_DMSO2h <- function_run_glmQLFTest(fit,contrast_DMSO2hrhypoxia_vs_DMSO2h, label_case,TRUE)
summary(A_res_DMSO2hrhypoxia_vs_DMSO2h$res)




##  =========== conditions21uM_pyPDS_2hr_hypoxia_vs_conditions2DMSO_control_2hr ================
contrast_1uM_pyPDS_2hr_hypoxia_vs_DMSO2h <- rep(0, length(colnames(fit)))
A <- grep("1uM_pyPDS_2hr_hypoxia",colnames(fit))
B <- grep("DMSO_control_2hr_normoxia",colnames(fit))
A
B
# A vs B
contrast_1uM_pyPDS_2hr_hypoxia_vs_DMSO2h[A] <- 1
contrast_1uM_pyPDS_2hr_hypoxia_vs_DMSO2h[B] <- -1
# [1] "conditions2DMSO_control_2hr"         "conditions21uM_pyPDS_2hr"            "conditions21uM_pyPDS_2hr_hypoxia"   
# [4] "conditions2DMSO_control_2hr_hypoxia" "biorep2rep2"                         "biorep2rep3"                        
# [7] "biorep2rep4" 
colnames(fit)[which(contrast_1uM_pyPDS_2hr_hypoxia_vs_DMSO2h!=0)]
label_case <-  "1uM_pyPDS_2hr_hypoxia_vs_DMSO2h"
A_res_1uM_pyPDS_2hr_hypoxia_vs_DMSO2h <- function_run_glmQLFTest(fit,contrast_1uM_pyPDS_2hr_hypoxia_vs_DMSO2h, label_case,TRUE)
summary(A_res_1uM_pyPDS_2hr_hypoxia_vs_DMSO2h$res)

# ## paper_recheck data## ============================
# qlf <- glmQLFTest(fit, contrast = contrast_1uM_pyPDS_2hr_hypoxia_vs_DMSO2h) 
# de_pairing <- topTags(qlf, n= nrow(y$counts))$table
# res <- decideTestsDGE(qlf,adjust.method="fdr", p.value = 0.05,lfc=0)
# summary(res)
# summary_res <- summary(res)
# res
# ## paper_recheck data## ============================

if (filter_pyPDSnormo==0){
##  =========== conditions21uM_pyPDS_2hr_vs_conditions2DMSO_control_2hr ================
contrast_1uM_pyPDS_2hr_vs_DMSO2h <- rep(0, length(colnames(fit)))
A <- grep("1uM_pyPDS_2hr_normoxia",colnames(fit))
B <- grep("DMSO_control_2hr_normoxia",colnames(fit))
A
B
# A vs B
contrast_1uM_pyPDS_2hr_vs_DMSO2h[A] <- 1
contrast_1uM_pyPDS_2hr_vs_DMSO2h[B] <- -1
# [1] "conditions2DMSO_control_2hr"         "conditions21uM_pyPDS_2hr"            "conditions21uM_pyPDS_2hr_hypoxia"   
# [4] "conditions2DMSO_control_2hr_hypoxia" "biorep2rep2"                         "biorep2rep3"                        
# [7] "biorep2rep4" 
colnames(fit)[which(contrast_1uM_pyPDS_2hr_vs_DMSO2h!=0)]
label_case <-  "1uM_pyPDS_2hr_vs_DMSO2h"
A_res_1uM_pyPDS_2hr_vs_DMSO2h <- function_run_glmQLFTest(fit,contrast_1uM_pyPDS_2hr_vs_DMSO2h, label_case, TRUE)
summary(A_res_1uM_pyPDS_2hr_vs_DMSO2h$res)

##  =========== conditions21uM_pyPDS_2hr_hypoxia_vs_conditions21uM_pyPDS_2hr ================
contrast_1uM_pyPDS_2hr_hypoxia_vs_1uM_pyPDS_2hr <- rep(0, length(colnames(fit)))
A <- grep("1uM_pyPDS_2hr_hypoxia",colnames(fit))
B <- grep("1uM_pyPDS_2hr_normoxia",colnames(fit))
# A vs B
contrast_1uM_pyPDS_2hr_hypoxia_vs_1uM_pyPDS_2hr[A] <- 1
contrast_1uM_pyPDS_2hr_hypoxia_vs_1uM_pyPDS_2hr[B] <- -1
# [1] "conditions2DMSO_control_2hr"         "conditions21uM_pyPDS_2hr"            "conditions21uM_pyPDS_2hr_hypoxia"   
# [4] "conditions2DMSO_control_2hr_hypoxia" "biorep2rep2"                         "biorep2rep3"                        
# [7] "biorep2rep4" 

label_case <-  "1uM_pyPDS_2hr_hypoxia_vs_1uM_pyPDS_2hr"
A_res_1uM_pyPDS_2hr_hypoxia_vs_1uM_pyPDS_2hr <- function_run_glmQLFTest(fit,contrast_1uM_pyPDS_2hr_hypoxia_vs_1uM_pyPDS_2hr, label_case, TRUE)
summary(A_res_1uM_pyPDS_2hr_hypoxia_vs_1uM_pyPDS_2hr$res)
}
##  =========== 1uM_pyPDS2hr_hypoxia_vs_DMSO2hr_hypoxia ================
contrast_1uM_pyPDS2hr_hypoxia_vs_DMSO2hr_hypoxia <- rep(0, length(colnames(fit)))
A <- grep("1uM_pyPDS_2hr_hypoxia",colnames(fit))
B <- grep("DMSO_control_2hr_hypoxia",colnames(fit))
A
B
# A vs B
contrast_1uM_pyPDS2hr_hypoxia_vs_DMSO2hr_hypoxia[A] <- 1
contrast_1uM_pyPDS2hr_hypoxia_vs_DMSO2hr_hypoxia[B] <- -1
# [1] "conditions2DMSO_control_2hr"         "conditions21uM_pyPDS_2hr"            "conditions21uM_pyPDS_2hr_hypoxia"   
# [4] "conditions2DMSO_control_2hr_hypoxia" "biorep2rep2"                         "biorep2rep3"                        
# [7] "biorep2rep4" 
label_case <-  "1uM_pyPDS2hr_hypoxia_vs_DMSO2hr_hypoxia"
A_res_1uM_pyPDS2hr_hypoxia_vs_DMSO2hr_hypoxia <- function_run_glmQLFTest(fit,contrast_1uM_pyPDS2hr_hypoxia_vs_DMSO2hr_hypoxia, label_case, TRUE)
summary(A_res_1uM_pyPDS2hr_hypoxia_vs_DMSO2hr_hypoxia$res)

# A_res_1uM_pyPDS2hr_hypoxia_vs_DMSO2hr_hypoxiasummary_res
# A_res_DMSO2hrhypoxia_vs_DMSO2h
# A_res_1uM_pyPDS_2hr_hypoxia_vs_DMSO2h
# A_res_1uM_pyPDS_2hr_vs_DMSO2h
# A_res_1uM_pyPDS_2hr_hypoxia_vs_1uM_pyPDS_2hr
# A_res_1uM_pyPDS2hr_hypoxia_vs_DMSO2hr_hypoxia

if (filter_pyPDSnormo==0){
SUMM <- cbind(A_res_DMSO2hrhypoxia_vs_DMSO2h$summary_res,
              A_res_1uM_pyPDS_2hr_hypoxia_vs_DMSO2h$summary_res,
              A_res_1uM_pyPDS_2hr_vs_DMSO2h$summary_res,
              A_res_1uM_pyPDS_2hr_hypoxia_vs_1uM_pyPDS_2hr$summary_res,
              A_res_1uM_pyPDS2hr_hypoxia_vs_DMSO2hr_hypoxia$summary_res)

}

if (filter_pyPDSnormo==1){
  SUMM <- cbind(A_res_DMSO2hrhypoxia_vs_DMSO2h$summary_res,
                A_res_1uM_pyPDS_2hr_hypoxia_vs_DMSO2h$summary_res,
                A_res_1uM_pyPDS2hr_hypoxia_vs_DMSO2hr_hypoxiasummary_res)

}
print(SUMM)
#write.table(SUMM,file='summary_all_cases_DBA.pyPDS_PDS_hypoxia_DMSO_FDR0.1.csv',quote = F, sep = ",", col.names = NA)
write.table(SUMM,file=paste0('summary_all_cases_DBA.pyPDS_PDS_hypoxia_DMSO_case',case,'.csv'),quote = F, sep = ",", col.names = NA)

write.table(table(conditions2),file=paste0('summary_all_cases_DBA.pyPDS_PDS_hypoxia_DMSO_case',case,'_Summary_table_Nlibraries_per_conditions.csv'),quote = F, sep = ",", col.names = NA)
write.table(cbind(as.character(conditions2),biorep2),file=paste0('summary_all_cases_DBA.pyPDS_PDS_hypoxia_DMSO_case',case,'_experimental_design.csv'),quote = F, sep = ",", col.names = NA)

