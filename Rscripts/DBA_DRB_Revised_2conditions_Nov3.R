#' @title Script for perform differential binding analysis
#' @description Script that perform differential analysis of samples esposed to DRB(in DMSO) vs DMSO (control).
#' @param coverage files where file names encode biological replicate and condition information
#' @return Plots (quality controls, non-supervised analysis plots), EdgeR table and plots with differential analysis outcome
#' @author angela simeone  \email{angela.simeone@gmail.com}
# 

#setwd("/Users/simeon01/Documents/Karen/Differential_analysis_pyPDS_hypoxia_normoxia")

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

path1='~/DRB_coverages'

setwd(path1)
pattern_file_extension=".consensus_noDRB_DRB_1hr.merged.bed"
list_files <- list.files(path = ".",pattern = pattern_file_extension)
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


stat5_files <- list.files(path = path1, pattern = "stat5")

norm_Matrix_cov_selected <- norm_sign_to_tot_map_reads(Matrix_cov_selected,stat5_files,".consensus_noDRB_DRB_1hr.merged.bed")
norm_Matrix_cov <- norm_sign_to_tot_map_reads(Matrix_cov,stat5_files,".consensus_noDRB_DRB_1hr.merged.bed")

setwd(path1)


## == clustering raw data == 
# # clustering - raw data
pdf('hclust_raw_counts.pdf')
h_Matrix_cov_selected <- hclust(dist(t(Matrix_cov_selected)),method="ward.D2")
p <- ggdendrogram(h_Matrix_cov_selected,rotate=T) +  ggtitle('clustering - raw counts')
print(p)
dev.off()

## == clustering - data norm to total number of reads == 
pdf('hclust_RPM.pdf')
h_norm_Matrix_cov <- hclust(dist(t(norm_Matrix_cov$norm_cpm_M)),method="ward.D2")
p <- ggdendrogram(h_norm_Matrix_cov,rotate=T) +  ggtitle('clustering - RPM')
print(p)  
dev.off()

# ## == clustering - data norm to total number of reads -- NO INPUT == 
pdf('hclust_RPM_no_input.pdf')
h_norm_Matrix_cov_selected_noInput <- hclust(dist(t(norm_Matrix_cov_selected$norm_cpm_M)),method="ward.D2")
p <- ggdendrogram(h_norm_Matrix_cov_selected_noInput,rotate=T) +  ggtitle('clustering - RPM - no input')
print(p)
dev.off()
print(p)


## == clustering - data norm to total number of reads -- NO INPUT == 
pdf('hclust_RPM_no_input_excluding_low10percent.pdf')
A <- norm_Matrix_cov_selected$norm_cpm_M
A_sum <- rowSums(A)
#hist(A_sum)
summary(A_sum)
quantile(A_sum, 0.1)
A_selected <- A[A_sum<quantile(A_sum, 0.1),]
h_A_selected <- hclust(dist(t(A_selected)),method="ward.D2")
p <- ggdendrogram(h_A_selected,rotate=T) +  ggtitle('clustering - RPM - no input')
plot(p)
dev.off()

pdf('hclust_RPM_no_input_excluding_NOlow10percent_NOTop10percent.pdf')
A <- norm_Matrix_cov_selected$norm_cpm_M
A_sum <- rowSums(A)
#hist(A_sum)
summary(A_sum)
quantile(A_sum, 0.1)
A_selected <- A[A_sum>quantile(A_sum, 0.1) & A_sum<quantile(A_sum, 0.9) ,]
#A_selected <- A[A_sum>quantile(A_sum, 0.5) ,]
h_A_selected <- hclust(dist(t(A_selected)),method="ward.D2")
p <- ggdendrogram(h_A_selected,rotate=T) +  ggtitle('clustering - RPM - no input')
plot(p)
dev.off()


#################################
# Differential analysis
#################################

# in the differential analysis we will account for:
# bio rep [4 different ones]

# condition [2 experimental conditions DMSO  and Garcinol]

Counts_to_use_all <- Matrix_cov_selected
#filter bed intervals of interest
if(flag_BG4_at_TSS==1){
  # here filter those regions that have an overlap with TSS
  #  file_subset_TSS : file containing the set of peaks overlapping TSS
  BG4_overlapping_TSS <- read.table(file_subset_TSS,stringsAsFactors = F)
  BG4_overlapping_TSS_id <- paste(BG4_overlapping_TSS$V1,BG4_overlapping_TSS$V2,BG4_overlapping_TSS$V3,sep = "_")
  Counts_to_use_all <- Counts_to_use_all[rownames(Counts_to_use_all)%in%BG4_overlapping_TSS_id ,]
  rm(BG4_overlapping_TSS,BG4_overlapping_TSS)
}

# create the input counts to use later for filtering step # 
Counts_input <- Matrix_cov[,grep('nput',colnames(Matrix_cov))]
curr_lib_size_input <- norm_Matrix_cov$lib_size[grep('nput',names(norm_Matrix_cov$lib_size))]

temp_names <- colnames(Counts_to_use_all)
num_DRB_1h <- length(grep('^DRB_1h',temp_names))
num_noDRB <- length(grep('^noDRB',temp_names))

                                          
conditions <- temp_names
conditions[grep('^DRB_1h',temp_names)] <- "DRB_1hr"
conditions[grep('^noDRB',temp_names)] <- "noDRB"

conditions 
# [1] "DRB_1hr" "DRB_1hr" "DRB_1hr" "DRB_1hr" "DRB_1hr" "DRB_1hr" "DRB_1hr"
# [8] "DRB_1hr" "DRB_1hr" "DRB_1hr" "DRB_1hr" "DRB_1hr" "noDRB"   "noDRB"  
#[15] "noDRB"   "noDRB"   "noDRB"   "noDRB"   "noDRB"   "noDRB"   "noDRB"  
#[22] "noDRB"   "noDRB"   "noDRB"  

conditions <- factor(conditions) # there are 2 factors
conditions <- relevel(conditions,ref="noDRB") # relevel as reference base the condition no_TPL_1hr_no_dm6
biorep <- temp_names
biorep[grep('.*biorep1.*.SLX-16406',temp_names)] <- "biorep1_16406"
biorep[grep('.*biorep1.*.16767',temp_names)] <- "biorep1_16767"
biorep[grep('biorep3',temp_names)] <- "biorep3"
biorep[grep('biorep5',temp_names)] <- "biorep5"
biorep
# [1] "biorep1_16406" "biorep1_16767" "biorep1_16406" "biorep1_16767"
# [5] "biorep1_16406" "biorep1_16767" "biorep3"       "biorep3"      
# [9] "biorep3"       "biorep5"       "biorep5"       "biorep5"      
#[13] "biorep1_16767" "biorep1_16767" "biorep1_16767" "biorep3"      
#[17] "biorep3"       "biorep3"       "biorep5"       "biorep5"      
#[21] "biorep5"       "biorep1_16406" "biorep1_16406" "biorep1_16406"

curr_lib_size <- norm_Matrix_cov_selected$lib_size
design <- model.matrix(~ conditions+biorep)
colnames(design)
#[1] "(Intercept)"         "conditionsDRB_1hr"   "biorepbiorep1_16767"
#[4] "biorepbiorep3"       "biorepbiorep5"               "biorepbiorep5"  

#dir_selected='/Users/simeon01/Documents/Karen/20201023_revise_TPL_DRB/DRB_coverages/only_2Cond_DRBL_updated_filter_over_input_thr2_TSS_plus_minus500'
dir.create(dir_selected)
setwd(dir_selected)

# create EdgeR object
library(edgeR)
y <- DGEList(counts = Counts_to_use_all, group= conditions)
stopifnot(rownames(y$samples) == names(curr_lib_size))
y$samples$lib.size<- curr_lib_size
cpm_y <- cpm(y, lib.size = curr_TPL_lib_size)

### filtering step with Input signal
y_input <- DGEList(counts = Counts_input)
stopifnot(rownames(y_input$samples) == names(curr_lib_size_input))
y_input$samples$lib.size<- curr_lib_size_input
cpm_y_input <- cpm(Counts_input,lib.size = curr_lib_size_input)
#threshold_over_input = 0#2 # ive used this for most cases
treshold_to_use <- threshold_over_input * quantile(rowMeans(cpm_y_input),0.99)
message("======================================")
message("======================================")
message(paste0('|| threshold_to_use:',treshold_to_use))
message(paste0('|| 99quantile:',quantile(rowMeans(cpm_y_input),0.99)))
message("======================================")


## check if genomic sites are above threshold
intermediate_cpm_y2 <- data.frame(cpm_y)
unique_set_conditions <- unique(conditions)
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
message(c("=====> Matrix dimentsions before filtering: ",dim(y)[1],' rows ',dim(y)[2],' cols <====='))
tmp_dim1 <- dim(y)
y = y[ind_selected,] # at least 1 cell type has a average cpm > 1
tmp_dim1 <- rbind(tmp_dim1,dim(y))
rownames(tmp_dim1) <- c('initial_set','filtered_set')
write.table(tmp_dim1,file='summary_all_cases_DBA.TPL_vs_control_size_dataset.csv',quote = F, sep = ",", col.names = NA)
message(c("=====> Matrix dimentsions after filtering: ",dim(y)[1],' rows ',dim(y)[2],' cols <====='))
message("=================================================================")


y <- calcNormFactors(y, method= 'none')
y <- estimateDisp(y,design)
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)

# check MDS plot
pdf('MDS_edger_mdsplot_g4.pdf', 7, 7)
par(las= 1, cex=1.2)
col_map <- conditions
levels(col_map) <- 1:length(levels(col_map))
col_map <- as.numeric(factor(conditions))
temp_rainbow <- rainbow(7)
col_map_rainbow <- temp_rainbow[col_map]
plotMDS(y,
        #pch=c(0,3,0,3,0,3),
        #xlim=c(-3,3), ylim= c(-0.5,0.5), 
        #xlim=c(-1.5,1.5), ylim= c(-1,0.7), 
        #xaxp= c(-40, 40, 8),  yaxp= c(-10,10, 8),
        labels= y$samples$group,
        col=col_map_rainbow,
        main= 'MDS plot ', cex= 0.7)
grid()
dev.off()

pdf('MDS_edger_mdsplot_g4_symb.pdf', 7, 7)
par(las= 1, cex=1.2)
col_map <- conditions
levels(col_map) <- 1:length(levels(col_map))
col_map <- as.numeric(factor(conditions))
temp_rainbow <- rainbow(7)
col_map_rainbow <- temp_rainbow[col_map]
plotMDS(y,
        pch=c(rep(0,5),rep(1,5)),
        #xlim=c(-3,3), ylim= c(-0.5,0.5), 
        #xlim=c(-1.5,1.5), ylim= c(-1,0.7), 
        #xaxp= c(-40, 40, 8),  yaxp= c(-10,10, 8),
        #labels= paste(y$samples$group,bio_rep, sep= '--'),
        col=col_map_rainbow,
        main= 'MDS plot ', cex= 0.7)
grid()
dev.off()

pdf('PCA2_edger_mdsplot_g4.pdf', 9, 7)
pcaResult<-prcomp(t(y$counts))
plot(pcaResult$x,
     main= 'Principal components, counts',
     xlab="PC1",
     ylab="PC2",
     #xlab= sprintf('PC1 (sd: %s%%)', round(100 * (pcaResult$sdev[1] / sum(pcaResult$sdev)))),
     #ylab= sprintf('PC2 (sd: %s%%)', round(100 * (pcaResult$sdev[2] / sum(pcaResult$sdev)))),
     type= 'n')
text(x= pcaResult$x[,1], y= pcaResult$x[,2], labels= gsub('.con.*.','',temp_names), cex= 0.9, 
     #col= c(rep('#CCFF00FF',5),rep('#00FF66FF',5),rep('#FF0000FF',5)))
     col = col_map_rainbow)
dev.off()

dir.create(dir_selected)
setwd(dir_selected)

# fit the generalized linear model (quasi-likelihood negative binomial generalized log-linear model)
fit <- glmQLFit(y, design)
colnames(fit)
#[1] "(Intercept)"              "conditionsB_DRB_1hr"      "conditionsC_postDRB_1hr"  "conditionsD_postDRB_24hr" "biorepbiorep1_16767"     
#[6] "biorepbiorep3"            "biorepbiorep5"  

##  =========== DRB_1hr_vs_noDRB ================
contrast_DRB_vs_noDRB <- c(0,1,0,0,0)
colnames(fit)[which(contrast_DRB_vs_noDRB==1)]
label_case <-  "_DRB_vs_noDRB_fdr0.01"

qlf.DRB_vs_noDRB <- glmQLFTest(fit, contrast = contrast_DRB_vs_noDRB) 

de_pairing_DRB_vs_noDRB <- topTags(qlf.DRB_vs_noDRB, n= nrow(y$counts))$table

write.table(de_pairing_DRB_vs_noDRB,file = paste0('EdgeRtable_differential_binding_analysis_',label_case,'.csv'), col.names = NA,sep = ",",quote = F)


res_DRB_vs_noDRB <- decideTestsDGE(qlf.DRB_vs_noDRB,adjust.method="fdr", p.value = 0.01,lfc=0)
#summary(res_DRB_vs_noDRB)

summary_res_DRB_vs_noDRB <- summary(res_DRB_vs_noDRB)

tss_res_PDS <- export_bed_and_TSS_in_bed_regions(tss_file,'./',label_case,res_DRB_vs_noDRB)

pdf(paste0('MD_MA_plot_differential_binding_analysis_',label_case,'.pdf'),8,7)
plotMD(qlf.DRB_vs_noDRB,status = res_DRB_vs_noDRB,values =c(1,0,-1),col=c("#ff000080","#00000030","#0000ff80"),pch = 19,cex=c(1,0.8,1))
title(paste0("\n\nDEup:",length(which(res_DRB_vs_noDRB==1)),
                           " - DEdown: ",length(which(res_DRB_vs_noDRB==-1)),
                           " - notDE: ",length(which(res_DRB_vs_noDRB==0))," - Tot: ",length(res_DRB_vs_noDRB)))
dev.off()

message("-------------------------- .. -----------------------------")
message('fdr=0.01')
print(cbind(summary(res_DRB_vs_noDRB), summary(res_DRB_vs_noDRB)/sum(summary(res_DRB_vs_noDRB))))
res_DRB_vs_noDRB_fdr0.01 <- res_DRB_vs_noDRB
message("-------------------------- .. -----------------------------")

res_DRB_vs_noDRB <- decideTestsDGE(qlf.DRB_vs_noDRB,adjust.method="fdr", p.value = 0.05,lfc=0)
label_case <-  "_DRB_vs_noDRB_fdr0.05"
#summary(res_DRB_vs_noDRB)

summary_res_DRB_vs_noDRB <- summary(res_DRB_vs_noDRB)

tss_res_PDS <- export_bed_and_TSS_in_bed_regions(tss_file,'./',label_case,res_DRB_vs_noDRB)

pdf(paste0('MD_MA_plot_differential_binding_analysis_',label_case,'.pdf'),8,7)
plotMD(qlf.DRB_vs_noDRB,status = res_DRB_vs_noDRB,values =c(1,0,-1),col=c("#ff000080","#00000030","#0000ff80"),pch = 19,cex=c(1,0.8,1))
title(paste0("\n\nDEup:",length(which(res_DRB_vs_noDRB==1)),
             " - DEdown: ",length(which(res_DRB_vs_noDRB==-1)),
             " - notDE: ",length(which(res_DRB_vs_noDRB==0))," - Tot: ",length(res_DRB_vs_noDRB)))
dev.off()

message("-------------------------- .. -----------------------------")
message('fdr=0.05')
print(cbind(summary(res_DRB_vs_noDRB), summary(res_DRB_vs_noDRB)/sum(summary(res_DRB_vs_noDRB))))
res_DRB_vs_noDRB_fdr0.05 <- res_DRB_vs_noDRB
message("-------------------------- .. -----------------------------")


res_DRB_vs_noDRB <- decideTestsDGE(qlf.DRB_vs_noDRB,adjust.method="none", p.value = 0.01,lfc=0)
label_case <-  "_DRB_vs_noDRB_pval0.01"
#summary(res_DRB_vs_noDRB)

summary_res_DRB_vs_noDRB <- summary(res_DRB_vs_noDRB)

tss_res_PDS <- export_bed_and_TSS_in_bed_regions(tss_file,'./',label_case,res_DRB_vs_noDRB)

pdf(paste0('MD_MA_plot_differential_binding_analysis_',label_case,'.pdf'),8,7)
plotMD(qlf.DRB_vs_noDRB,status = res_DRB_vs_noDRB,values =c(1,0,-1),col=c("#ff000080","#00000030","#0000ff80"),pch = 19,cex=c(1,0.8,1))
title(paste0("\n\nDEup:",length(which(res_DRB_vs_noDRB==1)),
             " - DEdown: ",length(which(res_DRB_vs_noDRB==-1)),
             " - notDE: ",length(which(res_DRB_vs_noDRB==0))," - Tot: ",length(res_DRB_vs_noDRB)))
dev.off()


message("-------------------------- .. -----------------------------")
message('pval=0.01')
print(cbind(summary(res_DRB_vs_noDRB), summary(res_DRB_vs_noDRB)/sum(summary(res_DRB_vs_noDRB))))
res_DRB_vs_noDRB_pval0.01 <- res_DRB_vs_noDRB
message("-------------------------- .. -----------------------------")

res_DRB_vs_noDRB <- decideTestsDGE(qlf.DRB_vs_noDRB,adjust.method="none", p.value = 0.05,lfc=0)
label_case <-  "_DRB_vs_noDRB_pval0.05"
#summary(res_DRB_vs_noDRB)

summary_res_DRB_vs_noDRB <- summary(res_DRB_vs_noDRB)

tss_res_PDS <- export_bed_and_TSS_in_bed_regions(tss_file,'./',label_case,res_DRB_vs_noDRB)

pdf(paste0('MD_MA_plot_differential_binding_analysis_',label_case,'.pdf'),8,7)
plotMD(qlf.DRB_vs_noDRB,status = res_DRB_vs_noDRB,values =c(1,0,-1),col=c("#ff000080","#00000030","#0000ff80"),pch = 19,cex=c(1,0.8,1))
title(paste0("\n\nDEup:",length(which(res_DRB_vs_noDRB==1)),
             " - DEdown: ",length(which(res_DRB_vs_noDRB==-1)),
             " - notDE: ",length(which(res_DRB_vs_noDRB==0))," - Tot: ",length(res_DRB_vs_noDRB)))
dev.off()

message("-------------------------- .. -----------------------------")
message('pval=0.05')
print(cbind(summary(res_DRB_vs_noDRB), summary(res_DRB_vs_noDRB)/sum(summary(res_DRB_vs_noDRB))))
res_DRB_vs_noDRB_pval0.05 <- res_DRB_vs_noDRB
message("-------------------------- .. -----------------------------")



## summary_numbers


SUMM <- cbind(summary(res_DRB_vs_noDRB_fdr0.01),summary(res_DRB_vs_noDRB_fdr0.05),summary(res_DRB_vs_noDRB_pval0.01),summary(res_DRB_vs_noDRB_pval0.05))
print(SUMM)
colnames(SUMM) <- c('DRB_vs_noDRB_fdr0.01','DRB_vs_noDRB_fdr0.05','DRB_vs_noDRB_pval0.01',
                    'DRB_vs_noDRB_pval0.05')
write.table(SUMM,file='summary_all_cases_DBA.DRB.csv',quote = F, sep = ",", col.names = NA)

