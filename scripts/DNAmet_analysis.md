DNAmet_analysis
================
M, Marcao
2024-06-20

.ppt presentation at
/Users/mmarcao/Library/CloudStorage/OneDrive-St.JudeChildren’sResearchHospital/Jasmine_group/LarissaFurtado_DNAmet_RNAseq_integration

Source to Differential Region methylation (DMR) at \# Source at
<https://github.com/hamidghaedi/Methylation_Analysis?tab=readme-ov-file>

# Process IDAT files into beta-value matrix

``` r
# Load packages ---------
library(BiocManager)
require(IlluminaHumanMethylationEPICmanifest)
require(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) 
library(shinyMethyl)
library(bumphunter)
library(minfi)
library(dplyr)

### PROCESSING IDAT FILES  ------------------
# Load metadata
load("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/met_meta.rda")
met # Felipe's DNAmet matrix
head(meta) # Metadata/labsheet
dim(meta)
# Creating group information
met$Group <- c(rep('ETV6_NTRK3_fused', 22), rep('Kinase_fused', 16))
target <- data.frame(sample = colnames(met), group = met$Group, source = 'StJude')
target$source[c(10:22)] <- 'PublicData'
met <- NULL # Delete it. I'm processing it on my own

# 1. Load IDAT files -----------
# ETV6_NTRK3 samples
ETV6_NTRK3_idat_dir <- '/media/ResearchHome/plummgrp/home/common/LFurtado-colab/idat_files_ETV6-NTRK3_fused_tumors'
RGC_data_etv6 <- read.metharray.exp(base = ETV6_NTRK3_idat_dir, targets = NULL, force = TRUE, recursive = T) 

summary.idat <- shinySummarize(RGC_data_etv6)
runShinyMethyl(summary.idat) 

# kinase samples
kinase_idat_dir <- '/media/ResearchHome/plummgrp/home/common/LFurtado-colab/idat_files_kinase_fused_tumors'
RGC_data_kinase <- read.metharray.exp(base = kinase_idat_dir, targets = NULL, force = TRUE, recursive = T) 

summary.idat <- shinySummarize(RGC_data_kinase)
runShinyMethyl(summary.idat) 

# GEO samples (GSE140686)
geo_idat_dir <- '/media/ResearchHome/plummgrp/home/common/LFurtado-colab/idat_files_GSE140686'
RGC_data_geo <- read.metharray.exp(base = geo_idat_dir, targets = NULL, force = TRUE, recursive = T) 

summary.idat <- shinySummarize(RGC_data_geo)
runShinyMethyl(summary.idat) 
summary.idat@sampleNames

# all samples 
all_idat_dir <- '/media/ResearchHome/plummgrp/home/common/LFurtado-colab'
RGC_data_all <- read.metharray.exp(base = all_idat_dir, targets = NULL, force = TRUE, recursive = T) 

summary.idat <- shinySummarize(RGC_data_all)
runShinyMethyl(summary.idat) 


# 2. pvalue detection
detP <- detectionP(RGC_data_all, type = "m+u") 
table(detP > 0.05)


# 3. Preprocess the data
proc_data <- preprocessRaw(RGC_data_all) 


# 4. Mask probes that failed p-value detection
proc_data_r <- ratioConvert(proc_data)
is.na(assays(proc_data_r)$Beta) <- (detP[rownames(proc_data_r), colnames(proc_data_r)] > 0.05)
beta <- getBeta(proc_data_r)
head(beta)
dim(beta) #865859     38

# 5. Remove mask probes
load("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/hm450.anno.Rda") # got it from Tathi # this is 450k but later I'm using EPIC annotation (integration, DM CpGs-to-genes etc)
probes_remove <- subset(hm450.anno, chrm_A %in% c("chrX","chrY", "chrM") & MASK_general == FALSE)$probeID 
beta <- beta[!rownames(beta) %in% probes_remove, ]
dim(beta) #856801     38
head(beta)
beta_meta <- target
# save(beta, beta_meta, file = '/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/DNAmet_V2/beta_and_meta_data.rda')


# QC PLOTS ------------------
qc <- getQC(proc_data)
plotQC(qc) # it gets good with the preprocessRaw() method
densityPlot(proc_data, sampGroups = target$group)
densityBeanPlot(proc_data, sampGroups = target$group)
controlStripPlot(RGC_data_all, controls="BISULFITE CONVERSION II")
```

# Differential Methylated Positions - For Clustering/Classifying Purpose (will be used on heatmap and correlation)

``` r
load('/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/DNAmet_V2/beta_and_meta_data.rda')
beta %>% dim() #856801     38
beta <- as.data.frame(beta)
beta <- beta[, colnames(beta) %in% beta_meta$sample]
identical(colnames(beta), beta_meta$sample) #TRUE
beta <- na.omit(beta)
dim(beta) #836811     38
# Run dmpFinder (method 1) -----
dmp_allspm <- dmpFinder(as.matrix(beta), pheno = beta_meta$group  , type = "categorical")
table(dmp_allspm$qval < 0.05) #99234 CpGs statistic significant  
dmp_sig <- dmp_allspm[dmp_allspm$qval < 0.05, ]
dmp_sig$CpG <- rownames(dmp_sig)
dmp_sig_method1 <- dmp_sig

# Wilcoxon test (method 2) ---------
identical(colnames(beta), beta_meta$sample) #TRUE
# Ordering the samples by group of comparison in beta
etv6 <- beta_meta[beta_meta$group %in% 'ETV6_NTRK3_fused', ]$sample
kinase <- beta_meta[beta_meta$group %in% 'Kinase_fused', ]$sample
etv6_mt <- beta[, colnames(beta) %in% etv6]
length(colnames(etv6_mt)) #22 samples
kinase_mt <- beta[, colnames(beta) %in% kinase]
length(colnames(kinase_mt)) #16 samples
beta <- cbind(kinase_mt, etv6_mt)

# Run Differential Methylation 
require(exactRankTests)
require(parallel)
values <- t(beta) #transpose beta
values <- data.frame(values)
# parallel processing ON for wilcoxon test
wpvalues <- unlist(mclapply(values,
                            function(CpG) {
                              zz <- wilcox.exact(CpG[1:  dim(kinase_mt)[2]],
                                                 CpG[c(dim(kinase_mt)[2]+1) : dim(beta)[2]], exact=T) 
                              z <- zz$p.value
                              return(z)
                            }, mc.cores= 20)) # set n of cores
# adjust pvalue 
wpvalues_adj <- p.adjust(wpvalues, method = "BH")
wpvalues_adj <- data.frame(wpvalues_adj)
wpvalues_adj$CpG <- rownames(wpvalues_adj)
hist(wpvalues_adj$wpvalues_adj)
table(wpvalues_adj$wpvalues_adj < 0.05) #75795 CpGs statistic significant 
wpvalues_adj_sig <- wpvalues_adj[wpvalues_adj$wpvalues_adj < 0.05, ]
dmp_sig_method2 <- wpvalues_adj_sig

# Merge both DMP results
DMP_CpGs <- merge(dmp_sig_method1, dmp_sig_method2, by = 'CpG')
dim(DMP_CpGs) #71094 CpGs statistic significant 
# saveRDS(dmp_sig_merged, '/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/DNAmet_V2/dmp_sig_merged.rds')


# Calculate Mean Differential DNA methylation ----------
DMP_CpGs <- readRDS('/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/DNAmet_V2/dmp_sig_merged.rds')
load('/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/DNAmet_V2/beta_and_meta_data.rda')

calculate_diff_mean <- function(data, 
                                metadata, 
                                condition_one, 
                                threshold) {
  
  list_diffmean_dfs <- list()
  data <- as.data.frame(data)
  condition_all_but_one <- group_names[!group_names %in% condition_one]
  ###condition_all_but_one <- 'Mesenchymal-like' #TESTING only. 
  for(i in 1:length(condition_all_but_one)) {
    print(paste0("Processing group: ", condition_all_but_one[i]))
    
    try({ #it ignores any error that might prevent the code of going forward
      #I did it because of the absence of hyper or hypo probes were making the loop stop
      loop_data <- data
      condition_1_sampleID <- metadata[metadata$group %in% condition_one, ]$sample
      condition_2_sampleID <- metadata[metadata$group %in% condition_all_but_one[i], ]$sample
      
      print(paste0("Length of condition_1_sampleID: ", length(condition_1_sampleID)))
      print(paste0("Length of condition_2_sampleID: ", length(condition_2_sampleID)))
      
      loop_data$meanM1 <- apply(loop_data[, condition_1_sampleID], 1, mean, na.rm = TRUE)
      loop_data$meanM2 <- apply(loop_data[, condition_2_sampleID], 1, mean, na.rm = TRUE)
      loop_data$DiffMean <- loop_data$meanM1 - loop_data$meanM2
      loop_data$Comparison <- paste0(condition_one, '_', condition_all_but_one[i], '_', condition_one, '_', 'Orientation')
      
      
      
      loop_data$DNAmet_orientation <- NA
      #loop_data$DNAmet_orientation <- as.character(loop_data$DNAmet_orientation) #this time it's not necessary
      if(dim(loop_data[loop_data$DiffMean > threshold, ])[1] > 0) {
        loop_data[loop_data$DiffMean > threshold, ]$DNAmet_orientation <- 'hyper'
      } else {
        # do nothing
      }
      
      if(dim(loop_data[loop_data$DiffMean < -threshold, ])[1] > 0) {
        loop_data[loop_data$DiffMean < -threshold, ]$DNAmet_orientation <- 'hypo'
      } else {
        # do nothing
      }
      
      loop_data[loop_data$DNAmet_orientation %in% NA, ]$DNAmet_orientation <- 'not_diff'
      
      
      
      loop_data$probeID <- rownames(loop_data)
      
      list_diffmean_dfs[[i]] <- loop_data[, c('DiffMean', 'Comparison', 'DNAmet_orientation', 'probeID')]
      print(paste0(list_diffmean_dfs[[i]]$Comparison[1], ' has been stored into the list.'))
    },  silent = FALSE)
  }
  return(list_diffmean_dfs) 
}

# Filtering beta to only sig. CpGs 
beta <- beta[rownames(beta) %in% DMP_CpGs$CpG, ]
dim(beta) #71094    38
group_names <- names(table(beta_meta$group))
list_diffmean_dfs <- calculate_diff_mean(data = beta,
                                        metadata = metadata,
                                        condition_one = group_names[2], #"kinese_fused_tumors" is the direction 
                                        threshold = 0.3)

dmp_diffmean <- do.call('rbind', list_diffmean_dfs)
head(dmp_diffmean)
length(names(table(dmp_diffmean$Comparison))) # it should be equal to the number if groups you are comparing - 1

names(DMP_CpGs)[names(DMP_CpGs) == "CpG"] <- "probeID"
DMP_CpGs <- DMP_CpGs[, c('probeID', 'wpvalues_adj')]
dmp_diffmean <- merge(dmp_diffmean, DMP_CpGs, by = 'probeID')
# saveRDS(dmp_diffmean, '/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/DNAmet_V2/dmp_diffmean.rds')
```

# Differential Methylated Regions - For Biological Relations (will be used on RNA-DNA data integration and CpG-to-promoter analysis)

``` r
# Load packages ------
library(dplyr)
library(minfi)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DMRcate")
library(DMRcate)
library(Gviz)
library(RColorBrewer)
library(GenomicRanges)
library(rtracklayer)
library(HelpersMG)
library(data.table)


# Step 1. Generate DMRs table (DNAmet orientation, CpGs regions, n of CpGs, CpGs probe_ID etc) -----  ----- ----- ----- ----- ----- ----- ----- ----- -----
# Load data  -------------
load("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/beta_and_meta_data.rda")
head(beta)
dim(beta) #856801     38
table(is.na(beta))
# FALSE     TRUE 
# 32506123    52315
beta <- na.omit(beta)
dim(beta) #36811     38
head(beta_meta)
dim(beta_meta) #38  3
rownames(beta_meta) <- beta_meta$sample

identical(rownames(beta_meta), colnames(beta)) #TRUE

# Prepare metadata as "design" for the DMR
beta_meta$group <- as.factor(beta_meta$group)
beta_meta$group <- relevel(beta_meta$group, ref = "ETV6_NTRK3_fused") # this set Kinase_fused as the direction of DNAmeth  
design <- model.matrix(~ group, data = beta_meta)
# setting some annotation
myAnnotation <- cpg.annotate(object = beta, datatype = "array", 
                             what = "Beta", 
                             analysis.type = "differential", 
                             design = design, 
                             contrasts = FALSE, 
                             coef = 2, # groupKinase_fused design column - the orientation we want to see (same orientation we set for RNAseq)
                             arraytype = "EPIC",
                             fdr = 0.001)
str(myAnnotation)

# DMR analysis
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2) #pcutoff = defaut
results.ranges <- extractRanges(DMRs, genome = 'hg38')
print(results.ranges) # DMR results 
dmr.table <- as.data.frame(results.ranges)

### Saving
# saveRDS(dmr.table, '/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/Round_2/reviewed_DM_direction/dmr_table.rds')
# dmr.table <- readRDS("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/Round_2/reviewed_DM_direction/dmr_table.rds")


### Retrieving CpGs from DMRs ----------
# Turn illumina manifest (CpG-gene annotation) into GRange object
# So we can get info. from indiviudal CpGs and gene annotations from the manifest
library(readr)
EPIC.hg38.anno <- read_csv("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/Round_2/infinium-methylationepic-v-1-0-b5-manifest-file.csv", skip = 7)
# Make it a genomic range object
EPIC.hg38.anno.gr <- makeGRangesFromDataFrame(EPIC.hg38.anno, keep.extra.columns=T, start.field = "Start_hg38", end.field = "End_hg38", seqnames.field = "CHR_hg38", strand.field="Strand_hg38", na.rm=TRUE) 

# Components from DMR 
chr <- dmr.table$seqnames
start <- dmr.table$start
end <- dmr.table$end
# Merging illumina manifest annotation to DMR output
filtered_gr <- subsetByOverlaps(EPIC.hg38.anno.gr, GRanges(seqnames = chr, ranges = IRanges(start = start, end = end), keep.extra.columns = TRUE))
df <- as.data.frame(filtered_gr)
length(df$Name) #426 CpGs from DMR output
df$UCSC_RefGene_Name
df_CpGs_DMR <- data.frame(df$Name)

### Saving 
# saveRDS(df_CpGs_DMR, '/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/Round_2/df_CpGs_DMR.rds')
```

# Volcano Plot (only DMP CpGs used)

``` r
### VOLCANO PLOT ---------
dmp_diffmean <- readRDS('/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/DNAmet_V2/dmp_diffmean.rds')
load('/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/DNAmet_V2/beta_and_meta_data.rda')


dim(beta) # 856801     38
beta <- as.data.frame(beta)
beta$probeID <- rownames(beta)
volcano <- merge(dmp_diffmean, beta, by = "probeID")
dim(volcano)
table(volcano$DNAmet_orientation)
# hyper   hypo  not_diff 
# 75      580    70439

png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/volcano_DMP_DNAmet.png", width = 8, height = 6, units = "in", res = 300)

library(ggplot2); theme_set(theme_classic()) 
ggplot(data=volcano[, ], aes(x=DiffMean, y=-1*log10(wpvalues_adj), colour=DNAmet_orientation)) +
  geom_point() +
  xlab("Diff Mean Methylation") + ylab("-1 * log10 Significance") + 
  scale_color_manual(breaks=c("hyper","hypo","not_diff"), # color scale (for points) 
                     values=c("red","blue","black"),
                     labels=c('Hyper methylated',"Hypo methylated","< |0.3| Diff Mean"),
                     name="Legend")  +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 20),  
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 16, color = "black", face = "bold")) +
  ggtitle("Kinase_fused vs ETV6_NTRK3_fused")

dev.off()
```

# PCA/tSNE Plot (whole array - 850k)

``` r
# PCA/tSNE plot ------------
# Load DNAmet matrix and preparing metadata 
dmp_diffmean <- readRDS('/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/DNAmet_V2/dmp_diffmean.rds')
DMP_CpGs_sig <- dmp_diffmean[dmp_diffmean$DNAmet_orientation %in% c("hyper", "hypo"), ]$probeID
load('/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/DNAmet_V2/beta_and_meta_data.rda')
beta[1:3,1:3]
beta_meta[1:3,]
names(beta_meta)[names(beta_meta) == "sample"] <- "smp_ID"
# Load $smp_type information about stjude patients - I got it from Larissa's spreadsheet at https://sjcrh-my.sharepoint.com/:x:/r/personal/mmarcao_stjude_org1/_layouts/15/Doc.aspx?sourcedoc=%7BE754BE85-9809-4FF6-B5B3-B6E6D0997BBF%7D&file=Metadata_lvf_1-31-24_maycon.xlsx&action=default&mobileredirect=true&DefaultItemOpen=1&web=1

df <- data.frame(
  smp_ID = c("207558610121_R01C01", NA, "205715820168_R06C01", "207495040012_R05C01", NA, "206467000107_R05C01", "207558610121_R02C01", "207558610121_R03C01", "207558610121_R04C01", "207558610121_R05C01", NA, "207558610121_R06C01", "207558610121_R07C01", "207179240158_R05C01", NA, "207558610121_R08C01", "207179240091_R05C01", "207343240018_R08C01", "206154070072_R03C01", "207558610086_R02C01", "206250280180_R02C01", NA, "207339140028_R01C01", "207558610127_R01C01", "207558610127_R02C01", "207558610127_R03C01", "207558610127_R04C01", "207558610127_R05C01", "207558610127_R06C01", "207558610127_R07C01", "201869680190_R01C01", "201332340140_R08C01", "200848860136_R01C01", "200848860134_R08C01", "200928190007_R08C01", "200928190007_R07C01", "200928190007_R06C01", "200928190007_R04C01", "200928190007_R03C01", "200928190007_R02C01", "200928190007_R01C01", "200925700210_R06C01", "200848860099_R05C01", "9553932008_R01C02"),
  
  fusion = c("NRF1::BRAF", "PLEKHH2::ALK", "PLEKHH2::ALK", "TPR::NTRK1", "TPM3::NTRK1", "TPM3::NTRK1", "TPR::NTRK1", "LMNA::NTRK1", "LMNA::NTRK1", "PDE4DIP::NTRK1", "PDE4DIP::NTRK1", "EML4::NTRK3", "RBPMS::NTRK3", "RCC1::ALK", "EML4::NTRK3", "TTYH3::BRAF", "TPM3::NTRK1", "EML4::ALK\nOutside lab", "TPM3::NTRK1\nOutside lab", "TPM3::NTRK1", "ETV6:NTRK3", "ETV6:NTRK3", "ETV6:NTRK3", "ETV6::NTRK3", "ETV6::NTRK3", "ETV6::NTRK3", "ETV6::NTRK3", "ETV6::NTRK3", "ETV6::NTRK3", "ETV6::NTRK3", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  
  histologic = c("Spindle cell neoplasm with NRF1::BRAF rearrangement", "Spindle cell neoplasm", "Spindle cell neoplasm", "Spindle cell neoplasm with evidence of TPR::NTRK1 rearrangement", "Spindle cell sarcoma with NTRK1 gene rearrangement", "NTRK-rearranged spindle cell neoplasm", "malignant spindle cell neoplasm with TPR::NTRK1 fusion, most in keeping with infantile fibrosarcoma", "low-grade sarcoma with LMNA-NTRK1 fusion and homozygous deletion of CDKN2A,", "Lipofibromatosis", "Recurrent/residual variant of infantile fibrosarcoma", "Infantile fibrosarcoma, variant with NTRK1 fusion, status post chemotherapy", "Infantile fibrosarcoma-like tumor", "Congenital NTRK-rearranged spindle cell neoplasm of the gastrointestinal tract with evidence of RBPMS::NTRK3 fusion transcript", "Spindle cell neoplasm with ALK gene rearrangement", "Low-grade fibroblastic neoplasm with neural differentiation suggestive of NTRK1-associated mesenchymal tumor.", "suspicious for infantile fibrosarcoma", "Spindle cell neoplasm with TPM3::NTRK1 fusion, most c/w infantile fibrosarcoma", "Spindle cell neoplasm with EML4::ALK fusion", "Spindled cell neoplasm with TPM3-NTRK1 fusion and meningioangiomatosis-like growth pattern;", "NTRK-rearranged sarcoma with evidence of TPM3::NTRK1 fusion", "IFS", "IFS", "IFS", "IFS", "IFS", "IFS", "IFS", "IFS", "IFS", "IFS", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  
  tumor_type = c("kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "kinase-fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor", "ETV6-NTRK3 fused tumor"),
  
  smp_type = c("FFPE", "FFPE", "Fresh/frozen", "FFPE", "FFPE", "FFPE", "FFPE", "FFPE", "FFPE", "Fresh/frozen", "Fresh/frozen", "FFPE", "FFPE", "FFPE", "FFPE", "FFPE", "FFPE", NA, NA, "FFPE", "Fresh/frozen", "Fresh/frozen", "Fresh/frozen", "Fresh/frozen", "FFPE", "FFPE", "FFPE", "FFPE", "FFPE", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  
  stringsAsFactors = FALSE
)

# Review $fusion labels 
beta_meta_st <- merge(beta_meta, df, by="smp_ID")
beta_meta_st[beta_meta_st$fusion %in% "ETV6:NTRK3", ]$fusion <- "ETV6::NTRK3"
beta_meta_st[beta_meta_st$fusion %in% "TPM3::NTRK1\nOutside lab", ]$fusion <- "TPM3::NTRK1"
beta_meta_st[beta_meta_st$fusion %in% "EML4::ALK\nOutside lab", ]$fusion <- "EML4::ALK"
beta_meta_GEO <- beta_meta[!beta_meta$smp_ID %in% beta_meta_st$smp_ID, ]
dim(beta_meta_GEO)#13 3
beta_meta <- plyr::rbind.fill(beta_meta_st, beta_meta_GEO)
dim(beta_meta) #38  7


# Create three subsets out of beta and beta_meta
# 1. StJude + PublicData - they're the objets their selves 
dim(beta_meta) #38  7
dim(beta) # 856801     38
# 2. Stjude only 
beta_meta_st <- beta_meta[beta_meta$source %in% "StJude", ]
rownames(beta_meta_st) <- beta_meta_st$smp_ID
beta_st <- beta[, colnames(beta) %in% beta_meta_st$smp_ID]
beta_meta_st <- beta_meta_st[colnames(beta_st), ]
identical(rownames(beta_meta_st), colnames(beta_st)) # TRUE
# 3. ETV6-StJude + PublicData (there's only ETV6 on PublicData)
beta_meta_ETV6 <- beta_meta[beta_meta$group %in% "ETV6_NTRK3_fused", ]
rownames(beta_meta_ETV6) <- beta_meta_ETV6$smp_ID
beta_ETV6 <- beta[, colnames(beta) %in% beta_meta_ETV6$smp_ID]
beta_meta_ETV6 <- beta_meta_ETV6[colnames(beta_ETV6), ]
identical(rownames(beta_meta_ETV6), colnames(beta_ETV6)) # TRUE


# Stjude only - PCA/tSNE  -----------
beta_st <- na.omit(beta_st)
pca <- prcomp(t(beta_st)) 
aux <- as.data.frame(pca$x[, 1:3]) 
scores <- merge(beta_meta_st, aux, by.y=0, by.x="smp_ID", all.x=T)
png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/PCA_WholeArray_st_only_DNAmet.png", width = 10, height = 6, units = "in", res = 300)
# Plot PCA
library(ggplot2); theme_set(theme_classic())
ggplot(scores, aes(x=PC1, y=PC2, colour=factor(group), 
                   #shape = smp_type
                   )) +
  geom_point(size = 4) +
  scale_color_manual(values=c('slateblue3', 'wheat3'), name="Group") +
  xlab(paste0("PC1 (", prettyNum(summary(pca)$importance[2,1]*100, digits = 2), "%)")) +
  ylab(paste0("PC2 (", prettyNum(summary(pca)$importance[2,2]*100, digits = 2), "%)")) +
  scale_x_continuous(labels = scales::scientific_format()) +
  scale_y_continuous(labels = scales::scientific_format()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16),  
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, color = "black", face = "bold")) +
  #geom_text_repel(aes(label = smpID)) +
  ggtitle("Whole Array - StJude samples only") 

dev.off()

library(Rtsne)
dim(beta_st)
set.seed(234234)
tsne_realData <- Rtsne(t(beta_st), perplexity=3, check_duplicates = FALSE) # #function to run t-sne
pdata.teste.tsne <- scores #renaming metadata
pdata.teste.tsne$tSNE1 <- tsne_realData$Y[,1]
pdata.teste.tsne$tSNE2 <- tsne_realData$Y[,2]
png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/tsne_WholeArray_st_only_DNAmet.png", width = 8, height = 6, units = "in", res = 300)
# Plot tSNE
library(ggplot2); theme_set(theme_classic())
ggplot(pdata.teste.tsne, aes(x=tSNE1, y=tSNE2, 
                             colour=factor(group), 
                             #shape = smp_type
                             )) +
  geom_point(size = 4) +
  scale_color_manual(values=c('slateblue3', 'wheat3'), name="Group") +
  # xlab(paste0("PC1 (", prettyNum(summary(pca)$importance[2,1]*100, digits = 2), "%)")) +
  # ylab(paste0("PC2 (", prettyNum(summary(pca)$importance[2,2]*100, digits = 2), "%)")) +
  scale_x_continuous(labels = scales::scientific_format()) +
  scale_y_continuous(labels = scales::scientific_format()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16),  
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, color = "black", face = "bold")) +
  #geom_text_repel(aes(label = smpID)) +
  ggtitle("Whole Array - StJude samples only") 

dev.off()



# Stjude + PublicData - PCA/tSNE -----------
beta <- na.omit(beta)
pca <- prcomp(t(beta)) 
aux <- as.data.frame(pca$x[, 1:3]) 
scores <- merge(beta_meta, aux, by.y=0, by.x="smp_ID", all.x=T)
png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/PCA_WholeArray_st_GEO_DNAmet.png", width = 10, height = 6, units = "in", res = 300)
# Plot PCA
library(ggplot2); theme_set(theme_classic())
ggplot(scores, aes(x=PC1, y=PC2, colour=factor(group), 
                   #shape = smp_type
                   )) +
  geom_point(size = 4) +
  scale_color_manual(values=c('slateblue3', 'wheat3'), name="Group") +
  xlab(paste0("PC1 (", prettyNum(summary(pca)$importance[2,1]*100, digits = 2), "%)")) +
  ylab(paste0("PC2 (", prettyNum(summary(pca)$importance[2,2]*100, digits = 2), "%)")) +
  scale_x_continuous(labels = scales::scientific_format()) +
  scale_y_continuous(labels = scales::scientific_format()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16),  
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, color = "black", face = "bold")) +
  #geom_text_repel(aes(label = smpID)) +
  ggtitle("Whole Array - StJude + PublicData")

dev.off()

library(Rtsne)
dim(beta_st)
set.seed(234234)
tsne_realData <- Rtsne(t(beta), perplexity=4, check_duplicates = FALSE) # #function to run t-sne
pdata.teste.tsne <- scores #renaming metadata
pdata.teste.tsne$tSNE1 <- tsne_realData$Y[,1]
pdata.teste.tsne$tSNE2 <- tsne_realData$Y[,2]
png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/tsne_WholeArray_st_GEO_DNAmet.png", width = 8, height = 6, units = "in", res = 300)
# Plot tSNE
library(ggplot2); theme_set(theme_classic())
ggplot(pdata.teste.tsne, aes(x=tSNE1, y=tSNE2, 
                             colour=factor(group), 
                             #shape = smp_type
                             )) +
  geom_point(size = 4) +
  scale_color_manual(values=c('slateblue3', 'wheat3'), name="Group") +
  # xlab(paste0("PC1 (", prettyNum(summary(pca)$importance[2,1]*100, digits = 2), "%)")) +
  # ylab(paste0("PC2 (", prettyNum(summary(pca)$importance[2,2]*100, digits = 2), "%)")) +
  scale_x_continuous(labels = scales::scientific_format()) +
  scale_y_continuous(labels = scales::scientific_format()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16),  
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, color = "black", face = "bold")) +
  #geom_text_repel(aes(label = smpID)) +
  ggtitle("Whole Array - StJude + PublicData") 

dev.off()

png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/tsne_WholeArray_st_GEO_source_DNAmet.png", width = 8, height = 6, units = "in", res = 300)
# Plot tSNE
library(ggplot2); theme_set(theme_classic())
ggplot(pdata.teste.tsne, aes(x=tSNE1, y=tSNE2, 
                             colour=factor(group), 
                             shape = source
                             )) +
  geom_point(size = 4) +
  scale_color_manual(values=c('slateblue3', 'wheat3'), name="Group") +
  # xlab(paste0("PC1 (", prettyNum(summary(pca)$importance[2,1]*100, digits = 2), "%)")) +
  # ylab(paste0("PC2 (", prettyNum(summary(pca)$importance[2,2]*100, digits = 2), "%)")) +
  scale_x_continuous(labels = scales::scientific_format()) +
  scale_y_continuous(labels = scales::scientific_format()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16),  
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, color = "black", face = "bold")) +
  #geom_text_repel(aes(label = smpID)) +
  ggtitle("Whole Array - StJude + PublicData") 

dev.off()


# ETV6-Stjude + ETV6-PublicData - PCA/tSNE (Supp. file probably) -----------
beta_ETV6 <- na.omit(beta_ETV6)
pca <- prcomp(t(beta_ETV6)) 
aux <- as.data.frame(pca$x[, 1:3]) 
scores <- merge(beta_meta_ETV6, aux, by.y=0, by.x="smp_ID", all.x=T)

scores$group3 <- scores$group
scores[scores$source %in% "StJude", ]$group3 <- "ST_ETV6"
scores[scores$source %in% "PublicData", ]$group3 <- "Public_ETV6"

png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/PCA_WholeArray_ETV6_st_GEO_DNAmet.png", width = 10, height = 6, units = "in", res = 300)
# Plot PCA
library(ggplot2); theme_set(theme_classic())
ggplot(scores, aes(x=PC1, y=PC2, colour=factor(group3))) +
  geom_point(size = 4) +
  #scale_color_manual(values=c('slateblue3', 'wheat3'), name="Group") +
  scale_color_manual(values=c('blue','red'), name="group3")+
  # xlab(paste0("PC1 (", prettyNum(summary(pca)$importance[2,1]*100, digits = 2), "%)")) +
  # ylab(paste0("PC2 (", prettyNum(summary(pca)$importance[2,2]*100, digits = 2), "%)")) +
  scale_x_continuous(labels = scales::scientific_format()) +
  scale_y_continuous(labels = scales::scientific_format()) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16),  
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, color = "black", face = "bold")) +
  #geom_text_repel(aes(label = smpID)) +
  ggtitle("Whole Array - ETV6 StJude and PublicData (Supp. file)") 

dev.off()
```

# Heatmap and correlation test (only DMP CpGs used)

``` r
### Heatmap  ------------
library(pheatmap)
library(RColorBrewer)
dim(beta) #
dim(beta_meta)#

my_sample_col <- data.frame(row.names = beta_meta$smp_ID, 
                            Group = beta_meta$group, 
                            #Histology = labsheet_hm$histologic_diagnosis,
                            Fusion = beta_meta$fusion) 

NTRK <- my_sample_col[grep('NTRK', my_sample_col$Fusion), ]
NTRK_butETV6 <- rownames(NTRK[!NTRK$Fusion %in% c('ETV6::NTRK3'), ])
NTRK_onlyETV6 <- rownames(NTRK[!rownames(NTRK) %in% NTRK_butETV6, ])

my_sample_col$Fusion_concat <- NA
my_sample_col[rownames(my_sample_col) %in% NTRK_butETV6, ]$Fusion_concat <- 'Others_NTRK'
my_sample_col[rownames(my_sample_col) %in% NTRK_onlyETV6, ]$Fusion_concat <- 'ETV6_NTRK3'

my_sample_col[my_sample_col$Fusion_concat %in% NA, ]$Fusion_concat[1:5] <- "Not_NTRK" 

my_sample_col[my_sample_col$Fusion_concat %in% NA, ]$Fusion_concat <- "NA"

my_sample_col[my_sample_col$Fusion %in% NA, ]$Fusion <- "NA"

my_sample_col$Source <- NA
my_sample_col[my_sample_col$Fusion %in% "NA", ]$Source <- "PublicData" 
my_sample_col[!my_sample_col$Fusion %in% "NA", ]$Source <- "StJude" 


ann_colors = list(
  Fusion_concat = c(ETV6_NTRK3 = 'black', Not_NTRK = 'blue', Others_NTRK = 'red', `NA` = "gray"))

png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/heatmap_DMP_st_GEO_DNAmet.png", width = 8, height = 8, units = "in", res = 300)
pheatmap(beta[rownames(beta) %in% DMP_CpGs_sig, ], 
         #annotation_row = my_probe_col, 
         annotation_col = my_sample_col[, ],
         show_rownames = FALSE,
         scale = "column",
         annotation_colors = ann_colors,
         main = "DMP CpGs (655) - Maybe we can still use this heatmap")

dev.off()


### Plot heatmap within only StJude samples 
# *** The Following Chunk of Code is not working ***
beta_meta_st <- beta_meta[beta_meta$source %in% "StJude", ]
my_sample_col_st <- data.frame(row.names = beta_meta_st$smp_ID,
                            Group = beta_meta_st$group,
                            #Histology = labsheet_hm$histologic_diagnosis,
                            Fusion = beta_meta_st$fusion)

NTRK <- my_sample_col_st[grep('NTRK', my_sample_col_st$Fusion), ]
NTRK_butETV6 <- rownames(NTRK[!NTRK$Fusion %in% c('ETV6::NTRK3'), ])
NTRK_onlyETV6 <- rownames(NTRK[!rownames(NTRK) %in% NTRK_butETV6, ])

my_sample_col_st$Fusion_concat <- NA
my_sample_col_st[rownames(my_sample_col_st) %in% NTRK_butETV6, ]$Fusion_concat <- 'Others_NTRK'
my_sample_col_st[rownames(my_sample_col_st) %in% NTRK_onlyETV6, ]$Fusion_concat <- 'ETV6_NTRK3'

my_sample_col_st[my_sample_col_st$Fusion_concat %in% NA, ]$Fusion_concat[1:5] <- "Not_NTRK"


ann_colors = list(
  Fusion_concat = c(ETV6_NTRK3 = 'black', Not_NTRK = 'blue', Others_NTRK = 'red'))

# my_sample_col_st <- my_sample_col_st[colnames(beta_st), ]
# identical(colnames(beta_st), rownames(my_sample_col_st)) # TRUE
png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/heatmap_DMP_st_only_DNAmet.png", width = 8, height = 8, units = "in", res = 300)
pheatmap(beta_st[rownames(beta_st) %in% DMP_CpGs_sig, ],
         #annotation_row = my_probe_col,
         annotation_col = my_sample_col_st,
         show_rownames = FALSE,
         scale = "column",
         annotation_colors = ann_colors,
         main = "DMP CpGs (655) Only StJude samples")

dev.off()


### Correlation test ------------
beta_DMP <- beta[DMP_CpGs$probeID, ]
library(ggplot2)
library(corrplot)
length(colnames(beta_DMP))
length(rownames(meta))
intersect(colnames(beta_DMP),
          rownames(meta))
identical(colnames(beta_DMP),
          rownames(meta)) #FALSE

meta <- meta[colnames(beta_DMP), ]
identical(colnames(beta_DMP),
          rownames(meta)) #TRUE

# colnames(beta_DMP) must be sample ID to calculate p-values 
cor_test_mat <- cor.mtest(beta_DMP, conf.level = 0.95)
cor_test_mat <- cor_test_mat$p

identical(colnames(beta_DMP), rownames(cor_test_mat)) #TRUE 
colnames(beta_DMP) <- meta$fusion_concat #doing that to be able to use fusion_concat as sample names. I've already checked and the samples order between cor_test_mat and beta_DMP are identical !

cor_matrix <- cor(beta_DMP)

df_color <- data.frame(smp_name = colnames(beta_DMP),
                       color = NA)

df_color$color <- c("blue", rep("black", 5), "blue", rep("red",4),
                    "black", "red", "red","blue", "blue", rep("red", 4),
                    "black", "blue", "red", "black", "black") # accordingly to the correlation plot sample display 
                    
smp_colors <- df_color$color


rownames(cor_test_mat)
png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/correlation_DMP_st_only_DNAmet.png", width = 8, height = 8, units = "in", res = 300)
corrplot(cor_matrix, 
         method = "color", 
         type = "full",
         tl.col = smp_colors, 
         tl.srt = 90, 
         order = "hclust",
         addrect = 2, 
         rect.col = 'blue', 
         rect.lwd = 3,
         p.mat = cor_test_mat, 
         sig.level = 0.001,  #insignificant values have an x sign 
         title = "DNAmet Correlation - StJude smps - DMP CpGs",
         mar = c(0,0,3,0),
         tl.cex = 1)
dev.off()
```

# Stemness score

``` r
### Stemness --------
# load("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/beta_and_meta_data.rda")
# head(beta_meta)
# head(beta)
load('/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/DNAmet_V2/beta_and_meta_data.rda')
names(beta_meta)[names(beta_meta) == "sample"] <- "smp_ID"
load("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/objects/DNAmet_stemness_model.Rda") 
w <- mm$w
w[1:5]
w_df <- as.data.frame(w)
w_df$probeID <- rownames(w_df)
X <- as.data.frame(beta)
X <- na.omit(X) 
# Make DNA mtx (X) and the stemness weigth (w_df) be the same dim and in the same order
X <- X[rownames(X) %in%  rownames(w_df) ,]
w_df <- w_df[rownames(X), ]
length(intersect(rownames(X),rownames(w_df))) # 206
identical(rownames(X),rownames(w_df)) # TRUE
X <- as.matrix(X)
dim(X) # 208  38
dim(w_df) #208   2
length(intersect(rownames(w_df),rownames(X))) # 206
identical(rownames(w_df),rownames(X)) #TRUE
w_df$probeID <- NULL 
# Apply stemness model into DNA mtx 
ss <- t(w_df) %*% X 
ss[1,1:3]
# Scale the scores into a ratio from 0 to 1 and store as data frame.
ss <- ss - min(ss)
ss <- ss / max(ss)
ss <- as.data.frame(t(ss))
colnames(ss) <- "Stemness_DNAmet" 
ss$smp_ID <- rownames(ss)

# Plot PCA/tSNE + Stemness
beta_meta <- merge(beta_meta, ss, by="smp_ID")
beta <- na.omit(beta)
pca <- prcomp(t(beta)) 
aux <- as.data.frame(pca$x[, 1:3]) 
scores <- merge(beta_meta, aux, by.y=0, by.x="smp_ID", all.x=T)
png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/PCA_WholeArray_stemness_st_GEO_DNAmet.png", width = 8, height = 8, units = "in", res = 300)
# PCA 
library(ggplot2); theme_set(theme_classic())
ggplot(scores, aes(PC1, PC2)) + geom_point(size=5, aes( fill=Stemness_DNAmet, color = group), pch = 21, stroke = 1.5) + 
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  scale_color_manual(values = c('slateblue3', 'wheat3')) +
  ylab(paste0('PC2 ', summary(pca)$importance[2, 2] * 100, '%')) +
  xlab(paste0('PC1 ', summary(pca)$importance[2, 1] * 100, '%') ) + theme_bw() +
  ggtitle('WholeArray - Stemness prediction')

dev.off()


library(Rtsne)
beta <- na.omit(beta)
set.seed(234234)
tsne_realData <- Rtsne(t(beta), perplexity=4, check_duplicates = FALSE) # #function to run t-sne
pdata.teste.tsne <- scores #renaming metadata
pdata.teste.tsne$tSNE1 <- tsne_realData$Y[,1]
pdata.teste.tsne$tSNE2 <- tsne_realData$Y[,2]
png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/tsne_WholeArray_stemness_st_only_DNAmet.png", width = 8, height = 6, units = "in", res = 300)
# tSNE
library(ggplot2); theme_set(theme_classic())
ggplot(pdata.teste.tsne, aes(tSNE1, tSNE2)) + geom_point(size=5, aes( fill=Stemness_DNAmet, color = group), pch = 21, stroke = 1.5) + 
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  scale_color_manual(values = c('slateblue3', 'wheat3')) +
  # ylab(paste0('PC2 ', summary(pca)$importance[2, 2] * 100, '%')) +
  # xlab(paste0('PC1 ', summary(pca)$importance[2, 1] * 100, '%') ) + theme_bw() +
  ggtitle('WholeArray - Stemness prediction')

dev.off()


png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/boxplot_stemness_st_GEO_DNAmet.png", width = 8, height = 8, units = "in", res = 300)
library(ggplot2); theme_set(theme_classic())
ggplot(scores[, ], aes(x=group, y=Stemness_DNAmet)) + 
  geom_boxplot(fill= c('slateblue3', 'wheat3'), 
    outlier.color = NA) + 
  geom_jitter (alpha=0.2)  +
  xlab("Groups") + 
  ylab("Stemness") + 
  theme(axis.title.x = element_text(colour="black", size = 12, vjust=-1), 
        axis.title.y = element_text(colour="black", size = 12), 
        axis.text.x = element_text(angle = 45, vjust= 1, size = 12, hjust = 1)) +
  #facet_wrap(~ sample_char) +
  ggtitle("Stemness diff. beetween groups (p = 0.01745)") 

dev.off()

t.test(scores[scores$group %in% "ETV6_NTRK3_fused", ]$Stemness_DNAmet,
       scores[scores$group %in% "Kinase_fused", ]$Stemness_DNAmet) # p-value = 0.01745

library(rstatix)
scores %>%
wilcox_test(Stemness_DNAmet ~ group, p.adjust.method = "none") %>%
  add_significance() #p-value = 0.0224
```

# Mapping Differential Methyalted Regions CpG-to-promoter

``` r
# DMR -------
# CpGs to promoters 
# only CpGs from dmr_table
df_CpGs_DMR <- readRDS("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/Round_2/reviewed_DM_direction/df_CpGs_DMR.rds")
names(df_CpGs_DMR) <- "CpG_ID"
# dmr_table output (not a Genomic Range obj)
dmr_table <- readRDS("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/Round_2/reviewed_DM_direction/dmr_table.rds")
# transforming both back in a Genomic Range obj to retrieve all dmr_table columns but with the specific CpGs in each genomic interval 
library(readr)
EPIC.hg38.anno <- read_csv("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/Round_2/infinium-methylationepic-v-1-0-b5-manifest-file.csv", skip = 7)
EPIC.hg38.anno <- as.data.frame(EPIC.hg38.anno)
df_CpGs_DMR <- EPIC.hg38.anno[EPIC.hg38.anno$Name %in% df_CpGs_DMR$CpG_ID, ]

### out of work vvvvvv
# Make both a genomic range object
# df_CpGs_DMR.anno.gr <- makeGRangesFromDataFrame(df_CpGs_DMR, keep.extra.columns=T, start.field = "Start_hg38", end.field = "End_hg38", seqnames.field = "CHR_hg38", strand.field="Strand_hg38", na.rm=TRUE) 
### out of work ^^^^^^^

# Reference Genome Annotation (illumina manifest annotation)
EPIC.hg38.anno.gr <- makeGRangesFromDataFrame(EPIC.hg38.anno, keep.extra.columns=T, start.field = "Start_hg38", end.field = "End_hg38", seqnames.field = "CHR_hg38", strand.field="Strand_hg38", na.rm=TRUE) 

# DMR output 
dmr_table.anno.gr <- makeGRangesFromDataFrame(dmr_table, keep.extra.columns=T, start.field = "start", end.field = "end", seqnames.field = "seqnames", strand.field="strand", na.rm=TRUE)

# Retrieving DMR CpGs - because in DMR ouput we don't have "probe_ID" information + Mapping CpGs top promoter regions 
library(ELMER)
promoter.gr <- get.feature.probe(promoter = TRUE, TSS.range = list(upstream = 2000, downstream = 500), rm.chr=c("chrX", "chrY"), met.platform = "EPIC")

# Interact to one genomic interval at a time
df_CpG_list <- list()
for(i in 1:length(dmr_table.anno.gr$no.cpgs)) {
  dmr_table_gr_1 <- dmr_table.anno.gr[i, ]
  # Merging them both - illumina manifest annotation to DMR output
  filtered_gr <- subsetByOverlaps(promoter.gr, dmr_table_gr_1)
  df_CpG <- as.data.frame(filtered_gr)
  if(nrow(df_CpG) != 0){ 
    df_dmr <- as.data.frame(dmr_table_gr_1)
    df_CpG$maxdiff <- NA
    df_CpG$meandiff <- NA
    df_CpG$meth_status <- NA
    df_CpG$maxdiff <- df_dmr$maxdiff
    df_CpG$meandiff <- df_dmr$meandiff
    
    table(df_CpG$meandiff > 0)
    if(any(df_CpG$meandiff > 0)){
      df_CpG[df_CpG$meandiff > 0,]$meth_status <- "hyper"
    }
    if(any(df_CpG$meandiff < 0)){
      df_CpG[df_CpG$meandiff < 0,]$meth_status <- "hypo"
    }
    df_CpG_list[[i]] <- df_CpG
  } 
  
}

df_CpG_promoter <- do.call(rbind, df_CpG_list) # basic way to get DMR CpGs + meth_status + promoter regions  


# Hyper meth CpGs 
hyper_CpG_promoter_df <- df_CpG_promoter[df_CpG_promoter$meth_status %in% "hyper", ]

df <- table(hyper_CpG_promoter_df$gene) %>% data.frame()
df <- df[order(df$Freq, decreasing = TRUE), ]

png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/hist_DMR_hyper_CpGtoPromoter_freq_st_GEO_DNAmet.png", width = 8, height = 8, units = "in", res = 300)
library(ggplot2); theme_set(theme_classic())
ggplot(df[, ], aes(reorder(x = factor(Var1), Freq), y = Freq, )) +
  geom_bar(stat = "identity", color = "black") +
  labs(#title = "Frequency of Cells by Sample ID",
    x = "Genes",
    y = "Frequency CpG on putative promoters") +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16),  
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(size = 16, color = "black", face = "bold")) +
  labs(fill = "") +
#facet_wrap(~ factor(Var1), ncol = 3) +

ggtitle("Freq of HYPER CpG from DMR and the promoter-genes annotation") 

dev.off()


# Hypo meth CpGs 
hypo_CpG_promoter_df <- df_CpG_promoter[df_CpG_promoter$meth_status %in% "hypo", ]

df <- table(hypo_CpG_promoter_df$gene) %>% data.frame()
df <- df[order(df$Freq, decreasing = TRUE), ]
png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/hist_DMR_hypo_CpGtoPromoter_freq_st_GEO_DNAmet.png", width = 8, height = 8, units = "in", res = 300)
library(ggplot2); theme_set(theme_classic())
ggplot(df[, ], aes(reorder(x = factor(Var1), Freq), y = Freq, )) +
  geom_bar(stat = "identity", color = "black") +
  labs(#title = "Frequency of Cells by Sample ID",
    x = "Genes",
    y = "Frequency CpG on putative promoters") +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16),  
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(size = 16, color = "black", face = "bold")) +
  labs(fill = "") +
  #facet_wrap(~ factor(Var1), ncol = 3) +
  
  ggtitle("Freq of HYPO CpG from DMR and the promoter-genes annotation") 

dev.off()

# Enrichment analysis from DMR-promoter-TFs 
list_chars <- list()
genes_concat <- hyper_CpG_promoter_df$gene
for(i in 1:length(genes_concat)) {
  chars <- strsplit(genes_concat, split = ";")[[i]]
  chars <- gsub(' ','', chars)
  list_chars[[i]] <- chars
}
genes_concat <- as.vector(do.call(rbind, list_chars))
genes_concat <- genes_concat[!genes_concat %in% NA]
DMR_prom_TFs_hyper <- unique(genes_concat)


list_chars <- list()
genes_concat <- hypo_CpG_promoter_df$gene
for(i in 1:length(genes_concat)) {
  chars <- strsplit(genes_concat, split = ";")[[i]]
  chars <- gsub(' ','', chars)
  list_chars[[i]] <- chars
  
}
genes_concat <- as.vector(do.call(rbind, list_chars))
genes_concat <- genes_concat[!genes_concat %in% NA]
DMR_prom_TFs_hypo <- unique(genes_concat)


library(clusterProfiler)
library(org.Hs.eg.db)

# Genes Ontology of annotated genes to Hyper methylated CpGs ------
ego2_hyper <- enrichGO(gene = DMR_prom_TFs_hyper,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 # ont  One of "BP", "MF", and "CC" subontologies,                                   # or "ALL" for all three
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

ego2_hyper@result[order(ego2_hyper@result$Count, decreasing = TRUE),][1:10,]
ego2_hyper@result[order(ego2_hyper@result$p.adjust, decreasing = FALSE),][1:10,]
# barplot(ego2_hyper, drop=TRUE, main = "")
#barplot(ego2_hyper, drop=FALSE, showCategory = 12, main = "")
#clusterProfiler::dotplot(ego2_hyper, showCategory=30) + ggtitle("")


# Plot Barplot independently to {barplot} function
png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/enrichmentBarplot_HyperMethylated_DMR_CpGtoPromoter.png", width = 8, height = 6, units = "in", res = 300)
pathways_to_plot <- ego2_hyper@result[ego2_hyper@result$p.adjust <= 0.05, ]
library(ggplot2); theme_set(theme_classic())
ggplot(pathways_to_plot, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
  geom_col() +
  labs(x = "Gene Count", y = "", title = "Kinase Hyper Methylated DMR CpGs-to-Promoters") +
  scale_fill_gradient(low = "red", high = "blue", name = "p.adjust")

dev.off()



# Genes Ontology of annotated genes to Hypo methylated CpGs -----
ego2_hypo <- enrichGO(gene = DMR_prom_TFs_hypo,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 # ont  One of "BP", "MF", and "CC" subontologies,                                   # or "ALL" for all three
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

ego2_hypo@result[order(ego2_hypo@result$Count, decreasing = TRUE),][1:10,]
ego2_hypo@result[order(ego2_hypo@result$p.adjust, decreasing = FALSE),][1:10,]

# barplot(ego2_hypo, drop=TRUE, main = "")

png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/enrichmentBarplot_Hypo_Methylated_DMR_CpGtoPromoter.png", width = 8, height = 6, units = "in", res = 300)
pathways_to_plot <- ego2_hypo@result[ego2_hypo@result$p.adjust <= 0.06, ]
library(ggplot2); theme_set(theme_classic())
ggplot(pathways_to_plot, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
  geom_col() +
  labs(x = "Gene Count", y = "", title = "Kinase Hypo Methylated DMR CpGs-to-Promoters") +
  scale_fill_gradient(low = "red", high = "blue", name = "p.adjust")

dev.off()


# write.csv(data.frame(DMR_prom_TFs_hyper), file = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/objects/DMR_prom_TFs_hyper.csv")
# write.csv(data.frame(DMR_prom_TFs_hypo), file = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/objects/DMR_prom_TFs_hypo.csv")
```
