RNAseq_analysis
================
M, Marcao
2024-06-18

\###.ppt presentation at
\###/Users/mmarcao/Library/CloudStorage/OneDrive-St.JudeChildren’sResearchHospital/Jasmine_group/LarissaFurtado_DNAmet_RNAseq_integration

### Align ETV6 smps

``` r
files <- dir("/media/ResearchHome/plummgrp/home/common/LFurtado-colab/RNA/FASTQ_files_ETV6_NTRK3_fused_tumors/", pattern = "fastq")
path = "/media/ResearchHome/plummgrp/home/common/LFurtado-colab/RNA/FASTQ_files_ETV6_NTRK3_fused_tumors/"
path.out = "/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/ETV6_smps/"
# #out.prefix = str_extract(files,"[^_]*_[^_]*")
# library(stringr)
# file_split_1 <- str_split_fixed(as.character(files), "[_][R][0-9][_][00]", 2)[,1]
# file_split_2 <-  str_split_fixed(as.character(file_split_1), "[_][L]", 2)[,1]
# unique_sample_ids <- unique(file_split_2)
# saveRDS(unique_sample_ids, '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/ETV6_smps/unique_sample_ids.rds')
# unique_sample_ids <- readRDS('/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/ETV6_smps/unique_sample_ids.rds')


unique_sample_ids <- c("SJBT032721_D1-S3", # done
                       "SJIFS031085_D1-S9", # done
                       "SJST030433_D1-S3", # done 
                       "SJST030567_D1-S10", # done
  "SJST032838_D1-S4", # done 
  "SJST032952_D1-S9", # NEED TO REVIEW THE FILE PAIR. MULTIPLE FASTQ.
  "SJST032952_D2-S5", # done
  "SJST032952_D4-S4", # done
  "SJST032952_D4-S9", # NEED TO REVIEW  - only one fastq
  "SJST033312_D1-S13", # done
  "SJST033312_D1-S2") # done - it's a way small file.BUT it's fine. It matches to the size of their fastq

# library(dplyr)
for(i in 1:length(unique_sample_ids)) {
  #i <- 1
  samples <- files[grep(unique_sample_ids[i], files)]
  samples_R1 <- samples[grep("R1_0", samples)]
  samples_R2 <- samples[grep("R2_0", samples)]
  
  if (length(samples) > 2) {
    # test code
    # print('4 fastqs')
  system(paste0("STAR --genomeDir /media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens/HS_STAR_dir_2 --runThreadN 20 --readFilesIn ",
         path, samples_R1[1],',',path,samples_R1[2],' ',path,samples_R2[1],',',path,samples_R2[2], " --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ",path.out, unique_sample_ids[i]))
  
} else { 
  # test code
  # print(paste0(length(samples), ' fastqs'))
  system((paste0("STAR --genomeDir /media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens/HS_STAR_dir_2 --runThreadN 20 --readFilesIn ",
                path, samples_R1[1],' ',path,samples_R2[1], " --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ",path.out, unique_sample_ids[i])))
}
}




# Check it if we ran all fastq samples ------------
# all the fastq samples
unique_sample_ids
# catch all .bam files 
library(stringr)
bam_out <- dir("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/ETV6_smps/", pattern = ".bam$")
file_split_1 <- str_split_fixed(as.character(bam_out), "[Aligned]", 2)[,1]
unique_sample_ids_check <- unique(file_split_1)
# compare them. If it's TRUE so you're all covered. All the samples has been aligned
identical(unique_sample_ids, unique_sample_ids_check) #FALSE 
# SJST032952_D1-S9 - NEED TO REVIEW THE FILE PAIR
# SJST032952_D4-S9 - NEED TO REVIEW - ONLY ONE FASTQ FILE. It has a .bam output but it's not right
```

### Align kinase smps

``` r
files <- dir("/media/ResearchHome/plummgrp/home/common/LFurtado-colab/RNA/FASTQ_files_kinase_fused_tumors/", pattern = "fastq")
path = "/media/ResearchHome/plummgrp/home/common/LFurtado-colab/RNA/FASTQ_files_kinase_fused_tumors/"
path.out = "/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/kinase_smps/"
#out.prefix = str_extract(files,"[^_]*_[^_]*")
unique_sample_ids <- c("SJST030375_R2-S1",
                       "SJST031362_D1-S21",
                       "SJST031362_D1-S8",
                       "SJST031690_D1-S1",
                       "SJST031792_D1-S19",
                       "SJST031920_D1-S3",
                       "SJST032767_D1-S1",
                       "SJST032767_D2-S8",
                       "SJST033308_D1-S8",
                       "SJST033389_D1-S3",
                       "SJST033491_D1-S2",
                       "SJST033791_D1-S3",
                       "SJST033835_D1-S10",
                       "SJST033983_D1-S6",
                       "SJST034036_D1-S5",
                       "SJST034397_D1-S6",
                       "SJST034534_D1-S3",
                       "SJST034815_D1-S3")



# library(dplyr)
for(i in 1:length(unique_sample_ids)) {
 
  samples <- files[grep(unique_sample_ids[i], files)]
  samples_R1 <- samples[grep("R1_0", samples)]
  samples_R2 <- samples[grep("R2_0", samples)]
  
  if (length(samples) > 2) {
    system(paste0("STAR --genomeDir /media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens/HS_STAR_dir_2 --runThreadN 20 --readFilesIn ", path, samples_R1[1],',',path,samples_R1[2],' ',path,samples_R2[1],',',path,samples_R2[2], " --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ",path.out, unique_sample_ids[i]))
    
  } else { system(paste0("STAR --genomeDir /media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens/HS_STAR_dir_2 --runThreadN 20 --readFilesIn ", path, samples_R1[1],' ', path, samples_R2[1], " --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ", path.out, unique_sample_ids[i]))
  }
}



# Check it if we ran all fastq samples ------------
# all the fastq samples
unique_sample_ids
# catch all .bam files 
library(stringr)
bam_out <- dir("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/kinase_smps/", pattern = ".bam$")
file_split_1 <- str_split_fixed(as.character(bam_out), "[Aligned]", 2)[,1]
unique_sample_ids_check <- unique(file_split_1)
# compare them. If it's TRUE so you're all covered. All the samples has been aligned
identical(unique_sample_ids, unique_sample_ids_check) #TRUE 
```

####### Obs: MultiQC screen shoots on ppt slides

### From .Bam to Count Matrix

``` r
# https://www.rdocumentation.org/packages/Rsubread/versions/1.22.2/topics/featureCounts
# package recommended by Felipe
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Rsubread")

library(Rsubread)
library(DESeq2)
library(sva)
library(dplyr)
library(stringr)

# ETV6 samples ---------------
# Define the paths to your BAM files and annotation file
annotation_file <- '/media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens/gencode.v19.chr_patch_hapl_scaff.annotation.gtf'
bam_files <- dir('/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/ETV6_smps/', patter = '.bam')
# Run featureCounts on the BAM files
setwd("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/ETV6_smps")
bam_files <- bam_files[-c(8)] # the SJST032952_D4-S9 is not properly aligned

start <- Sys.time()

counts <- featureCounts(files = bam_files[],
                        annot.ext = annotation_file,
                        isGTFAnnotationFile = TRUE,
                        strandSpecific = "2",
                        isPairedEnd=TRUE,
                        nthreads = 32)
end <- Sys.time()
print(end - start) # 5.331839 mins on nthreads = 32

# saveRDS(counts, '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/ETV6_counts.rds')

# Tidying count matrix up
ETV6_counts <- readRDS("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/ETV6_counts.rds")
head(ETV6_counts$annotation)
head(ETV6_counts$counts)
ETV6_counts <- as.matrix(ETV6_counts$counts)
library(stringr)
first_split <- str_split_fixed(as.character(colnames(ETV6_counts)), "[Aligned]", 2)[,1]
colnames(ETV6_counts) <- first_split
head(ETV6_counts)
ETV6_smp_ID <- colnames(ETV6_counts)
metadata_1 <- data.frame(ETV6_smp_ID, 
           rep('ETV6', length(ETV6_smp_ID)))

names(metadata_1) <- c('smpID', 'group')


# Kinase samples ---------------
# Define the paths to your BAM files and annotation file
annotation_file <- '/media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens/gencode.v19.chr_patch_hapl_scaff.annotation.gtf'
bam_files <- dir('/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/kinase_smps/', patter = '.bam')
# Run featureCounts on the BAM files
setwd("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/kinase_smps")

start <- Sys.time()

counts <- featureCounts(files = bam_files[], 
                        annot.ext = annotation_file, 
                        isGTFAnnotationFile = TRUE,
                        strandSpecific = "2",
                        isPairedEnd=TRUE,
                        nthreads = 32)
end <- Sys.time()
print(end - start) #8.907089 mins on nthreads = 32

# saveRDS(counts, '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/kinase_counts.rds')

# Tidying count matrix up
kinase_counts <- readRDS("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/kinase_counts.rds")
head(kinase_counts$annotation)
head(kinase_counts$counts)
kinase_counts <- as.matrix(kinase_counts$counts)
library(stringr)
first_split <- str_split_fixed(as.character(colnames(kinase_counts)), "[Aligned]", 2)[,1]
colnames(kinase_counts) <- first_split
head(kinase_counts)
kinase_smp_ID <- colnames(kinase_counts)
metadata_2 <- data.frame(kinase_smp_ID, 
                       rep('kinase', length(kinase_smp_ID)))

names(metadata_2) <- c('smpID', 'group')

# Combine metadata (there's only group info here)
metadata <- plyr::rbind.fill(metadata_1, metadata_2)
rownames(metadata) <- metadata$smpID
metadata$group <- plyr::mapvalues(metadata$group, from = c("ETV6", "kinase"), to = c("ETV6_NTRK3_fused", "Kinase_fused")) #changing group names
# saveRDS(metadata, '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/metadata.rds')

dim(ETV6_counts) #63568     9
dim(kinase_counts) #63568    18
identical(rownames(ETV6_counts), rownames(kinase_counts)) # TRUE
count_matrix <- cbind(ETV6_counts, kinase_counts)
head(count_matrix)

# saveRDS(count_matrix, '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/count_matrix.rds')
```

### Resolving sample IDs

``` r
library(stringr)
# Matching samples ID from aligment to Larissa's labsheet ------------------
# bam files/featureCounts to count matrix checkout --------------
# ETV6 samples
bam_files_ETV6 <- dir('/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/ETV6_smps/', patter = '.bam')
library(stringr)
first_split_ETV6 <- str_split_fixed(as.character(bam_files_ETV6), "[Aligned]", 2)[,1]
# bam_files_ETV6 <- str_split_fixed(as.character(first_split_ETV6), "[-]", 2)[,1]
# kinase samples
library(Rsubread)
bam_files_kinase <- dir('/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/kinase_smps/', patter = '.bam')
library(stringr)
first_split_kinase <- str_split_fixed(as.character(bam_files_kinase), "[Aligned]", 2)[,1]
# bam_files_kinase <- str_split_fixed(as.character(first_split_kinase), "[-]", 2)[,1]

# length(bam_files_ETV6) #10 samples
# length(bam_files_kinase) #18 samples
length(first_split_ETV6) #10
first_split_ETV6 <- first_split_ETV6[!first_split_ETV6 %in% 'SJST033312_D1-S13']
length(first_split_ETV6) #9
length(first_split_kinase)#18


# .bam files to count matrix
count_matrix <- readRDS("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/count_matrix.rds")
count_matrix <- as.data.frame(count_matrix)
dim(count_matrix) #27 samples
table(duplicated(colnames(count_matrix))) # no duplicates
count_matrix <- count_matrix[, !colnames(count_matrix) %in% 'SJST033312_D1-S13'] # multiqc fall out

# fixing the sample ID
# library(stringr)
# colnames(count_matrix) <- str_split_fixed(as.character(colnames(count_matrix)), "[-]", 2)[,1]
# from featureCounts/count matrix
# table(colnames(count_matrix))
# table(duplicated(colnames(count_matrix))) # 2 duplicates: SJST031362_D1 and SJST033312_D1

# from .bam files
# table(c(bam_files_ETV6,bam_files_kinase)) 
# table(duplicated(c(bam_files_ETV6,bam_files_kinase))) 
# 3 duplicates
# SJST032952_D4 (had only one fastq. So it started aligning but then stopped.), 
# SJST033312_D1 ()
# SJST031362_D1
table(duplicated(c(first_split_ETV6, first_split_kinase))) # no duplicates

setdiff(colnames(count_matrix), c(first_split_ETV6, first_split_kinase)) #nada
setdiff(c(first_split_ETV6, first_split_kinase), colnames(count_matrix)) #SJST032952_D4-S9
# SJST032952_D4-S9 has only one fastq. So it had been aligning but then stopped it. That's why we see this sample ID on the .bam files

#@@@ Conclusion: the count matrix is good. There're all sampleID as it should be. We don't need to change anything there.


# labsheet/metadata checkout -------------
library(readxl)
labsheet <- read_excel("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis/metadata_lvf_sub.xlsx")
labsheet <- data.frame(labsheet)
labsheet <- labsheet[!labsheet$fastq_ID %in% "NA", ]
names(labsheet)[names(labsheet) == 'fastq_ID'] <- 'smpID'
length(labsheet$smpID) #27 samples
table(duplicated(labsheet$smpID)) # no duplicates


# Take note: Labsheet is shorten on sample ID. eg: they have "SJST030375_D1" but there's still missing a piece of it (eg: SJST030375_D1_XX). So it might get things a bit confusing in the future analysis. 


# Match count matrix to labsheet -----------
# First, make their sampleID match
library(stringr)
colnames(count_matrix) <- str_split_fixed(as.character(colnames(count_matrix)), "[-]", 2)[,1]
# Comparing both sampleID
length(colnames(count_matrix)) #27 samples
table(duplicated(colnames(count_matrix))) # 1 duplicate
table(colnames(count_matrix)) # SJST031362_D1 
# Solve duplicate (SJST031362_D1) in count matrix 
colnames(count_matrix) <- make.names(colnames(count_matrix), unique = TRUE)
colnames(count_matrix)

length(labsheet$smpID)#27 samples
table(duplicated(labsheet$smpID)) # 0 duplicates now

setdiff(colnames(count_matrix), labsheet$smpID) #"SJST030375_R2" "SJST031362_D1.1"
setdiff(labsheet$smpID, colnames(count_matrix)) #"STST030375_R2" "SJST030375_D1" "SJST032952_D1"


#@@@ Conclusion: STST030375_R2 is misspelled in labshet; labsheet has two sample ("SJST030375_D1", "SJST032952_D1") which count_matrix doesn't. Total: Those 3 samples must fall out of our analysis; we must duplicate SJST031362_D1 in the labsheet

# editing misspelling 
labsheet$smpID[10] <- 'SJST030375_R2'
# deleting samples there're only in labsheet
labsheet <- labsheet[!labsheet$smpID %in% c("SJST030375_D1", "SJST032952_D1"), ]
# duplicate SJST033312_D1 in the labsheet
row1 <- labsheet[labsheet$smpID %in% 'SJST031362_D1', ]
row1$smpID <- 'SJST031362_D1.1'
labsheet <- rbind(labsheet, row1)



length(labsheet$smpID) #26 samples (1 uplicate) but still unique smpIDs
length(colnames(count_matrix))#26 samples (1 duplicate) but still unique smpIDs
length(intersect(labsheet$smpID, colnames(count_matrix))) #25
setdiff(labsheet$smpID, colnames(count_matrix))

# save(count_matrix, labsheet, file = '/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/analysis/reviewed_smpID_countmtx_labsheet.rda')
```

# Downstream analysis

### Spotting outlier (SJST034534_D1)

``` r
load("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/bulk_RNA/objects/Round1/reviewed_smpID_countmtx_labsheet.rda")
labsheet
count_matrix
labsheet$group<- gsub(' ', '_', labsheet$group)

# Convert gene ID
rownames(count_matrix) <- str_split_fixed(as.character(rownames(count_matrix)), "[.]", 2)[,1] # good it doesn't have duplicates ENSEMBL here
count_matrix <- data.frame(count_matrix)
ens <- rownames(count_matrix) #genes from DE
library(org.Hs.eg.db)
symbols <- mapIds(org.Hs.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = 'ENSEMBL')

symbols <- symbols[!is.na(symbols)]

symbols <- data.frame(ENSEMBL = names(symbols),
                      SYMBOL = as.vector(symbols))

count_matrix$ENSEMBL <- rownames(count_matrix)

count_matrix <- merge(count_matrix, symbols, by='ENSEMBL')
count_matrix <- count_matrix[!duplicated(count_matrix$SYMBOL), ]
rownames(count_matrix) <- count_matrix$SYMBOL
count_matrix$ENSEMBL <- NULL
count_matrix$SYMBOL <- NULL
head(count_matrix); dim(count_matrix) #33909    26


pca <- prcomp(t(count_matrix)) 
aux <- as.data.frame(pca$x[, 1:3]) 
scores <- merge(labsheet, aux, by.y=0, by.x="smpID", all.x=T)

png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/PCA_SJST034534_D1_outlier_RNAseq_counts.png", width = 8, height = 6, units = "in", res = 300)
# Plot it 
library(ggplot2); theme_set(theme_classic())
library(ggrepel)
ggplot(scores, aes(x=PC1, y=PC2, colour=factor(group), shape = smp_type)) +
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
    geom_text_repel(aes(label = smpID), size = 3, box.padding = 0.5, point.padding = 0.5) +
  ggtitle("Raw Counts - whole transcriptome") 

dev.off()
```

# OBS: the “batch correction through matrix design” affects only the Differential Expression Analysis. It’s only the DESeq2 normalization we can see on PCA/tSNE projections

``` r
# Load & Prepare Data  -----
library(DESeq2)
library(sva)
library(dplyr)
library(stringr)

load("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/bulk_RNA/objects/Round1/reviewed_smpID_countmtx_labsheet.rda")
labsheet
count_matrix
labsheet$group<- gsub(' ', '_', labsheet$group)

# Convert gene ID
rownames(count_matrix) <- str_split_fixed(as.character(rownames(count_matrix)), "[.]", 2)[,1] # good it doesn't have duplicates ENSEMBL here
count_matrix <- data.frame(count_matrix)
ens <- rownames(count_matrix) #genes from DE
library(org.Hs.eg.db)
symbols <- mapIds(org.Hs.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = 'ENSEMBL')

symbols <- symbols[!is.na(symbols)]

symbols <- data.frame(ENSEMBL = names(symbols),
                      SYMBOL = as.vector(symbols))

count_matrix$ENSEMBL <- rownames(count_matrix)

count_matrix <- merge(count_matrix, symbols, by='ENSEMBL')
count_matrix <- count_matrix[!duplicated(count_matrix$SYMBOL), ]
rownames(count_matrix) <- count_matrix$SYMBOL
count_matrix$ENSEMBL <- NULL
count_matrix$SYMBOL <- NULL
head(count_matrix); dim(count_matrix) #33909    26

# remove outlier
count_matrix <- count_matrix[, !colnames(count_matrix) %in% "SJST034534_D1"]
dim(count_matrix) #25 smp
labsheet <- labsheet[!labsheet$smpID %in% "SJST034534_D1", ]
dim(labsheet)#25 smp

# Create histologic_for_outlier column 
labsheet$histologic_for_outlier <- 'Others'
outlier_smp <- c('SJST033491_D1',
                 'SJST032952_D2',
                 'SJST032952_D4',
                 'SJST032767_D2',
                 'SJST033312_D1',
                 'SJST033835_D1',
                 'SJST031920_D1')
his_diag <- labsheet[labsheet$smpID %in% outlier_smp, ]$histologic_diagnosis

labsheet[labsheet$smpID %in% outlier_smp[1], ]$histologic_for_outlier <- his_diag[1]
labsheet[labsheet$smpID %in% outlier_smp[2], ]$histologic_for_outlier <- his_diag[2]
labsheet[labsheet$smpID %in% outlier_smp[3], ]$histologic_for_outlier <- his_diag[3]
labsheet[labsheet$smpID %in% outlier_smp[4], ]$histologic_for_outlier <- his_diag[4]
labsheet[labsheet$smpID %in% outlier_smp[5], ]$histologic_for_outlier <- his_diag[5]
labsheet[labsheet$smpID %in% outlier_smp[6], ]$histologic_for_outlier <- his_diag[6]
labsheet[labsheet$smpID %in% outlier_smp[7], ]$histologic_for_outlier <- his_diag[7]


# make it identical (same sequence)
count_matrix <- count_matrix[, labsheet$smpID]
identical(colnames(count_matrix), labsheet$smpID) #TRUE
# add rownames into labsheet
rownames(labsheet) <- labsheet$smpID


# Differential Expression -----
# make it identical (same sequence)
count_matrix <- count_matrix[, labsheet$smpID]
identical(colnames(count_matrix), labsheet$smpID) #TRUE
# add rownames into labsheet
rownames(labsheet) <- labsheet$smpID
# leave labsheet with only variables we want in the matrix design 
labsheet <- labsheet[, c("group", "smp_type") ]


head(labsheet)
# Create DESEq object 
dds <- DESeqDataSetFromMatrix(count_matrix, colData = labsheet, design = ~ smp_type + group)

# Important step for further normalization during DESEq()
dds <- estimateSizeFactors(dds)
# sizeFactors(dds) #Idk what it does ... neither if we need to run it
# Pre-filtering 
smallestGroupSize <- min(as.vector(table(dds$group)))
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
# indicate direction for comparison
# dds$condition <- relevel(dds$group, ref = "kinase-fused_tumor")
# dds$group <- relevel(dds$group, ref = "ETV6-NTRK3_fused_tumor")
# dds$group <- relevel(dds$group, ref = "kinase-fused_tumor")

# Run Differential Expression analysis
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", 
                                 "kinase-fused_tumor", #direction
                                 "ETV6-NTRK3_fused_tumor"))
res_kinaseDirection <- as.data.frame(res)

# Filter results for significant genes
# Subset stat. significant genes
res_pos <- subset(res_kinaseDirection, res_kinaseDirection$padj < 0.05 & res_kinaseDirection$log2FoldChange >= 0.5)
res_neg <- subset(res_kinaseDirection, res_kinaseDirection$padj < 0.05 & res_kinaseDirection$log2FoldChange <= -0.5)
res_sig <- rbind(res_pos, res_neg)
dim(res_sig) #617   6
# PLOT - FoldChange distribution - histogram
hist(res_sig$log2FoldChange, 
     main = "FoldChange dist. (s.sig DEGs)", 
     xlab = "FoldChange",
     ylab = "N of genes")



# PCA -----

# PCA - input: count matrix
png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/PCA_RNAseq_counts.png", width = 8, height = 6, units = "in", res = 300)

counts <- counts(dds, normalized=FALSE)
pca <- prcomp(t(counts)) 
aux <- as.data.frame(pca$x[, 1:3]) 
labsheet$smpID <- rownames(labsheet)
scores <- merge(labsheet, aux, by.y=0, by.x="smpID", all.x=T)
# Plot it 
library(ggplot2); theme_set(theme_classic())
ggplot(scores, aes(x=PC1, y=PC2, colour=factor(group), shape = smp_type)) +
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
  ggtitle("Raw Counts - whole transcriptome") 

dev.off()
```

``` r
# t-SNE 
png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/tsne_RNAseq_counts.png", width = 8, height = 6, units = "in", res = 300)

library(Rtsne)
set.seed(234234)
tsne_realData <- Rtsne(t(counts), perplexity=2, check_duplicates = FALSE) # #function to run t-sne
pdata.teste.tsne <- scores #renaming metadata
pdata.teste.tsne$tSNE1 <- tsne_realData$Y[,1]
pdata.teste.tsne$tSNE2 <- tsne_realData$Y[,2]

library(ggplot2); theme_set(theme_classic())
ggplot(pdata.teste.tsne, aes(x=tSNE1, y=tSNE2, colour=factor(group), shape = smp_type)) +
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
  ggtitle("Raw Counts - whole transcriptome") 

dev.off()
```

``` r
# PCA - input: normalized_counts matrix
png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/PCA_RNAseq_normalized_batchcorrected.png", width = 8, height = 6, units = "in", res = 300)

normalized_counts <- counts(dds, normalized=TRUE)
pca <- prcomp(t(normalized_counts)) 
aux <- as.data.frame(pca$x[, 1:3]) 
labsheet$smpID <- rownames(labsheet)
scores <- merge(labsheet, aux, by.y=0, by.x="smpID", all.x=T)
# Plot it 
library(ggplot2); theme_set(theme_classic())
ggplot(scores, aes(x=PC1, y=PC2, colour=factor(group), shape = smp_type)) +
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
  ggtitle("Normalized / batch corrected (mtx design) - whole transcriptome") 

dev.off()
```

``` r
# t-SNE 
png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/tsne_RNAseq_normalized_batchcorrected.png", width = 8, height = 6, units = "in", res = 300)

library(Rtsne)
set.seed(234234)
tsne_realData <- Rtsne(t(normalized_counts), perplexity=2, check_duplicates = FALSE) # #function to run t-sne
pdata.teste.tsne <- scores #renaming metadata
pdata.teste.tsne$tSNE1 <- tsne_realData$Y[,1]
pdata.teste.tsne$tSNE2 <- tsne_realData$Y[,2]

library(ggplot2); theme_set(theme_classic())
ggplot(pdata.teste.tsne, aes(x=tSNE1, y=tSNE2, colour=factor(group), shape = smp_type)) +
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
  ggtitle("Normalized / batch corrected (mtx design) - whole transcriptome") 

dev.off()
```

``` r
# Volcano plot  -----
# All DEGs
png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/Volcano_RNAseq_AllDEGs.png", width = 8, height = 6, units = "in", res = 300)

library(EnhancedVolcano)
#rownames(res_sig) <- res_sig$SYMBOL # don't run it if it doesn't need to
EnhancedVolcano(#res_sig,
                res_kinaseDirection, #DEG output not filtered yet
                lab = "", 
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = c("FLG"  , "SLC6A15"    ,  "ALK"   ,  "CTXND1"  ,    "PAX3"),
                pCutoff = 10e-2, #pvalue cutoff line
                FCcutoff = 0.5, #foldChange cutoff line
                xlab = bquote(~Log[2]~ 'fold change'),
                pointSize = 4.0,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                #boxedLabels = TRUE,
                #colAlpha = 4/5,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                #drawConnectors = TRUE,
                #widthConnectors = 1.0,
                #colConnectors = 'black',
                title = "kinase-fused_tumor direction") + theme_classic()

dev.off()
```

``` r
# PCA
# PCA - input: normalized_counts_DEGs matrix
png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/PCA_RNAseq_normalized_batchcorrected_DEGs.png", width = 8, height = 6, units = "in", res = 300)

normalized_counts_DEGs <- normalized_counts[rownames(res_sig), ]
pca <- prcomp(t(normalized_counts_DEGs))
aux <- as.data.frame(pca$x[, 1:3])
labsheet$smpID <- rownames(labsheet)
scores <- merge(labsheet, aux, by.y=0, by.x="smpID", all.x=T)
# Plot it
library(ggplot2); theme_set(theme_classic())
ggplot(scores, aes(x=PC1, y=PC2, colour=factor(group), shape = smp_type)) +
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
  ggtitle("Normalized / batch corrected (mtx design) - only DEGs (617)")

dev.off()
```

``` r
# Biological Pathways Analysis -----
# IMPROVIMENTS TO BE DONE: i) add Wnt signaling of Down regulated genes into plots (to ontology and betwork plots); ii) add foldChange into plots 
library(clusterProfiler)
library(org.Hs.eg.db)
DEGs <- res_sig
DEGs$gene <- rownames(DEGs)
up_genes_kinase <- DEGs[DEGs$log2FoldChange > 0 & DEGs$padj < 0.05, ]$gene
down_genes_kinase <- DEGs[DEGs$log2FoldChange < 0 & DEGs$padj < 0.05, ]$gene

# Up regulated genes - Gene Ontology
ego2_up <- enrichGO(gene = up_genes_kinase,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 # ont  One of "BP", "MF", and "CC" subontologies,                                   # or "ALL" for all three
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

ego2_up@result[order(ego2_up@result$Count, decreasing = TRUE),][1:10,]
barplot(ego2_up, drop=TRUE, main = "")
#clusterProfiler::dotplot(ego2_up, showCategory=30) + ggtitle("")
ego2_up_plot <- ego2_up@result[ego2_up@result$Description %in% 
                                     c("cilium organization", "cilium assembly", "axoneme assembly", 
           "cilium movement", "microtubule-based movement", 
           "microtubule bundle formation", "axonemal dynein complex assembly",
           "cilium or flagellum-dependent cell motility"), ]

png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/enrichmentBarplot_upDEGs_RNAseq.png", width = 8, height = 6, units = "in", res = 300)

library(ggplot2); theme_set(theme_classic())
ggplot(ego2_up_plot, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
  geom_col() +
  labs(x = "Gene Count", y = "", title = "Tumor-Kinase Up Regulated Pathways") +
  scale_fill_gradient(low = "red", high = "blue", name = "p.adjust") 

dev.off()
```

``` r
# Up regulated genes - Gene-Concept Network
# Trying to catch redundant pathways 
table(unlist(strsplit(ego2_up_plot$geneID, "/")))
gene_list <- strsplit(ego2_up_plot$geneID, "/")
names(gene_list) <- ego2_up_plot$Description
UpSetR::upset(fromList(gene_list), 
              order.by = "freq", 
              nsets = 8)
# Remove c("cilium assembly", "axome assembly", "microtubele bundle formation", "axonemal dynein complex assembly", "cillium moviment")
# ego2_up_plot[!ego2_up_plot$Description %in% c("cilium assembly",), ]
ego2_up_sub <- filter(ego2_up, Description %in% ego2_up_plot$Description)
edox <- setReadable(ego2_up_sub, 'org.Hs.eg.db', 'ENTREZID')
geneList <- up_genes_kinase

png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/enrichmentCnetplot_upDEGs_RNAseq.png", width = 8, height = 6, units = "in", res = 300)

p1 <- cnetplot(edox, foldChange=geneList, 
               #circular = TRUE, colorEdge = TRUE, 
               #showCategory = "", 
               ) + ggtitle("Tumor-Kinase Up Regulated Pathways ") 


p1
dev.off()

# NOT INFORMATIVE PLOT vvvvvvvvv
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(enrichplot)
# library(GOSemSim)
# library(DOSE)
# ego2_par <- pairwise_termsim(ego2_up)
# p2 <- emapplot(ego2_par, cex_category=1,layout="kk") 
# p2
# NOT INFORMATIVE PLOT ^^^^^^^^^^
```

``` r
# Down regulated genes - Gene Ontology
ego2_down <- enrichGO(gene = down_genes_kinase,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 # ont  One of "BP", "MF", and "CC" subontologies,                                   # or "ALL" for all three
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

ego2_down@result[order(ego2_down@result$Count, decreasing = TRUE),][1:10,]
barplot(ego2_down, drop=TRUE, main = "")
# clusterProfiler::dotplot(ego2_down, showCategory=30) + ggtitle("")
ego2_down_plot <- ego2_down@result[ego2_down@result$Description %in% c("Wnt signaling pathway",
                                     "epithelial cell proliferation",
                                     "skeletal system morphogenesis",
                                     "embryonic skeletal system development",
                                     "bone development"), ]

png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/enrichmentBarplot_downDEGs_RNAseq.png", width = 8, height = 6, units = "in", res = 300)

library(ggplot2); theme_set(theme_classic())
ggplot(ego2_down_plot, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
  geom_col() +
  labs(x = "Gene Count", y = "", title = "Tumor-Kinase Down Regulated Pathways") +
  scale_fill_gradient(low = "red", high = "blue", name = "p.adjust") 

dev.off()
```

``` r
# Down regulated genes - Gene-Concept Network
ego2_down_sub <- filter(ego2_down, Description %in% ego2_down_plot$Description) #trying to add wnt pathways into network plot
ego2_down_sub@result$p.adjust <- 0.001 # this is only to make it go to the plot !
ego2_down_sub@result$qvalue <- 0.001 # this is only to make it go to the plot !
edox <- setReadable(ego2_down_sub, 'org.Hs.eg.db', 'ENTREZID')
geneList <- down_genes_kinase

png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/enrichmentCnetplot_downDEGs_RNAseq.png", width = 8, height = 6, units = "in", res = 300)

p1 <- cnetplot(edox, foldChange=geneList, 
               #circular = TRUE, colorEdge = TRUE, 
               #showCategory = "", 
               ) 
p1
dev.off()

# NOT INFORMATIVE PLOT vvvvvvvvv
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(enrichplot)
# library(GOSemSim)
# library(DOSE)
# ego2_par <- pairwise_termsim(ego2_down)
# p2 <- emapplot(ego2_par, cex_category=1,layout="kk") 
# p2
# NOT INFORMATIVE PLOT ^^^^^^^^
```

``` r
# Heatmap -----
### Load expression matrix 
library(dplyr)
library(DESeq2)
library(stringr)

# Load expression matrix 
load("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/bulk_RNA/objects/Round1/reviewed_smpID_countmtx_labsheet.rda")
labsheet
labsheet$group<- gsub(' ', '_', labsheet$group)

# Convert gene ID
rownames(count_matrix) <- str_split_fixed(as.character(rownames(count_matrix)), "[.]", 2)[,1] # good it doesn't have duplicates ENSEMBL here
count_matrix <- data.frame(count_matrix)
ens <- rownames(count_matrix) #genes from DE
library(org.Hs.eg.db)
symbols <- mapIds(org.Hs.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = 'ENSEMBL')

symbols <- symbols[!is.na(symbols)]

symbols <- data.frame(ENSEMBL = names(symbols),
                      SYMBOL = as.vector(symbols))

count_matrix$ENSEMBL <- rownames(count_matrix)

count_matrix <- merge(count_matrix, symbols, by='ENSEMBL')
count_matrix <- count_matrix[!duplicated(count_matrix$SYMBOL), ]
rownames(count_matrix) <- count_matrix$SYMBOL
count_matrix$ENSEMBL <- NULL
count_matrix$SYMBOL <- NULL
head(count_matrix); dim(count_matrix) #33909    26

# remove outlier
count_matrix <- count_matrix[, !colnames(count_matrix) %in% "SJST034534_D1"]
dim(count_matrix) #25 smp
labsheet <- labsheet[!labsheet$smpID %in% "SJST034534_D1", ]
dim(labsheet)#25 smp
# make it identical (same sequence)
count_matrix <- count_matrix[, labsheet$smpID]
identical(colnames(count_matrix), labsheet$smpID) #TRUE
# add rownames into labsheet
rownames(labsheet) <- labsheet$smpID

# Add fusion variable from the labsheet (excel file) into our metadata (the "labsheet" object)
# from the current labsheet shared via One Drive:
smpID <- c("SJST033983_D1", "SJST032767_D1", "SJST032767_D2", "SJST034397_D1", "SJST031920_D1",
           "SJST033491_D1", "SJST033389_D1", "SJST031690_D1", "SJST031362_D1", "STST030375_R2",
           "SJST030375_D1", "SJST031792_D1", "SJST033791_D1", "SJST033835_D1", "SJST033308_D1",
           "SJST034534_D1", "SJST034036_D1", NA, NA, "SJST034815_D1", "SJST032952_D2",
           "SJST032952_D1", "SJST032952_D4", "SJST033312_D1", "SJST032838_D1", "SJBT032721_D1",
           "SJIFS031085_D1", "SJST030567_D1", "SJST030433_D1", NA)

fusion_type <- c("NRF1::BRAF", "PLEKHH2::ALK", "PLEKHH2::ALK", "TPR::NTRK1", "TPM3::NTRK1",
                 "TPM3::NTRK1", "TPR::NTRK1", "LMNA::NTRK1", "LMNA::NTRK1", "PDE4DIP::NTRK1",
                 "PDE4DIP::NTRK1", "EML4::NTRK3", "RBPMS::NTRK3", "RCC1::ALK", 
                 "EML4::NTRK3","TTYH3::BRAF","TPM3::NTRK1","EML4::ALK Outside lab",
                 'TPM3::NTRK1 Outside lab',"TPM3::NTRK1","ETV6:NTRK3","ETV6:NTRK3","ETV6:NTRK3",
                 'ETV6::NTRK3','ETV6::NTRK3','ETV6::NTRK3','ETV6::NTRK3','ETV6::NTRK3',
                 'ETV6::NTRK3','ETV6::NTRK3')

# Create a dataframe
fusion_df <- data.frame(smpID = smpID, fusion_type = fusion_type)
fusion_df <- fusion_df[!fusion_df$smpID %in% NA, ]
fusion_df$smpID[fusion_df$smpID %in% 'STST030375_R2'] <- 'SJST030375_R2'

table(duplicated(fusion_df$smpID))
setdiff(rownames(labsheet), fusion_df$smpID) 
length(intersect(rownames(labsheet), fusion_df$smpID)) #24 #SJST031362_D1.1 we have its annotation for SJST031362_D1. So it's okay move on with these 24 sample ID 

# Merge fusion_type variable into labsheet
labsheet_hm <- merge(labsheet, fusion_df, by = 'smpID')


# Load normalized expression matrix 
dds <- readRDS("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/objects/DESeq2_output.rds")
normalized_counts <- counts(dds, normalized=TRUE)
# Filter normalized_counts by labsheet
normalized_counts_hm <- normalized_counts[, colnames(normalized_counts) %in% labsheet_hm$smpID]

### PLOT - heatmap - whole RNAseq - focus on NTRK label
library(pheatmap)
library(RColorBrewer)
normalized_counts_hm # exp matrix
labsheet_hm # metadata

my_sample_col <- data.frame(row.names = labsheet_hm$smpID, 
                            Group = labsheet_hm$group, 
                            Histology = labsheet_hm$histologic_diagnosis,
                            Fusion = labsheet_hm$fusion_type) # patient ID in row names; any columns are for sample information
my_sample_col[my_sample_col$Fusion %in% "ETV6:NTRK3", ]$Fusion <- "ETV6::NTRK3" #correcting label typo

NTRK <- my_sample_col[grep('NTRK', my_sample_col$Fusion), ]
NTRK_butETV6 <- rownames(NTRK[!NTRK$Fusion %in% c('ETV6::NTRK3'), ])
NTRK_onlyETV6 <- rownames(NTRK[!rownames(NTRK) %in% NTRK_butETV6, ])

my_sample_col$Fusion_concat <- NA
my_sample_col[rownames(my_sample_col) %in% NTRK_butETV6, ]$Fusion_concat <- 'Others_NTRK'
my_sample_col[rownames(my_sample_col) %in% NTRK_onlyETV6, ]$Fusion_concat <- 'ETV6_NTRK3'
my_sample_col[my_sample_col$Fusion_concat %in% NA, ]$Fusion_concat <- 'Not_NTRK'


ann_colors = list(
  Fusion_concat = c(ETV6_NTRK3 = 'black', Not_NTRK = 'blue', Others_NTRK = 'red'))

my_sample_col$Histology <- NULL
# pheatmap(normalized_counts_hm, 
#          #annotation_row = my_probe_col, 
#          annotation_col = my_sample_col[, ], # "group" [1] label is redundant here. Leave it out
#          show_rownames = FALSE,
#          scale = "row",
#          annotation_colors = ann_colors,
#          main = 'Whole RNAseq - normalized (18,777 features)')


# Load exp matrix 
DEGs <- readRDS("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/objects/DEGs_Volcanoplot_input_Kinase_direction.rds")
head(DEGs)
dim(DEGs)

# heatmap based only on DEGs
png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/heatmap_allDEGs_RNAseq.png", width = 8, height = 6, units = "in", res = 300)

pheatmap(normalized_counts_hm[rownames(normalized_counts_hm) %in% rownames(DEGs), ], 
         #annotation_row = my_probe_col, 
         annotation_col = my_sample_col[, ],
         show_rownames = FALSE,
         scale = "row",
         annotation_colors = ann_colors,
         main = 'DEGs - normalized (634 features)')

dev.off()
```

``` r
# Gene Exp. Correlation -----
library(ggplot2)
library(corrplot)
length(colnames(normalized_counts_hm))
length(rownames(my_sample_col))
intersect(colnames(normalized_counts_hm),
          rownames(my_sample_col))
identical(colnames(normalized_counts_hm),
          rownames(my_sample_col)) #FALSE

my_sample_col <- my_sample_col[colnames(normalized_counts_hm), ]
identical(colnames(normalized_counts_hm),
          rownames(my_sample_col)) #TRUE

colnames(normalized_counts_hm) <- my_sample_col$Fusion_concat

cor_matrix <- cor(normalized_counts_hm)

df_color <- data.frame(smp_name = colnames(normalized_counts_hm),
                       color = NA)
# # df_color[df_color$smp_name %in% "ETV6_NTRK3", ]$color <- "black"
# # df_color[df_color$smp_name %in% "Not_NTRK", ]$color <- "blue"
# # df_color[df_color$smp_name %in% "Others_NTRK", ]$color <- "red"
# df_color$color <- c("black", "black", "blue", "red", "black", "blue", "red",
#                     "black", "black", rep("red", 5), "black", "red", "red", "black", "blue", "black", "blue", rep("red", 3)) # accordingly to the correlation plot sample display 
# smp_colors <- df_color$color
# 
# 
# cor_test_mat <- cor.mtest(normalized_counts_hm, conf.level = 0.95)
# cor_test_mat <- cor_test_mat$p
# 
# corrplot(cor_matrix, method = "color", type = "full",
#          tl.col = smp_colors, tl.srt = 45, order = "hclust",
#          addrect = 3, rect.col = 'blue', rect.lwd = 3,
#          p.mat = cor_test_mat, sig.level = 0.001)



df_color$color <- c("black", "black", "red", rep("black", 4), "blue",
                    "black", "black", "blue", "blue", "red", "blue",
                    rep("red", 10)) # accordingly to the correlation plot sample display 
smp_colors <- df_color$color


cor_test_mat <- cor.mtest(normalized_counts_hm[rownames(normalized_counts_hm) %in% rownames(DEGs), ], conf.level = 0.95)
cor_test_mat <- cor_test_mat$p

cor_matrix_degs <- cor(normalized_counts_hm[rownames(normalized_counts_hm) %in% rownames(DEGs), ]
)

png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/correlationPlot_DEGs_RNAseq.png", width = 8, height = 6, units = "in", res = 300)

corrplot(cor_matrix_degs, method = "color", type = "full",
         tl.col = smp_colors, tl.srt = 45, order = "hclust",
         addrect = 3, rect.col = 'blue', rect.lwd = 3,
         p.mat = cor_test_mat, sig.level = 0.001)

dev.off()
```

``` r
dim(normalized_counts) #18777    25
head(normalized_counts)
dim(labsheet) #25  2
head(labsheet)
labsheet$smp_ID <- rownames(labsheet)
# Stemness (no stat. sig difference) -----
load("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/objects/RNA_stemness_model.Rda") # mm = modelo
w = mm$w
w[1:5]
length(w) #12953 (number of genes on the model)
# Filter matrix expression by the genes on stemenss model
matrix = normalized_counts
length(intersect(rownames(matrix), names(w))) #12206 
predict.DATA = matrix[rownames(matrix) %in% names(w) ,]
length(rownames(predict.DATA)) # 12206
w = w[ rownames(predict.DATA) ]
length(intersect(names(w),rownames(predict.DATA))) # 12206 
w[1:5]
length(names(w)) # 12206
is.vector(w) #TRUE

# Score the Matrix `X` using Spearman correlation
ss = apply( predict.DATA, 2, function(z) {cor( z, w, method="sp", use="complete.obs" )} )
ss[1:5]


# Scale the scores into a ratio from 0 to 1 and store as data frame.
ss <- ss - min(ss)
ss <- ss / max(ss)
ss <- as.data.frame(ss)
colnames(ss) <- "Stemness_RNAseq" 
ss$smp_ID <- rownames(ss)

# Plot PCA/tSNE + Stemness
labsheet <- merge(labsheet, ss, by="smp_ID")
normalized_counts <- na.omit(normalized_counts)
pca <- prcomp(t(normalized_counts)) 
aux <- as.data.frame(pca$x[, 1:3]) 
scores <- merge(labsheet, aux, by.y=0, by.x="smp_ID", all.x=T)
png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/PCA_normalized_stemness_st_RNAseq.png", width = 8, height = 8, units = "in", res = 300)

# PCA 
library(ggplot2); theme_set(theme_classic())
ggplot(scores, aes(PC1, PC2)) + geom_point(size=5, aes( fill=Stemness_RNAseq, color = group), pch = 21, stroke = 1.5) + 
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  scale_color_manual(values = c('slateblue3', 'wheat3')) +
  ylab(paste0('PC2 ', summary(pca)$importance[2, 2] * 100, '%')) +
  xlab(paste0('PC1 ', summary(pca)$importance[2, 1] * 100, '%') ) + theme_bw() +
  ggtitle('Normalized counts - Stemness prediction')

dev.off()


library(Rtsne)
normalized_counts <- na.omit(normalized_counts)
set.seed(234234)
tsne_realData <- Rtsne(t(normalized_counts), perplexity=3, check_duplicates = FALSE) # #function to run t-sne
pdata.teste.tsne <- scores #renaming metadata
pdata.teste.tsne$tSNE1 <- tsne_realData$Y[,1]
pdata.teste.tsne$tSNE2 <- tsne_realData$Y[,2]
png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/tsne_normalized_stemness_st_only_RNAseq.png", width = 8, height = 6, units = "in", res = 300)
# tSNE
library(ggplot2); theme_set(theme_classic())
ggplot(pdata.teste.tsne, aes(tSNE1, tSNE2)) + geom_point(size=5, aes( fill=Stemness_RNAseq, color = group), pch = 21, stroke = 1.5) + 
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  scale_color_manual(values = c('slateblue3', 'wheat3')) +
  # ylab(paste0('PC2 ', summary(pca)$importance[2, 2] * 100, '%')) +
  # xlab(paste0('PC1 ', summary(pca)$importance[2, 1] * 100, '%') ) + theme_bw() +
  ggtitle('Normalized counts - Stemness prediction')

dev.off()


png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_3/scripts/Documentation/images/boxplot_stemness_st_RNAseq.png", width = 8, height = 8, units = "in", res = 300)
library(ggplot2); theme_set(theme_classic())
ggplot(scores[, ], aes(x=group, y=Stemness_RNAseq)) + 
  geom_boxplot(fill= c('slateblue3', 'wheat3'), 
    outlier.color = NA) + 
  geom_jitter (alpha=0.2)  +
  xlab("Groups") + 
  ylab("Stemness") + 
  theme(axis.title.x = element_text(colour="black", size = 12, vjust=-1), 
        axis.title.y = element_text(colour="black", size = 12), 
        axis.text.x = element_text(angle = 45, vjust= 1, size = 12, hjust = 1)) +
  #facet_wrap(~ sample_char) +
  ggtitle("Stemness diff. beetween groups (p = 0.7634)") 

dev.off()

t.test(scores[scores$group %in% "ETV6-NTRK3_fused_tumor", ]$Stemness_RNAseq,
       scores[scores$group %in% "kinase-fused_tumor", ]$Stemness_RNAseq) # p-value = 0.7634

library(rstatix)
scores %>%
wilcox_test(Stemness_RNAseq ~ group, p.adjust.method = "none") %>%
  add_significance() #p-value = 0.315
```
