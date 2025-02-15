library(dplyr)
library(minfi)
library(DMRcate)
library(Gviz)
library(RColorBrewer)
library(GenomicRanges)
library(rtracklayer)
library(HelpersMG)
library(data.table)


### Loading Differential Methylated Regions (DMRs) and Differential Expression Genes
dmr.table <- readRDS("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/DNAmet/Round_2/reviewed_DM_direction/dmr_table.rds")
DEGs <- readRDS("/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/bulk_RNA/objects/Round2/reviewed_DEGs_direction/Volcano_Kinase_direction_DEGs_plot.rds")
DEGs$gene <- rownames(DEGs)




# Find intersecting genes within DMR-gene annotations and DEGs
# Get DMR-genes list  (making it a flat gene list)
dmr_genes <- dmr.table$overlapping.genes
list_chars <- list()
for(i in 1:length(dmr_genes)) {
  chars <- strsplit(dmr_genes, split = ",")[[i]]
  chars <- gsub(' ','', chars)
  list_chars[[i]] <- chars
  
}

dmr_genes <- as.vector(do.call(rbind, list_chars))
dmr_genes <- dmr_genes[!dmr_genes %in% NA]

# Get DEGs genes lits (making it a flat gene list)
deg_genes <- DEGs$gene
list_chars <- list()
for(i in 1:length(deg_genes)) {
  chars <- strsplit(deg_genes, split = ",")[[i]]
  chars <- gsub(' ','', chars)
  list_chars[[i]] <- chars
  
}

deg_genes <- as.vector(do.call(rbind, list_chars))
deg_genes <- deg_genes[!deg_genes %in% NA]

length(dmr_genes) #1556 genes (DMR-CpGs annoated to genes)
length(deg_genes) #617 genes (the DEGs them selves)

length(intersect(dmr_genes,
                 deg_genes)) # 22 common genes  

inters_genes <- intersect(dmr_genes,
                          deg_genes)



# Add regulation directions to the intersected genes
##### eg: gene A is in a hypo methylated region AND it's expression is up regulated

# inters_genes in DEGs - adding directions (up and down regulated)
DEGs_inter <- DEGs[DEGs$gene %in% inters_genes, ]
DEGs_inter$regulation <- NA
DEGs_inter[DEGs_inter$log2FoldChange > 0, ]$regulation <- "up"
DEGs_inter[DEGs_inter$log2FoldChange < 0, ]$regulation <- "down"

# inters_genes in DMRs - adding directions (up and down regulated)
subset_dmr_tables <- list()
for (i in seq_along(inters_genes)) {
  subset_dmr_tables[[i]] <- dmr.table[grepl(inters_genes[i], dmr.table$overlapping.genes), ]
}

dmr.table_inter <- do.call(rbind, subset_dmr_tables)
dmr.table_inter$regulation <- NA
dmr.table_inter[dmr.table_inter$meandiff > 0, ]$regulation <- "up"
dmr.table_inter[dmr.table_inter$meandiff < 0, ]$regulation <- "down"

# up-down DNAmet-RNAseq relation 
# Up methylated CpGs
up_1 <- dmr.table_inter[dmr.table_inter$regulation %in% "up", ]$overlapping.genes
# Down regulated DEGs
down_1 <- DEGs_inter[DEGs_inter$regulation %in% "down", ]$gene #right way
# down_1 <- DEGs_inter[DEGs_inter$regulation %in% "up", ]$gene #testing it - this way we will get the same DMRs-DEGs we got at first. It happened because I also updated the DEGs direction 

# up_1 - adjusting elements to 1 gene: 1 element
list_chars <- list()
for(i in 1:length(up_1)) {
  chars <- strsplit(up_1, split = ",")[[i]]
  chars <- gsub(' ','', chars)
  list_chars[[i]] <- chars
  
}

up_1 <- as.vector(do.call(rbind, list_chars))
up_1 <- up_1[!up_1 %in% NA]

# down_1 - adjusting elements to 1 gene: 1 element
list_chars <- list()
for(i in 1:length(down_1)) {
  chars <- strsplit(down_1, split = ",")[[i]]
  chars <- gsub(' ','', chars)
  list_chars[[i]] <- chars
  
}

down_1 <- as.vector(do.call(rbind, list_chars))
down_1 <- down_1[!down_1 %in% NA]

length(intersect(up_1, down_1)) # 5
intersect(up_1, down_1) # "KSR1"     "WNK1"     "SMURF2"   "ARHGEF10" "PROM1"



# down-up DNAmet-RNAseq relation
# Down methylated CpGs
down_2 <- dmr.table_inter[dmr.table_inter$regulation %in% "down", ]$overlapping.genes
# Up regulated DEGs
up_2 <- DEGs_inter[DEGs_inter$regulation %in% "up", ]$gene #right way
# up_2 <- DEGs_inter[DEGs_inter$regulation %in% "down", ]$gene #testing it - this way we will get the same DMRs-DEGs we got at first. It happened because I also updated the DEGs direction 

# up_2
list_chars <- list()
for(i in 1:length(up_2)) {
  chars <- strsplit(up_2, split = ",")[[i]]
  chars <- gsub(' ','', chars)
  list_chars[[i]] <- chars
  
}

up_2 <- as.vector(do.call(rbind, list_chars))
up_2 <- up_2[!up_2 %in% NA]

# down_2
list_chars <- list()
for(i in 1:length(down_2)) {
  chars <- strsplit(down_2, split = ",")[[i]]
  chars <- gsub(' ','', chars)
  list_chars[[i]] <- chars
  
}

down_2 <- as.vector(do.call(rbind, list_chars))
down_2 <- down_2[!down_2 %in% NA]

length(intersect(up_2, down_2)) # 4
intersect(up_2, down_2)  # "NOS1"   "GPT2"   "EFCAB6" "CFAP54" - after reviewing DMRs 

# Subset dmr_tables by interseted genes + inverse DNAmet-GeneExp correlation
# Eg: dmr region 1 that contains high expressed gene A (from RNA seq) and is hypomethylated. That way it would be reasonable to take that gene A could be under DNAmet regulation

subset_dmr_tables <- list()
inters_genes <- c(intersect(up_1, down_1), intersect(up_2, down_2))
inters_genes <- unique(inters_genes)
for (i in seq_along(inters_genes)) {
  subset_dmr_tables[[i]] <- dmr.table[grepl(inters_genes[i], dmr.table$overlapping.genes), ]
}

dmr.table_inter <- do.call(rbind, subset_dmr_tables)
dmr.table_inter$meth_status <- NA
dmr.table_inter[dmr.table_inter$meandiff > 0, ]$meth_status <- "hyper_methylated"
dmr.table_inter[dmr.table_inter$meandiff < 0, ]$meth_status <- "down_methylated"

# From dmr.table_inter$overlapping.genes, get the genes also present in DEG from RNAseq
DEGs[DEGs$gene%in% c('KSR1',
                     'WNK1',
                     'SMURF2',
                     'AC019257.8',
                     'ARHGEF10',
                     'FGFBP2',
                     'PROM1',
                     'NOS1',
                     'GPT2',
                     'EFCAB6-AS1',
                     'EFCAB6',
                     'CFAP54'), ]$gene # "NOS1"     "GPT2"     "EFCAB6"   "CFAP54"   "PROM1"    "WNK1"  "ARHGEF10" "SMURF2"   "KSR1"  

# Correct multiple genes names from  overlapping.genes
dmr.table_inter$gene <- NA
dmr.table_inter$gene <- dmr.table_inter$overlapping.genes
dmr.table_inter$gene[4] <- "ARHGEF10"
dmr.table_inter$gene[5] <- "ARHGEF10"
dmr.table_inter$gene[6] <- "PROM1"
dmr.table_inter$gene[9] <- "EFCAB6"

# add the FC values manually because there's one gene duplicated cause it has two DMRs methylated. So we can't merge it directly 
DEGs_in_DMR <- DEGs[DEGs$gene%in% c('KSR1',
                                    'WNK1',
                                    'SMURF2',
                                    'ARHGEF10',
                                    'PROM1',
                                    'NOS1',
                                    'GPT2',
                                    'EFCAB6',
                                    'CFAP54'), ]



FC_info <- DEGs_in_DMR[, c("log2FoldChange", "gene")]
dmr.table_inter <- merge(dmr.table_inter, FC_info, by = 'gene')

# dmr.table_inter$log2FoldChange[1] <- -0.7067505
# dmr.table_inter$log2FoldChange[2] <- -0.9193195
# dmr.table_inter$log2FoldChange[3] <- -0.8706579
# dmr.table_inter$log2FoldChange[4] <- -0.8719283 
# dmr.table_inter$log2FoldChange[5] <- -0.8719283 #same FC, yes. It's right
# dmr.table_inter$log2FoldChange[6] <- -4.4897163 
# dmr.table_inter$log2FoldChange[7] <- 3.8963369 
# dmr.table_inter$log2FoldChange[8] <- 1.4430567
# dmr.table_inter$log2FoldChange[9] <- 1.5324998
# dmr.table_inter$log2FoldChange[10] <- 2.5872476





### Visualization
my_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462")

# Plot DMR-DEGs integration scatter plot


library(ggplot2); theme_set(theme_classic())
p <- ggplot(dmr.table_inter, aes(x = no.cpgs, y = width, col = gene, size = log2FoldChange)) +
  geom_point() +
  labs(
    title = "N of CpG by genome interval (bp)",
    x = "Frequency of CpGs",
    y = "Width - genomic interval (base pair)"
  ) +
  theme_classic() +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ meth_status)


## Quciclky check up on DNAmet-GeneExp relationship [sanity check]
dmr.table_inter[dmr.table_inter$overlapping.genes %in% c('NOS1',
                                                         'GPT2',
                                                         'CFAP54'), ]$meandiff # meandiff negative

DEGs[DEGs$gene%in% c('NOS1',
                     'GPT2',
                     'CFAP54'), ]$log2FoldChange #log2FoldChange positive

# So the directions are making sense. 




## Setting a minimum cutoff for $meandiff ~ 0.2 (by eye)

dmr.table_inter_filtered <- dmr.table_inter[dmr.table_inter$meandiff >= 0.16 |
                                              dmr.table_inter$meandiff <= -0.16, ]

png(filename = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_4/Plots/RNA_DNAmet_related_genes_summary_plot.png",  width = 8, height = 6, units = "in", res = 300)

library(ggplot2); theme_set(theme_classic())
p <- ggplot(dmr.table_inter_filtered, aes(x = no.cpgs, y = width, col = gene, size = log2FoldChange)) +
  geom_point() +
  labs(
    title = "N of CpG by genome interval (bp)",
    x = "Frequency of CpGs",
    y = "Width - genomic interval (base pair)"
  ) +
  theme_classic() +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ meth_status)

print(p)
dev.off()


## A summary table for the RNA-DNAmet integration 
dmr.table_inter_filtered[, c("gene",
                    "meandiff",
                    "maxdiff",
                    "meth_status",
                    "log2FoldChange")]

# Export table with the most reasonable targets to have their exp. regulated by DNA methylation
summ_table <- dmr.table_inter_filtered[, c("gene",
                             "meandiff",
                             "maxdiff",
                             "meth_status",
                             "log2FoldChange")]

summ_table$expression_status <- NA
summ_table[summ_table$log2FoldChange < 0, ]$expression_status <- "down_expressed"
summ_table[summ_table$log2FoldChange > 0, ]$expression_status <- "up_expressed"

summ_table$Group_direction <- "Kinase"

names(summ_table)[names(summ_table) == "meandiff"] <- "meandiff_DNAmet"
names(summ_table)[names(summ_table) == "maxdiff"] <- "maxdiff_DNAmet"


# # Save .csv file with target genes to send to Larissa
# write.csv(summ_table, file = "/mnt/scratch1/maycon/Larissa_inffibrosarcoma/scripts_git/round_4/Objects/RNA_DNAmet_related_genes.csv")

# END.

