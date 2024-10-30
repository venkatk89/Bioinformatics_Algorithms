library(DESeq2)
library(airway) # data is from this package
library(tidyverse)


# Load data

data("airway")

# get info about the samples from the package
sample_info <- as.data.frame(colData(airway))

# subset only the celline dexamethosone_treatment_status for this tutorial
sample_info <- sample_info[,c(2,3)]

# simple data cleaning
sample_info$dex <- gsub('trt', 'treated', sample_info$dex)
sample_info$dex <- gsub('untrt', 'untreated', sample_info$dex)

# assign column names. The row names should be the sample IDs. These sample_ids should be colum names of Counts data
names(sample_info) <- c('cellLine', 'dexamethasone')

# write the sample info to disk
# write.table(sample_info, file = "sample_info.csv", sep = ',', col.names = T, row.names = T, quote = F)

# get count of reads/fragments mapping to each gene in a matrix. rownames are gene_ids and columnames are sample_ids
countsData <- assay(airway)
head(countsData)

# write the counts info to disk
# write.table(countsData, file = "counts_data.csv", sep = ',', col.names = T, row.names = T, quote = F)


# DESeq2 requires colnames of countsdata to match rownames of sampledata. Check that
sum(rownames(sample_info) == colnames(countsData)) == length(rownames(sample_info))


# construct a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = countsData, colData = sample_info, design = ~dexamethasone)

dds


# pre-filtering: removing rows with low gene counts for reducing compute time
# keeping rows that have at least 10 reads total
dds <- dds[(rowSums(counts(dds)) >= 10),]


# set the factor level. If unset, factors are named alphabetically (making "treated" as reference)
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")


# as there are no technical replicates, they are not collapsed


# perform differential gene expression analysis
dds <- DESeq(dds)

# get the results of dge in a matrix
deseq_results <- results(dds)
deseq_results
summary(deseq_results)


# Get MA plot
plotMA(deseq_results_005)


# Get volcano plot

# Get -log(p_value) and group genes based on regulation type
deseq_results_tidy <- as.data.frame(deseq_results) %>%
  arrange(padj) %>%
  dplyr::mutate(logp = -log10(padj)) %>%
  dplyr::mutate(group = case_when(
    log2FoldChange >= 1 & padj <= 0.05 ~ "Up-regulated",
    log2FoldChange <= -1 & padj <= 0.05 ~ "Down-regulated",
    TRUE ~ "Not-significant"
  )) 


# Summarizing the count of DEGs in each category
table(deseq_results_tidy$group)

ggplot(deseq_results_tidy, aes(x = log2FoldChange, y = logp, colour = group, alpha = 0.5))+
  geom_point()+
  ggtitle("Differential gene expression", subtitle = "Volcano plot")+
  xlab("log2FoldChange")+
  ylab("-log10(Adjusted P-value)")+
  theme_classic()+
  theme(legend.position = "none")+
  geom_hline(yintercept=-log10(0.05),linetype="dashed")+
  geom_vline(xintercept= c(-1,1),linetype="dashed") 
