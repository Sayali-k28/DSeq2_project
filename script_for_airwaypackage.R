# Script to get data from airway package

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("airway")

library("airway")

data(airway)
airway

sample_info <- as.data.frame(colData(airway))
sample_info <- sample_info[,c(2,3)]

sample_info$dex <- gsub('trt', 'treated', sample_info$dex)
sample_info$dex <- gsub('untrt', 'untreated', sample_info$dex)
names(sample_info) <- c('cellLine', 'dexamethasone')
write.table(sample_info, file = "sample_info.csv", sep = ',', col.names = T, row.names = T, quote = F)

countsData <- assay(airway)
write.table(countsData, file = "counts_data.csv", sep = ',', col.names = T, row.names = T, quote = F)



# Step 1: Install BiocManager if you haven't already
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Step 2: Install DESeq2
BiocManager::install("DESeq2")

# Step 3: Verify the installation by loading the package
library(DESeq2)
