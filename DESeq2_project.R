# script to get data from airway package

#load libraries
library(DESeq2)
library(tidyverse)
library(airway)

#step 1: prepapring count data-----------

# read in counts data
counts_data <- read.csv('counts_data.csv')
head(counts_data)

# read in sample info
colData <- read.csv('sample_info.csv')

#making sure the row names in colData matches to the column names in count_data
all(colnames(counts_data) %in% rownames(colData))

#are they in the same order?
all(colnames(counts_data) == rownames(colData))

#step2: construct a DESeqDataSet object-----------
dds <- DESeqDataSetFromMatrix(countData  = counts_data,
                       colData = colData,
                       design = ~ dexamethasone)

dds

#pre-filtering: removing rows with low gene counts (this is an optional step)
#keeping rows that have atleast 10 rows total

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds

#set the factor level
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")

# NOTE- collapse technical replicates (and never the biological sets)

#step 3: Run DESwq -------------
dds <- DESeq(dds)
res <- results(dds)

res

#explore results ------------

summary(res)
res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

#contrasts
resultsNames(dds)

#e.g. : treated_4hrs, treated_8hrs, untreated

#results(dds, contrast = c("dexamethasone", "treated_4hrs", "untreated"))

# MA plot
plotMA(res)
