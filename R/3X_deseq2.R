library(DESeq2)       
library(dplyr)       
library(tibble) 

countData <- as.matrix(read.csv("genes.counts.csv", row.names = 1))

new_col <- c("3X01_POST", "3X01_PRE", "3X04_POST", "3X04_PRE", "3X05_POST", "3X05_PRE", "3X07_POST", "3X07_PRE",
             "3X09_POST", "3X09_PRE", "3X10_POST", "3X10_PRE", "3X11_POST", "3X11_PRE", "3X12_POST", "3X12_PRE",
             "3X14_POST", "3X14_PRE", "3X15_POST", "3X15_PRE", "3X17_POST", "3X17_PRE", "3X18_POST", "3X18_PRE", 
             "3X19_POST", "3X19_PRE", "3X20_POST", "3X20_PRE", "3X22_POST", "3X22_PRE", "3X23_POST", "3X23_PRE",
             "3X24_POST", "3X24_PRE", "3X25_POST", "3X25_PRE", "3X26_POST", "3X26_PRE", "3X28_POST", "3X28_PRE",
             "3X29_POST", "3X29_PRE", "3X30_POST", "3X30_PRE", "3X32_POST", "3X32_PRE", "3X35_POST", "3X35_PRE",
             "3X36_POST", "3X36_PRE", "3X37_POST", "3X37_PRE", "3X38_POST", "3X38_PRE", "3X39_POST", "3X39_PRE",
             "3X41_POST", "3X41_PRE", "3X42_POST", "3X42_PRE", "3X43_POST", "3X43_PRE", "3X44_POST", "3X44_PRE",
             "3X46_POST", "3X46_PRE", "3X47_POST", "3X47_PRE", "3X48_POST", "3X48_PRE", "3X49_POST", "3X49_PRE",
             "3X50_POST", "3X50_PRE", "3X51_POST", "3X51_PRE", "3X52_POST", "3X52_PRE", "3X53_POST", "3X53_PRE",
             "3X54_POST", "3X54_PRE", "3X55_POST", "3X55_PRE", "3X56_POST", "3X56_PRE", "3X57_POST", "3X57_PRE",
             "3X58_POST", "3X58_PRE")
colnames(countData) <- new_col

#Metadata
colData <- read.csv("metadata.csv")

# Making ID, condition (time), and group into factors.
colData$ID <- factor(colData$ID)
colData$condition <- factor(colData$condition)
colData$group <- factor(colData$group)
colData$sex <- factor(colData$sex)

# re-level condition. This sets the contrast as post - pre (positive logFC = upregulated at post)
colData$condition <- relevel(colData$condition, ref = "pre")

# Add a column 'ID.n' with unique numbering within each group
colData <- colData %>%
  group_by(group) %>%
  mutate(ID.n = match(ID, unique(ID)))

colData <- colData %>%
  group_by(group) %>%
  mutate(ID.n = dense_rank(ID)) %>%
  ungroup()

# Make ID.n a factor as well. 
colData$ID.n <- factor(colData$ID.n)

# Step 3: Create a DESeqDataSet Object
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ group + group:ID.n + group:condition)  #This is the design suggested by Bioconductor where subjects are nested within groups.  


# Pre-filtering low count genes. It is suggested gene count = 10 is a reasonable cutoff to filter out low count genes. 
smallestGroupSize <- 15
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,] # This resulted in reducing the detected genes from ~63000 to 17466
# Step 4: Perform Normalization and Differential Expression Analysis
dds <- DESeq(dds)
norm_count <- counts(dds, normalized = TRUE) # If normalized count data is needed. 

res <- results(dds)

resultsNames(dds) 

# Define contrast
low_post <- results(dds, contrast=list("groupLOW.conditionpost")) # LOW Post-Pre
mod_post <- results(dds, contrast=list("groupMOD.conditionpost")) # MOD Post-Pre
high_post <- results(dds, contrast=list("groupHIGH.conditionpost")) # HIGH Post-Pre
high_low_inter <- results(dds, contrast=list("groupHIGH.conditionpost", "groupLOW.conditionpost")) # HIGH - LOW
high_mod_inter <- results(dds, contrast=list("groupHIGH.conditionpost", "groupMOD.conditionpost")) # HIGH - MOD
mod_low_inter <- results(dds, contrast=list("groupMOD.conditionpost", "groupLOW.conditionpost")) # MOD - LOW

# Extract as dataframes for each contrast
low_post <- as.data.frame(low_post[order(low_post$pvalue),])
mod_post <- as.data.frame(mod_post[order(mod_post$pvalue),])
high_post <- as.data.frame(high_post[order(high_post$pvalue),])

high_low_inter <- as.data.frame(high_low_inter[order(high_low_inter$pvalue),])
high_mod_inter <- as.data.frame(high_mod_inter[order(high_mod_inter$pvalue),])
mod_low_inter <- as.data.frame(mod_low_inter[order(mod_low_inter$pvalue),])


