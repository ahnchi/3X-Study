library(WGCNA)

#Acquire dds from 3X_deseq2.R
vsd <- vst(dds, blind = FALSE)
expression_data <- vsd
rownames(expression_data) <- gsub(".*\\|", "", rownames(expression_data))
expression_data <- expression_data[, !grepl("POST", colnames(expression_data))] # Use only baseline data

#filter for ENTS that are detected in over 25%
expression_data <- assay(expression_data)
expression_data <- as.data.frame(expression_data)  # Convert matrix to data frame
expression_data$Zero_count <- rowSums(expression_data == 0)
hist(expression_data$Zero_count)
expression_data1 = expression_data[expression_data$Zero_count<0.25*45,] # filtered out low count genes
expression_data <- expression_data1 %>%
  dplyr::select(-Zero_count)### 2.0 - Check for outlier/bad genes and samples ####

datExpr_3x = data.frame(t(expression_data)) 

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# check if all genes/samples are good
gsg = goodSamplesGenes(datExpr_3x, verbose = 3) 
gsg$allOK #[1] TRUE

# Cluster samples and determine cuttoff for outliers
sampleTree = hclust(dist(datExpr_3x), method = "average")
plot(sampleTree, main = paste0("Sample clustering to detect outliers"), sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = 80, col = "red")
dev.off()

clust = cutreeStatic(sampleTree, cutHeight = 80, minSize = 4)
table(clust)
# no removal
#clust
#1 
#45 

### 3.0 - Soft thresholding and network topology analysis ####

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr_3x, powerVector = powers, verbose = 5)

# Plot the results:
sft.data = sft$fitIndices

# Scale-free topology fit index as a function of the soft-thresholding power
g1 = ggplot(sft.data, aes(Power, SFT.R.sq, label = Power) ) +
  geom_text(nudge_y = 0.05) +
  geom_point() +
  geom_hline(yintercept = 0.8, color = "red") +
  labs(
    x="Soft Threshold (power)",y="Scale Free Topology Model Fit,signed R^2",
  ) +
  theme_classic()

# Mean connectivity as a function of the soft-thresholding power
g2 = ggplot(sft.data, aes(Power, mean.k., label = Power) ) +
  geom_text(nudge_y = 100) +
  geom_point() +
  geom_hline(yintercept = 0.0, color = "red") +
  labs(
    x="Soft Threshold (power)", y="Mean Connectivity",
  ) +
  theme_classic()

plot_grid(g1, g2)

### 4.0 - Module Generation ####
# net colors is n modules
allowWGCNAThreads()

# set power decided from soft thresholding and min genes per modules
pwr = 6
ngenes = 100

# takes several hours
net = blockwiseModules(datExpr_3x,
                       power = pwr, #determined above
                       maxBlockSize = 20000, # we can set this high because we have a lot of ram. We have ~17000 genes so this will do them all at onec
                       TOMType = "signed", # unsigned was default
                       minModuleSize = ngenes,
                       reassignThreshold = 0.3,
                       mergeCutHeight = 0.3,
                       numericLabels = TRUE, # names modules as numbers isntead of colors
                       pamRespectsDendro = FALSE,
                       saveTOMs = F,
                       nThreads = 24,
                       verbose = 3)

table(net$colors)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs

######################ORA on modules #######################
module_membership = as.data.frame(net$colors)
colnames(module_membership) = 'module'
module_membership$ID = row.names(module_membership)
mm = module_membership %>%
  mutate( module = paste0("ME", module) ) %>%
  filter(module != "ME0" )
signatures = split(mm$ID, mm$module)

library(hypeR)
library(gtools)
genesets <- msigdb_gsets("Homo sapiens", "C5", "BP", clean=TRUE)
bckg <- colnames(datExpr_3x)
mhyp <- hypeR(signatures, genesets, test="hypergeometric", background=bckg)