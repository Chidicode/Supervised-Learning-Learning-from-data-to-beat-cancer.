# Part1: Setup and Installation........


## Install and load Required Packages


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
BiocManager::install("edgeR")
BiocManager::install("EDASeq")
BiocManager::install("SummarizedExperiment")
BiocManager::install("biomaRt")
install.packages("dplyr")
install.packages('gplots')



## Load required libraries


library(TCGAbiolinks)            #For accessing and querying TCGA data 
library(edgeR)                   #For differential expression analysis
library(EDASeq)                  #Exploratory data analysis and normalization
library(SummarizedExperiment)    #Access assay data
library(biomaRt)                 #Access BioMart databases      
library(dplyr)                   #Manage and manipulate data
library(gplots)                  #Data visualize and heatmap generate



#Part 2: Data Preprocessing........


## Collect the list of cancer projects from TCGA databases

  -Lung Adenocarcinoma (LUAD) dataset select

gdcprojects <- getGDCprojects()
getProjectSummary('TCGA-LUAD')
?GDCquery



##Download data

### Build a query to retrieve gene expression data 



luadQ <- GDCquery(project = 'TCGA-LUAD',
                  data.category = 'Transcriptome Profiling',
                  experimental.strategy = "RNA-Seq",
                  workflow.type = "STAR - Counts",
                  access = "open",
                  data.type = "Gene Expression Quantification",
                  sample.type = c("Primary Tumor", "Solid Tissue Normal"))


###Downlaod query data

```{r}
GDCdownload(luadQ)
luad.data <- GDCprepare(luadQ)

```

## Data structure

```{r}
head(luad.data)
View(luad.data)
colnames(luad.data)       
nrow(luad.data)             
ncol(luad.data)
```


## Explore metadata information

```{r}
luad.data$race
luad.data$tumor_descriptor
luad.data$barcode
luad.data$sample_type
luad.data$sample_id

```


## Create a simple metadata for this analysis

```{r}

Metadata <- data.frame("barcode"= luad.data$barcode,
                       "race" = luad.data$race,
                       "tumor_type" = luad.data$tumor_descriptor,
                       "sample" = luad.data$sample_type,
                       "sample_id" = luad.data$sample_id)

```


## Save a metadata to a CSV file

```{r}
# Save the metadata to a CSV file
write.csv( Metadata, "TCGA_LUAD_metadata.csv", row.names = FALSE)
```


## Select unstranded dataset

```{r}

luad.raw.data <- assays(luad.data)
dim(luad.raw.data$unstranded)

```



##Downsize dataset 

- 20 select for primary tumor and 20 for Solid Tissue Normal data

```{r}

selectedBarcodes <- c(subset(Metadata, sample == "Primary Tumor")$barcode[c(1:20)],
                      subset(Metadata, sample == "Solid Tissue Normal")$barcode[c(1:20)])
```


## Select unstranded dataset

```{r}
selectedData <- luad.raw.data$unstranded[, c(selectedBarcodes)]
dim(selectedData)
```


#Part 3: Normalization........

## Data normalization and filtering

```{r}
normData <- TCGAanalyze_Normalization(tabDF = selectedData, geneInfo = geneInfoHT, method = "geneLength")



```



##Filtering

```{r}

fildata <- TCGAanalyze_Filtering(tabDF = normData,
                                 method = "quantile",
                                 qnt.cut = 0.25)


```



##Dimension of filtering data

```{r}
dim(fildata)
```



#Part 4: Differential expression analysis........


```{r}

results <- TCGAanalyze_DEA(mat1 = fildata[, c(selectedBarcodes)[1:20]],
                           mat2 = fildata[, c(selectedBarcodes)[21:40]],
                           Cond1type = "Primary Tumor",
                           Cond2type = "Solid Tissue Normal",
                           pipeline = "edgeR",
                           fdr.cut = 0.01,
                           logFC.cut = 1)

results.level <- TCGAanalyze_LevelTab(results, "Primary Tumor", "Solid Tissue Normal",
                                      fildata[, c(selectedBarcodes)[1:20]],
                                      fildata[, c(selectedBarcodes)[21:40]])

```



#Part 5: Data visualization........


##Volcano plot

```{r}
# Already have 'results' with logFC and p-values (p.adjust)

# Set thresholds for significance
logFC_cutoff = 1      # Log Fold Change threshold
pvalue_cutoff = 0.01  # Adjusted p-value threshold

# Create a new column to classify significance
results$threshold <- ifelse(results$logFC > logFC_cutoff & results$PValue < pvalue_cutoff, "Upregulated",
                            ifelse(results$logFC < -logFC_cutoff & results$PValue < pvalue_cutoff, "Downregulated", "Not Significant"))

# Assign colors based on thresholds
results$color <- ifelse(results$threshold == "Upregulated", "red",
                        ifelse(results$threshold == "Downregulated", "blue", "gray"))

# Negative log10 p-values for the y-axis
results$log10_pvalue <- -log10(results$PValue)

# Volcano plot using base R's plot function
plot(results$logFC, results$log10_pvalue,
     pch = 16,                      # Solid points
     col = results$color,            # Colors based on significance
     xlab = "Log Fold Change (logFC)",
     ylab = "-log10 Adjusted P-value",
     main = "Volcano Plot: Upregulated vs Downregulated",
     cex = 1.2)                      # Adjust point size

# Add a legend
legend("topright",
       legend = c("Upregulated", "Downregulated", "Not Significant"),
       col = c("red", "blue", "gray"),
       pch = 16)

```


##Heatmap Generation


### DEA with treatment levels

```{r}

#Differential expression analysis with treatment levels
results.level <- TCGAanalyze_LevelTab(results, "Primary Tumor", "Solid Tissue Normal",
                                      fildata[, c(selectedBarcodes)[1:20]],
                                      fildata[, c(selectedBarcodes)[21:40]])
```


###Data structure

```{r}
head(results.level)
dim(results.level)
```


###Fildata with result levels

```{r}
heat.data <- fildata[rownames(results.level),]
```


###Select plot color

```{r}

cancer.type <- c(rep("Primary Tumor", 20), rep("Solid Tissue Normal", 20))
ccodes <- c()

for (i in cancer.type) {
  if(i == "Primary Tumor") {
    ccodes <- c(ccodes, "red")
  } else{
    ccodes <- c(ccodes, "blue")
  }
}
```



##Generate heatmap with column clustering only

```{r}

# Save the plot to a PNG file with larger dimensions
png("heatmap_output.png", width = 1000, height = 800)
heatmap.2(
  x = as.matrix(heat.data),
  col = hcl.colors(10, palette = "Blue-Red 2"),
  Rowv = FALSE,              # Disable row clustering
  Colv = TRUE,               # Enable column clustering
  scale = "row",
  sepcolor = "black",
  trace = "none",
  key = TRUE,
  dendrogram = "col",         # Cluster columns only (no rows)
  cexRow = 0.5, 
  cexCol = 0.5,              # Define column text size
  main = "Heatmap of Primary Tumor vs Solid Normal Tissue",
  na.color = "black",
  ColSideColors = ccodes
)
```


#Part 6: Subset genes........



##Preview logFC values

```{r}
head(results.level$logFC)  
```


##Set cutoff thresholds for fold change and p-value

```{r}
# Set thresholds for significance
logFC_cutoff = 1      # Log Fold Change threshold
pvalue_cutoff = 0.01  # Adjusted p-value threshold
```


##Create a new column to classify significance

```{r}
results$threshold <- ifelse(results$logFC > logFC_cutoff & results$PValue < pvalue_cutoff, "Upregulated",
                            ifelse(results$logFC < -logFC_cutoff & results$PValue < pvalue_cutoff, "Downregulated", "Not Significant"))
```


##Subset significantly upregulated genes

```{r}
upregulated_genes <- results %>%
  filter(logFC > logFC_cutoff & PValue < pvalue_cutoff)
```


##Subset significantly downregulated genes

```{r}
downregulated_genes <- results %>%
  filter(logFC < logFC_cutoff & PValue < pvalue_cutoff)
```


##Print the results

```{r}
print(upregulated_genes)
print(downregulated_genes)
```

##Subset significantly upregulated and downregulated genes from results with levels

```{r}
upreg_genes <- rownames(subset(results.level, logFC > 1))
downreg_genes <- rownames(subset(results.level, logFC< -1))
```


##Print the results

```{r}
print(upreg_genes)
print(downreg_genes)
```


##Length 

```{r}
length(upreg_genes)
length(downreg_genes)
```


##Save file

```{r}
write.csv(upreg_genes, "upregulated_genes.levels.csv")
write.csv(downreg_genes, "downregulated_genes.levels.csv")
```


#Part 7: Gene Annotation........


- BioMart is a powerful, web-based data management and analysis tool that integrates data from multiple biological databases, including Ensembl (genome data) and others


##Connect to the Ensembl BioMart database for human genes

```{r}
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
```


##Convert ensembel IDs to gene IDs (upregulated genes) using biomaRt

```{r}


mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

upreg_genes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                     filters = "ensembl_gene_id",
                     values = upreg_genes,
                     mart = mart)$hgnc_symbol
```


##Data structure

```{r}
head(upreg_genes)
length(upreg_genes)
```


####Convert ensembel IDs to gene IDs (downregulated genes) using biomaRt

```{r}
downreg_genes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                       filters = "ensembl_gene_id",
                       values = downreg_genes,
                       mart = mart)$hgnc_symbol


```


##Data structure

```{r}
head(downreg_genes)
length(downreg_genes)
```


##Print results

```{r}
print(upreg_genes)
print(downreg_genes)
```


##Save file

```{r}
write.csv(upreg_genes, "upregulated_genes_ei.csv", row.names = FALSE)
write.csv(downreg_genes, "downregulated_genes_ei.csv", row.names = FALSE)
```






#Part 8: Functional Enrichment Analysis........

##Enrichment analysis

```{r}
upreg_EA <- TCGAanalyze_EAcomplete(TFname = "upregulated",
                                   upreg_genes)
downreg_EA <- TCGAanalyze_EAcomplete(TFname = "downregulated",
                                     downreg_genes)
```



#Visualize results with barplot (upreguated genes)

```{r}
TCGAvisualize_EAbarplot(tf = rownames(upreg_EA$ResBP),
                        GOBPTab = upreg_EA$ResBP,
                        GOCCTab = upreg_EA$ResCC,
                        GOMFTab = upreg_EA $ ResMF,
                        PathTab = upreg_EA $ResPat,
                        nRGTab = upreg_genes,
                        nBar =5,
                        text.size = 2,
                        fig.width =30,
                        fig.height =15)
```


#Visualize results with barplot (downreguated genes)

```{r}
TCGAvisualize_EAbarplot(tf = rownames(downreg_EA$ResBP),
                        GOBPTab = downreg_EA$ResBP,
                        GOCCTab = downreg_EA$ResCC,
                        GOMFTab = downreg_EA $ ResMF,
                        PathTab = downreg_EA $ResPat,
                        nRGTab = downreg_genes,
                        nBar =5,
                        text.size = 2,
                        fig.width =30,
                        fig.height =15)
```
