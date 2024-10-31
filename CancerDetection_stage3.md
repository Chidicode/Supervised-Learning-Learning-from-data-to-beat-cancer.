
# Identify Potential Biomarkers for Cancer Detection by Using Differential Expression and Machine Learning Models



## Objective

The main objective of this study is to identify potential  biomarkers for early lung cancer detection by analyzing the TCGA dataset with differential gene expression data and predict cancer by using machine learning models. The aim is to find key genes that are significantly expressed in early-stage lung cancers compared to normal tissue. 


## Introduction

Lung adenocarcinoma represents a histopathological subtype of non-small cell lung cancer (NSCLC) and has various forms, such as solid, acinar, lepidic, papillary, and micropapillary subtypes [1]. This subtype comprises nearly half of all lung cancer cases [2]. It originates from the glandular epithelial cells that line the alveoli of the lungs and is characterized by the formation of glandular structures or mucin production[3].

It has many risk factors, including smoking, air pollution, genetic predisposition, occupational hazards, and exposure to substances such as silica, asbestos, diesel exhaust, and heavy metals [4]. In the recent advancements of bioinformatics, novel approaches have been developed for identifying potential drug targets for this cancer subtype [5].


## Dataset

The dataset is taken from The Cancer Genome Atlas (TCGA) website that consists of more than 30 different types of cancer datasets (https://portal.gdc.cancer.gov/).

### Download data

TCGAbiolinks is an R/Bioconductor package designed to facilitate the retrieval, analysis, and integration of data from The Cancer Genome Atlas (TCGA) and other Genomic Data Commons (GDC) resources.


### Data type

The query searches the  TCGA-LUAD (Lung Adenocarcinoma) dataset with open-access RNA-seq gene expression data using  STAR - Counts refers to the Spliced Transcripts Alignment to a Reference (STAR) alignment method, which is a popular tool for aligning RNA-Seq reads to a reference genome and pulls raw reads counts from primary tumor tissue (tumor) and solid tissue normal (normal tissue adjacent to the tumor). 


## Main Workflow

![Snip20240924_3](https://github.com/user-attachments/assets/508de50f-4f0c-4afc-b3f9-dba0ac4e8c14)






## Conclusion 

 - DEG analysis revealed a set of upregulated and downregulated genes in lung adenocarcinoma (LUAD), providing insight into the molecular mechanisms underlying cancer progression

 - The subsequent enrichment analysis of these genes offers valuable information on the biological processes, molecular functions, and pathways that are potentially disrupted during tumor development

 - ML models can be trained on gene expression data to  identify key biomarkers for lung cancer that are linked to prognosis and  targeted therapy


## References

1. Liang, Y., Xie, Y., Yu, H., Zhu, W., Yin, C., Zhang, X., & Dong, Z. (2023). Clinical significance of TMEM229A Q200del mutation in lung adenocarcinoma.
2. Shiba‐Ishii, A. (2021). Significance of stratifin in early progression of lung adenocarcinoma and its potential therapeutic relevance. Pathology international, 71(10), 655-665.
3. Herbst, R. S., Morgensztern, D., and Boshoff, C. "The biology and management of non-small cell lung cancer." Nature, vol. 553, 2018, pp. 446-454.
4. Verma, G., Ravichandar, S., Abraham, E. A., & Sampath, S. (2024). Myriad Presentation of Adenocarcinoma Lung. Journal of Evolution of Medical and Dental Sciences, 55-57.
5. Xie, H., Zhang, J. F., & Li, Q. (2021). Identification and analysis of genes associated with lung adenocarcinoma by integrated bioinformatics methods. Annals of human genetics, 85(3-4), 125-137.
