TCGA-BRCA RNA-seq Differential Expression Analysis

A. Overview - 
For this project, an RNA sequencing differential gene expression analysis of breast cancer (TCGA-BRCA) was carried out using the R programming language.
Differential expression analysis was performed using DESeq2 on primary breast tumour and solid tissue normal samples from The Cancer Genome Atlas (TCGA).
Downstream analysis included volcano plot and top-50 gene heatmap visualization as well as Gene Ontology (GO) Biological Process enrichment analysis of differentially expressed genes (DEGs).
The project serves as an illustration of a bioinformatics workflow starting from publicly available cancer transcriptomic data, leading to differential expression analysis, visualization, and GO enrichment analysis.

B. Research Question - "Which genes and biological processes demonstrate significant differences in expression between breast tumor and normal breast tissue samples in the TCGA-BRCA dataset?"

C. Dataset - The project utilizes RNA sequencing gene expression data from:
1. The Cancer Genome Atlas – Breast Invasive Carcinoma (TCGA-BRCA)
2. The data were downloaded from the Genomic Data Commons (GDC) using the R package TCGAbiolinks.

D. Sample Selection - The following sample types were selected for the balanced exploratory analysis:
1. 50% Primary Tumour
2. 50% Solid Tissue Normal
A total of 100 samples were analyzed across two groups. A fixed random seed of '123` was used during sampling to make a reproducible analysis .

E. Analysis Workflow - The analysis followed the workflow below:

TCGA-BRCA RNA-seq Data -> Sample Selection (50 Tumor + 50 Normal) -> Gene Count Matrix (60,660 genes) -> Low-count Filtering (26,284 genes retained) -> DESeq2 Differential Expression -> Significant DEGs (padj < 0.05) (|log2FC| > 1) -> Visualization (Volcano Plot - Top-50 Gene Heatmap) -> Ensembl → Entrez ID Mapping -> GO Biological Process Enrichment -> Biological Interpretation

F. Tools and Packages - The analysis was performed in R using the following packages:
1. TCGAbiolinks – TCGA/GDC data retrieval and preparation
2. SummarizedExperiment – Genomic assay data handling
3. DESeq2 – Differential gene expression analysis
4. EnhancedVolcano – Volcano plot visualization
5. pheatmap – Heatmap visualization
6. clusterProfiler – Functional enrichment analysis
7. org.Hs.eg.db – Human gene annotation
8. enrichplot – Enrichment visualization
   
G. Differential Expression Analysis
Gene expression counts were analysed using the DESeq2 workflow.
Genes with very low read counts were filtered before differential expression analysis. A gene was retained if it had a count of at least 10 in at least 10 samples.

H. Gene filtering - Initial Stage had 60,660 genes and After low-count filtering stage there were	26,284 genes
Differential expression was evaluated using the contrast: Tumor vs Normal
Genes were considered significant using: Adjusted p-value (padj) < 0.05 AND |log2 Fold Change| > 1

I. Differential expression results - A total of 6,853 significant differentially expressed genes (DEGs) were identified using the above thresholds.
1. The complete DESeq2 results are provided in: "results/DEG_results.csv"
2. The filtered list of significant DEGs is provided in: "results/significant_DEGs.csv"

J. Visualizations
a. Volcano Plot - The volcano plot displays the relationship between statistical significance and magnitude of gene expression change.
1. X-axis: log2 Fold Change
2. Y-axis: adjusted p-value
3. Significance threshold: padj < 0.05
4. Fold-change threshold: |log2FC| > 1
The plot labels the top 20 genes ranked by adjusted p-value.

b. Heatmap of Top Differentially Expressed Genes - A heatmap was generated using the 50 most significant genes based on adjusted p-value.
The expression values were variance-stabilized using the DESeq2 vst() transformation and row-scaled before visualization.
Sample annotations indicate:
1. Normal tissue
2. Primary tumour tissue

K. Gene Ontology Enrichment Analysis - To investigate the biological functions represented by the significant DEGs, Gene Ontology enrichment analysis was performed using:
1. clusterProfiler
2. org.Hs.eg.db

Significant Ensembl gene identifiers were mapped to Entrez Gene IDs before enrichment analysis.

a. Gene ID mapping - Of the 6,853 significant DEGs, 5,467 genes were successfully mapped to Entrez Gene IDs and used for GO enrichment analysis.

b. GO Biological Process results - The analysis identified 1,068 enriched Gene Ontology Biological Process terms using: "Benjamini-Hochberg adjusted p-value < 0.05"
Among the most strongly enriched biological processes were:

1. Extracellular structure organization
2. External encapsulating structure organization
3. Extracellular matrix organization
4. Regulation of membrane potential
5. Regulation of hormone levelsBlood circulation
6. Signal release
7. Regulation of trans-synaptic signaling
8. Skeletal system development
9. Vascular processes in the circulatory system
10. Humoral immune response
11. Antimicrobial humoral response
12. Axon development
13. Axonogenesis

The complete enrichment results are available in: "results/GO_BP_enrichment.csv"
GO enrichment visualization

L. Key Findings - The exploratory analysis identified substantial differences in gene expression between breast tumour and normal tissue samples.
Main results
1. 100 total samples analysed
2. 50 tumour samples
3. 50 normal samples
4. 60,660 genes initially available
5. 26,284 genes retained after low-count filtering
6. 6,853 significant DEGs
7. 5,467 significant genes mapped to Entrez IDs
8. 1,068 enriched GO Biological Process terms

The strongest enrichment was observed in processes related to extracellular structure and extracellular matrix organization, along with processes involving signalling, vascular activity, hormone regulation, and immune-related functions.
These findings provide a functional overview of biological processes represented among genes showing differential expression between tumour and normal breast tissue.

M. Repository Structure
TCGA-BRCA-RNA_seq-Differential-Expression/

│

├── README.md

├── analysis.R

│

├── figures/

│  ├── volcano_plot.png

│  ├── heatmap_top50.png

│  └── GO_BP_dotplot.png

│

└── results/

├── DEG_results.csv

├── significant_DEGs.csv

├── GO_BP_enrichment.csv

└── sessionInfo.txt

N. Reproducibility - The complete R analysis pipeline is provided in: "analysis.R"
The script contains the main steps used for:
1. Querying TCGA-BRCA data
2. Selecting tumour and normal samples
3. Preparing the RNA-seq count matrix
4. Filtering low-count genes
5. Performing DESeq2 differential expression analysis
6. Generating the volcano plot
7. Generating the top-gene heatmap
8. Performing GO Biological Process enrichment
9. Saving analysis results and figures
The R package versions and system information used during the analysis are provided in: "results/sessionInfo.txt"

O. Limitations - This project is intended as an exploratory bioinformatics analysis and has several limitations.
1. Only 100 TCGA-BRCA samples were analysed rather than the complete available cohort.
2. Samples were selected as a balanced subset of 50 tumour and 50 normal samples.
3. The analysis is based on observational transcriptomic data.
4. Differential expression does not by itself establish biological causation.
5. The identified genes should not be interpreted as clinically validated biomarkers.
6. Further validation using independent datasets and experimental approaches would be required to establish biological or clinical significance.
7. GO enrichment provides functional associations but does not demonstrate that individual biological processes directly cause breast cancer.

P. Future Directions - The current workflow can be extended in several ways:
1. Analyse the complete TCGA-BRCA cohort.
2. Validate significant genes using an independent breast cancer dataset.
3. Investigate breast cancer molecular subtypes.
4. Perform pathway enrichment using KEGG, Reactome, or MSigDB.
5. Construct protein-protein interaction networks.
6. Perform survival analysis for candidate genes.
7. Apply feature-selection and machine-learning approaches for candidate biomarker identification.
8. Extend the analysis toward multi-omics integration.
9. Investigate tumour heterogeneity using single-cell RNA-seq datasets.

Q. Conclusion
This project demonstrates a complete introductory RNA-seq bioinformatics workflow using TCGA-BRCA data. Starting from publicly available transcriptomic data, the analysis applies DESeq2-based differential expression, statistical filtering, visualization, gene annotation, and Gene Ontology enrichment analysis to identify genes and biological processes associated with differences between breast tumor and normal tissue. The project provides a foundation for further computational analysis and demonstrates practical skills in R programming, cancer genomics, transcriptomic data analysis, statistical analysis, data visualization, and functional enrichment analysis.

R. Author

Ishita Singh

B.Tech Bioinformatics

GitHub: @ishitasinghbiocode

