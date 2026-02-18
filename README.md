# scRNAseq-Analysis
Comprehensive single-cell RNA-seq analysis of GSE279086 investigating altered kidney oxidative metabolism in young adults with Type 1 Diabetes.

# Step 1: Raw Data Download
01_download.sh script automatically downloads raw single-cell gene expression matrices for dataset GSE279086 from the GEO FTP server. Each GSM accession is mapped to its corresponding sample ID. All associated matrix, feature, and barcode files are retrieved in a structured directory format. The script ensures reproducibility by organizing files per sample and skipping already downloaded files.

# Step 2: Raw File Cleanup and Standardization
02_cleanup_geo_files.sh cleans the downloaded GEO sample folders by removing pre-processed files, ensuring that only raw count matrices are used for downstream analysis. It standardizes file naming conventions to the 10x Genomics format (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz) for seamless compatibility with Scanpy and Seurat. This step guarantees consistency and reproducibility across all samples.

# Step 3: Creation of Seurat Objects and Quality Metric Annotation (03_10XtoSeurat.Rmd)
In this step, raw single-cell gene expression matrices for each GSM sample are imported into Seurat, a widely used framework for scRNA-seq analysis. Each sample directory containing 10x Genomics-formatted files (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz) is read independently to preserve sample-level identity prior to downstream integration.

For every sample, a Seurat object is constructed with minimal filtering thresholds (min.cells = 3, min.features = 200) to remove lowly expressed genes and likely empty droplets, ensuring biologically meaningful cellular profiles are retained. Sample identifiers are stored as metadata, enabling traceability across quality control, integration, and differential expression analyses.

To assess cell quality, mitochondrial gene content (percent.mt) and ribosomal gene content (percent.rb) are computed for each cell. Elevated mitochondrial transcript abundance is commonly associated with stressed or dying cells, while ribosomal gene proportions provide insight into translational activity and potential technical bias. These metrics form the foundation for downstream quality filtering and normalization decisions.

Finally, missing or undefined QC values are explicitly handled to maintain data integrity and prevent downstream computational errors. This step results in a list of fully annotated Seurat objects, one per sample, ready for systematic quality control, normalization, and batch-aware integration.

# Step 4: Quality Control, Integration, and Dimensionality Reduction (04_qc_batchcorrection_GSE279086)
This step performs stringent quality control (QC), sample integration, and dimensionality reduction to ensure that downstream biological interpretation is driven by true cellular signal rather than technical noise or batch effects.

# a) Input Management and Library Setup
All sample-level Seurat objects are loaded from the input directory and managed collectively to enable multi-sample analysis. A comprehensive set of bioinformatics packages is used, covering:
Data handling (Seurat, SingleCellExperiment)
Visualization (ggplot2, dittoSeq)
Batch correction (Harmony, batchelor)
Format interoperability (SeuratDisk, zellkonverter)
This ensures compatibility across R- and Python-based single-cell ecosystems.

# b) Sample Integrity Checks and Merging
Before merging, each Seurat object is examined for duplicate gene identifiers, which can cause artifacts during normalization and integration. As no duplicated genes were detected, all samples are safely merged into a single Seurat object, while preserving original GSM/sample identities in metadata.The merged dataset initially contains 181,842 cells across 28,317 genes, representing the full raw cellular landscape prior to QC.

# c) Metadata Harmonization
To enable biologically meaningful comparisons, phenotypic metadata (Control vs Type 1 Diabetes) is retrieved directly from GEO using GEOquery and mapped to each cell via GSM identifiers. This allows all QC diagnostics, clustering, and downstream analyses to be interpreted in a disease-aware context rather than purely technical space.

# d) Quality Metric Computation (Cell-Level Health Assessment)
For every cell, core QC metrics are computed:
nCount_RNA: total RNA molecules captured (library size)
nFeature_RNA: number of genes detected (cellular complexity)
percent.mt: mitochondrial RNA fraction (cell stress / apoptosis indicator)
percent.rb: ribosomal RNA fraction (technical dominance / low information content)

These metrics collectively distinguish healthy single cells from_Empty droplets, Dying or stressed cells, Doublets or multiplets

# e) Visual QC Diagnostics (Pre-Filtering)
QC thresholds are critical because inappropriate filtering can either retain technical artifacts or remove biologically valid cells. Violin and density scatter plots are used to visualize the distribution of QC metrics across samples, allowing clear identification of outliers such as empty droplets, stressed cells, or potential doublets. Thresholds are chosen at points where cells deviate sharply from the main population, ensuring that filtering is data-driven, reproducible, and biologically justified

# f) Multi-Layered Cell Filtering Strategy
Cells are filtered using multi-criteria approach:
Gene count filtering removes empty droplets (<200 genes) and likely doublets (>7000 genes)
Mitochondrial filtering removes metabolically compromised cells (>15% mtRNA)
Ribosomal filtering removes low-information cells (>10% rRNA)
Mahalanobis distance filtering detects subtle global outliers in RNA content, capturing abnormal cells missed by univariate thresholds

This conservative strategy dramatically improves data quality, reducing the dataset to 6,327 high-confidence cells.

# g) Post-QC Validation and Reporting
Cell retention is quantified per sample, and a QC summary table is exported. Post-QC violin and scatter plots confirm successful removal of technical artifacts while preserving biological diversity.
Only protein-coding genes are retained, eliminating non-informative transcripts prior to normalization.

# h) Normalization, Feature Selection, and PCA
The filtered dataset undergoes:
Log-normalization to account for library size differences.
Highly variable gene selection (2,500 genes) to focus on biologically informative variation.
Scaling and PCA (100 components) to compress high-dimensional gene expression into interpretable latent space.
An elbow plot and cumulative variance analysis guide objective PC selection, with the top 35 PCs explaining ~83% of total variance.

# i) Clustering and UMAP (Unintegrated Data)
Cells are clustered using a graph-based approach and visualized via UMAP using PCA space. This provides a baseline view of biological structure before batch correction, highlighting sample-driven segregation.

# j) Batch Correction Using Harmony
To correct for sample-specific technical effects while preserving biological signal, Harmony integration is applied using sample identity as the batch variable. Clustering and UMAP are repeated in Harmony space, producing:
Improved cross-sample mixing
Reduced batch-driven separation
More biologically coherent clusters

# k) Quantitative Batch Effect Evaluation
Batch correction effectiveness is assessed using a k-nearest neighbor batch-mixing metric, quantifying how frequently neighboring cells originate from the same sample. Comparison before and after Harmony demonstrates a substantial reduction in batch bias, validated both visually and quantitatively.

# l) Interoperability and Output Export
Final objects are saved in multiple formats:
.rds for Seurat workflows
.h5seurat and .h5ad for cross-platform compatibility
UMAP coordinates for external visualization or downstream modeling
This ensures full reproducibility, portability, and extensibility.

# Step 5: Cell Type Annotation, Composition Analysis, and Marker Gene Identification (05_celltypist_GSE279086)
In this step, the batch-corrected single-cell dataset is annotated at the cell-type level using CellTypist, followed by cell composition analysis and identification of cluster-specific marker genes. Cell-type annotation using CellTypist.

# a) Cell-type annotation using CellTypist
The Harmony-integrated AnnData object is first aligned with Seurat-derived UMAP coordinates to ensure consistent visualization across R and Python workflows. Log-normalized expression values are explicitly assigned to the AnnData.X matrix, as required by CellTypist’s logistic regression–based classifiers. A two-step Leiden clustering strategy is applied to generate fine-grained clusters, which are then used as input for over-clustered majority voting. This approach improves annotation robustness by reducing noise from individual cell-level misclassification and assigning biologically coherent cell-type labels using a pretrained Adult Human Kidney reference model.

# b) Condition-specific visualization
UMAPs are generated separately for each experimental condition, allowing direct comparison of cell-type distributions between control and Type 1 Diabetes samples. This helps identify condition-specific enrichment or depletion of kidney cell populations while preserving the shared embedding space.

# c) Cell-type proportion analysis
Annotated cells are quantified per condition to calculate relative cell-type proportions, enabling a compositional comparison between experimental groups. Cell types are ordered based on average abundance, and bar plots are generated to visually highlight shifts in cellular composition that may reflect disease-associated remodeling rather than transcriptional changes alone.

# d) Identification of unique cluster marker genes
To characterize transcriptional identity at high resolution, differential expression analysis is performed using a Wilcoxon rank-sum test. Marker genes are filtered based on adjusted p-values, log fold-change thresholds, and minimum cell counts to ensure statistical reliability. Only genes uniquely enriched in a single cell-type–condition group are retained, mimicking Seurat’s FindAllMarkers logic while emphasizing specificity. This enables precise biological interpretation of cell states across conditions.

# e) Marker gene visualization and final data export
Top unique marker genes are visualized using dot plots, summarizing both expression level and detection frequency across annotated groups. Finally, the fully annotated AnnData object is saved, providing a reproducible and analysis-ready dataset for downstream differential expression, pathway analysis, or integration with spatial and clinical data.

# Step 6: Cell-type–resolved Differential Expression and Pathway Analysis (06_DEGs_Pathways_GSE279086)
# a) Cross-platform data integration (Python → R)
The CellTypist-annotated AnnData object is imported into R using zellkonverter and converted into a Seurat object. This enables seamless interoperability between Python-based annotation workflows and R-based statistical modeling, while preserving raw counts, log-normalized expression, metadata, and cell-type labels. 

# b) Cell-type–specific differential expression using MAST
Differential gene expression (DGE) analysis is performed independently for each annotated cell type, comparing Type 1 Diabetes versus Control cells. The MAST model is used because it explicitly accounts for:
Zero inflation inherent to scRNA-seq data
Differences in cellular detection rates
Continuous expression among detected genes

Technical covariates (nCount_RNA, percent.mt) are included as latent variables to control for sequencing depth and mitochondrial stress, improving biological specificity. Cell types with insufficient representation are excluded to ensure statistical robustness.

# c) DEG summary across cell types
For each cell type, significantly upregulated and downregulated genes are quantified, producing a global overview of which kidney cell populations exhibit the strongest transcriptional dysregulation in Type 1 Diabetes. This enables prioritization of biologically responsive cell types for downstream interpretation.

# d) Pathway-level interpretation using Reactome GSEA
Gene Set Enrichment Analysis (GSEA) is applied using Reactome pathways. Genes are ranked by log₂ fold change, preserving directionality and avoiding arbitrary significance cutoffs. This approach captures coordinated pathway-level shifts, even when individual genes show modest effects.Gene symbols are mapped to ENTREZ IDs to ensure compatibility with curated pathway databases. GSEA is performed separately for each cell type, enabling identification of cell-type–specific metabolic, signaling, and disease-associated pathways altered in Type 1 Diabetes.

# e) Visualization and export
Top enriched pathways per cell type are visualized using dot plots, integrating normalized enrichment score (NES), statistical significance, and pathway size. All DEG tables, pathway results, plots, and the fully annotated Seurat object are saved, ensuring reproducibility and enabling downstream mechanistic or comparative analyses.

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# RESULTS and INTREPRETATION
# UMAP
Refer: umap_celltypes_ct_Type 1 Diabetes and Control.pdf
UMAP visualization of integrated single-cell RNA-seq data revealed clear segregation of major kidney cell types in both Control and Type 1 Diabetes (T1D) samples, including proximal tubules (PT), thick and thin ascending limbs (TAL/DTL), distal convoluted tubules (DCT), endothelial subtypes, podocytes, pericytes, and immune cells.

The preservation of distinct and well-separated clusters across conditions indicates successful data integration and robust cell-type annotation. The overall cellular architecture is maintained in T1D kidneys, suggesting that early diabetic pathology does not involve widespread loss of kidney cell lineages. Instead, disease-associated effects are more likely driven by transcriptional and functional changes within existing cell types rather than the emergence or disappearance of specific populations.

# Cell-type proportion changes between Control and T1D kidneys
Refer: Cellproportions_barplot and Markergenes.pdf
Comparative analysis of cell-type proportions revealed marked compositional remodeling in T1D kidneys. While endothelial peritubular capillary cells (EC-PTC) remained the most abundant endothelial population in both conditions. Notable shifts were observed across specific tubular compartments, including the proximal convoluted tubule (PC), descending thin limb (DTL), and medullary thick ascending limb (M-TAL), as well as stromal and vascular-associated populations such as vascular smooth muscle cells/pericytes and endothelial subtypes (EC-DVR, EC-AEA, EC-GC). In contrast, podocytes and NK cells were reduced in Type 1 Diabetes.

The persistence of EC-PTC abundance in T1D is biologically significant, as these cells are critical for oxygen delivery and metabolic support of highly energy-demanding tubular segments. Their retention suggests that impaired kidney oxidative metabolism in T1D is not due to capillary loss but likely reflects altered metabolic function within resident cells. Concurrent expansion of tubular segments and stromal cells points toward adaptive or stress-associated remodeling, whereas the reduction in podocytes is consistent with early glomerular vulnerability in diabetic kidney disease.

Marker gene expression confirmed robust cell-type annotation across Control and Type 1 Diabetes kidneys.

# Cell type specific DEG
Refer: DEGsummary_table.csv
Cell-type–specific differential expression analysis using MAST showed that gene expression changes in Type 1 Diabetes kidneys occur in selected cell populations rather than across all cell types. Principal cells displayed the largest number of differentially expressed genes, indicating strong transcriptional responses to diabetic conditions. Endothelial cell subtypes, particularly peritubular capillary and descending vasa recta endothelial cells, showed moderate but consistent changes, suggesting vascular adaptation in diabetes. In contrast, most tubular segments and podocytes showed little to no differential expression, indicating that their transcriptional programs remain largely stable at this stage of disease. Overall, these results highlight that molecular changes in T1D kidneys are cell-type–specific and justify cell-type–resolved analysis.

# REACTOME GSEA
Refer: REACTOME GSEA.pdf (PC, EC_PTC, EC_DVR, VSMC/P)
Reactome GSEA identified significant pathway enrichment in a subset of cell types, including principal cells (PC), peritubular capillary endothelial cells (EC_PTC), descending vasa recta endothelial cells (EC_DVR), and vascular smooth muscle/pericytes (VSMC/P). These cell types exhibited sufficient and coordinated transcriptional changes to support pathway-level enrichment. In contrast, other tubular segments and podocytes showed minimal differential expression, resulting in insufficient signal for robust pathway enrichment. This indicates that early Type 1 Diabetes–associated molecular alterations are concentrated within specific epithelial and vascular compartments rather than uniformly affecting all kidney cell types.

The Reactome GSEA results demonstrate that Type 1 Diabetes–associated molecular changes are spatially and functionally organized, beginning in metabolically active epithelial cells (PCs), extending to adjacent microvascular endothelial compartments (EC_PTC, EC_DVR), and finally involving vascular support cells (VSMC/P). Importantly, these changes reflect adaptive transcriptional remodeling linked to altered oxidative metabolism and vascular signaling, rather than widespread cell damage or loss. This reinforces the concept that early diabetic kidney disease is driven by cell-type–specific metabolic and microvascular stress.

