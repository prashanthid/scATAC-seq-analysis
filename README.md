# scATAC-seq-analysis
## Typical data analysis workflow
### 1. Pre-processing of sequencing reads
- Purpose: converting the sequencing machine format and removing adaptors.
- Deduplexing: Illumina's bcl2fastq
- Adaptor removal: Trimmomatic
- Mapping: Bowtie2 or BWA
- Sort: Samtools
### 2. Quality control
- Purpose: barcodes corresponding to low quality cells or doublets must be filtered out.
- QC associated with single cell techniques: low count depth (low quality reads) and high count depth (doublets).
- QC linked with scATAC-seq: fraction of reads in peaks (FRiP), ratio of reads in promoter regions, ratio of reads in blacklist sites, or enrichment of transcription start sites (TSS) are often used for barcode selection. Barcodes that do not show nucleosomal banding patterns that are unique to high-quality ATAC-seq data and features (e.g., peaks) that are located in blacklist regions or house-keeping genes are filtered out.
### 3. Cell-by-feature matrix formation
- Purpose: QC passed cells are converted to cell-by-feature matrix for downstream analysis. The genomic regions from raw peak reads and annotation of the defined regions using regulatory elements varies method to method and hence, the matrix.
- Defining genomic regions by the use of sample-specific information: utilization of bulk ATAC-seq peaks from public data or analyzing aggregated or merged peaks from scATAC-seq data. MACS2 is utilized for peak identification.
- Defining genomic regions by fixed-size bins or window.
- cell-by-feature matrix constructed by considering regulatory elements, such as TF motifs and TSS. motifs and k-mers for TF binding are specific for cell types, cell type annotation based on the information. The genomic regions are annotated with either the known TF motifs from public databases, such as cisBP, JASPAR, and HOMER, or k-mers for unsupervised annotation using motifmatchr.
### 4. Batch correction and data integration
- Purpose: Batch effects can occur by differences in experimenters, sample preparation protocols, sample harvest time, sequencing lanes, and sequencing technologies.
- Batch effects of single-cell omics data can be more systematically corrected with data integration approaches based on non-linear algorithms. These methods assume that all batches share at least one cell type with another and differences between batches are smaller than those between cell types.
- Harmony, Seurat and scVI: best trade-off between batch effect removal and conservation of biological variation.
### 5. Data transformation
- Purpose: Peak reads from a single cell have been reported to represent only about 1~10% of overall detectable peaks in scATAC-seq analysis. Therefore, instead of using the initial cell-to-feature matrices directly for downstream analysis, data transformation can be applied to compensate for the limitation from data sparsity.
- Term-frequency inverse-document-frequency (TF-IDF): transforms a cell-to-feature matrix to give more weight to rarer peaks in the cell population.
- Other methods: Jaccard distance
### 6. Dimension reduction, visualization and clustering
- Purpose: To overcome inherit sparsity, the cell-by-feature matrix undergoes DR which can mitigate redundant information and potential noise of high dimensional data, and may reduce the computational time for downstream analysis.
- PCA and number of components selected by the elbow of scree plot analysis or Jackstraw test.
- Other DR methods: topic modeling (cisTopic) generated by latent Dirichlet allocation (LDA), Latent Semantic indexing (LSI), Multidimensional scaling and diffusion map.
- Visualization: t-SNE and/or UMAP
- Clustering: Cells with similar accessibility profiles can be organized into clusters. hierarchical, k-means, k-medoids, and Louvain algorithm.
## Reviews, methods and other interesting reads
- Accelerating single-cell genomic analysis with GPUs. [Nolet et al, 2022](https://www.biorxiv.org/content/10.1101/2022.05.26.493607v1) [Github](https://github.com/NVIDIA-Genomics-Research/rapids-single-cell-examples)
- An open-source interactive pipeline tutorial for differential ATAC-seq footprint analysis - Google cloud. [Github](https://github.com/NIGMS/ATAC-Seq-and-Single-Cell-ATAC-Seq-Analysis/tree/main)
- Single-cell ATAC sequencing analysis: From data preprocessing to hypothesis generation. [Baek and Lee, 2020](https://www.sciencedirect.com/science/article/pii/S2001037020303019)
