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
### 5. Data transformation
### 6. Dimension reduction, visualization and clustering
## Downstream analysis for hypothesis generation

## Reviews, methods and other interesting reads
- Accelerating single-cell genomic analysis with GPUs. [Nolet et al, 2022](https://www.biorxiv.org/content/10.1101/2022.05.26.493607v1) [Github](https://github.com/NVIDIA-Genomics-Research/rapids-single-cell-examples)
- Single-cell ATAC sequencing analysis: From data preprocessing to hypothesis generation. [Baek and Lee, 2020](https://www.sciencedirect.com/science/article/pii/S2001037020303019)
