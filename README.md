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
### 4. Batch correction and data integration
### 5. Data transformation
### 6. Dimension reduction, visualization and clustering
## Downstream analysis for hypothesis generation

## Reviews, methods and other interesting reads
- Single-cell ATAC sequencing analysis: From data preprocessing to hypothesis generation. [Baek and Lee, 2020](https://www.sciencedirect.com/science/article/pii/S2001037020303019)
