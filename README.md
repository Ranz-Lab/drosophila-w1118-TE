# drosophila-w1118-TE

This repository contains all code used in the analysis of transposable element (TE) expression from single-nucleus RNA-seq data of *Drosophila melanogaster* *w<sup>1118<sup>* testis. The associated study is entitled “Heterochromatic loci drive the transposable element expression burst during *Drosophila* spermatogenesis.”

## 📁 Repository Structure

| Level | Folder                       | What it contains                                                                                                     |
| ----- | ---------------------------- | -------------------------------------------------------------------------------------------------------------------- |
| **0** | **`00_TE_quantification/`**  | Quantification and processing of TE expression at the locus and family levels                                        |
|       | ├── `00_SoloTE/`             | Locus-level TE quantification using SoloTE                                                                           |
|       | └── `01_scTE/`               | TE-family-level quantification using scTE                                                                            |
| **1** | **`01_genomic_context/`**    | Analysis of TE expression across genomic compartments and distances between expressed TE loci and neighbouring genes |                         
| **3** | **`02_H3K9me2/`**            | Integration of TE expression with H3K9me2 profiles across genomic compartments                                       |
| **4** | **`03_small_RNA_pathways/`** | Relationships between TE expression and piRNA- and RNAi-pathway gene expression                                      |
| **5** | **`04_hdWGCNA/`**            | TE co-expression network analysis and identification of cell-type-associated TE modules      |

---

## 📦 Data Access

Single-nucleus data was obtained from the Fly Cell Atlas. The *w<sup>1118<sup>* genome assembly and raw Nanopore and Illumina sequencing reads will be made available through the associated NCBI BioProject and data repositories. Accession information will be added upon public release.
