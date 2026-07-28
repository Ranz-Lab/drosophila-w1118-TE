# drosophila-w1118-TE

This repository contains all code used in the analysis of transposable element (TE) expression from single-nucleus RNA-seq data of *Drosophila melanogaster* **w1118** testis. The associated study is entitled “Heterochromatic loci drive the transposable element expression burst during *Drosophila* spermatogenesis.”

## 📁 Repository Structure

| Level | Folder                       | What it contains                                                                                                     |
| ----- | ---------------------------- | -------------------------------------------------------------------------------------------------------------------- |
| **0** | **`00_TE_quantification/`**  | Quantification and processing of TE expression at the locus and family levels                                        |
|       | ├── `00_SoloTE/`             | Locus-level TE quantification using SoloTE                                                                           |
|       | └── `01_scTE/`               | TE-family-level quantification using scTE                                                                            |
| **1** | **`01_TE_expression/`**      | Processing, integration, pseudobulk aggregation, and comparison of TE expression across testis cell types            |                        |
| **2** | **`02_genomic_context/`**    | Analysis of TE expression across genomic compartments and distances between expressed TE loci and neighbouring genes |
| **3** | **`03_permutation_tests/`**  | Permutation-based tests of associations between TE loci and neighbouring genes                                       |
| **4** | **`04_H3K9me2/`**            | Integration of TE expression with H3K9me2 profiles across genomic compartments                                       |
| **5** | **`05_small_RNA_pathways/`** | Relationships between TE expression and piRNA- and RNAi-pathway gene expression                                      |
| **6** | **`06_hdWGCNA/`**            | TE co-expression network analysis and identification of cell-type-associated TE modules      |

---

## 📦 Data Access

The analyses use single-nucleus RNA-seq data from *D. melanogaster* **w1118** testes together with a strain-matched de novo genome assembly and manually curated TE annotation.

Raw sequencing data, the genome assembly, TE annotation, and processed expression matrices will be made available through the associated NCBI BioProject and data repositories. Accession information will be added upon public release.
