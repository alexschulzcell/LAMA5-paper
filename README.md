# LAMA5-dependent organization of chondrogenic transcriptional programs

This repository contains all R scripts used for transcriptomic and qPCR analyses associated with the LAMA5 manuscript.

Analyses include quality control, differential gene expression across multiple experimental contexts, weighted gene co-expression network analysis (WGCNA), reference dataset analyses, and qPCR data processing.  
Most analyses are independent and can be run separately.

---

## Repository organization

The repository is organized into the following main folders:

- `QC/` – quality control analyses  
- `differential_expression/` – differential gene expression analyses  
- `WGCNA/` – weighted gene co-expression network analysis  
- `analysis_of_references/` – analyses using external reference datasets  
- `qPCR/` – quantitative PCR analyses
- `Imaging/` – Macros used for spheroid size measurements

---

## QC

| Folder / Scripts | Description |
|---|---|
| `QC/*.R` | Quality control of RNA-seq datasets, including filtering, exploratory analyses, and sample-level metrics. |

---

## Differential expression

Differential gene expression analyses were performed for multiple experimental contexts using DESeq2-based workflows.

### **all_celltypes**

| Folder / Scripts | Description |
|---|---|
| `differential_expression/all_celltypes/*.R` | Differential expression analysis across all profiled cell types (Undifferentiated, chondrogenic, osteogenic). |

---

### **chondrogenic_USCs**

| Folder / Scripts | Description |
|---|---|
| `differential_expression/chondrogenic_USCs/*.R` | Differential expression analysis of chondrogenically induced urine-derived stem cells (USCs). |

---

### **undifferentiated_USCs**

| Folder / Scripts | Description |
|---|---|
| `differential_expression/undifferentiated_USCs/*.R` | Differential expression analysis of undifferentiated USCs. |

---

### **differentiation**

| Folder / Scripts | Description |
|---|---|
| `differential_expression/differentiation/*.R` | Differential expression analyses across differentiation conditions. |

---

### **LiCl_rescue**

| Folder / Scripts | Description |
|---|---|
| `differential_expression/LiCl_rescue/*.R` | Differential expression analysis of LiCl rescue experiments assessing modulation of WNT signaling. |

---

## WGCNA

| Folder / Scripts | Description |
|---|---|
| `WGCNA/*.R` | Weighted gene co-expression network analysis to identify gene modules associated with experimental conditions and phenotypes. |

---

## Analysis of references

| Folder / Script | Description |
|---|---|
| `analysis_of_references/*.R` | Analyses using external reference datasets for contextualization and comparison with the study data. |

---

## qPCR

| Folder / Scripts | Description |
|---|---|
| `qPCR/*.R` | qPCR data processing, normalization, and statistical analysis. |

---
## qPCR

| Folder / Scripts | Description |
|---|---|
| `Imaging/*.ijm` | Size measurements of chondrogenic spheroids. |

---

## Repository tree
<pre>
.
├── QC/
├── differential_expression/
│   ├── all_celltypes/
│   ├── chondrogenic_USCs/
│   ├── differentiation/
│   ├── LiCl_rescue/
│   └── undifferentiated_USCs/
├── WGCNA/
├── analysis_of_references/
├── qPCR/
└── Imaging/
</pre>


---

## Required packages

Analyses were performed in R (≥4.2).  
Key packages include (not exhaustive):

| Package / Tool | Purpose |
|---|---|
| DESeq2 | Differential gene expression |
| WGCNA | Co-expression network analysis |
| ggplot2 | Data visualization |
| ggrepel | Plot labeling |
| clusterProfiler | Functional enrichment |
| org.Hs.eg.db | Gene annotation |

Exact package versions are provided via `sessionInfo()` outputs or `renv.lock` (if applicable).

---

## Data availability

The RNA-seq raw data are available at ArrayExpress: accession E-MTAB-16566.

---

## Code availability

All analysis scripts are provided to support transparency and reproducibility of the analyses reported in the manuscript.  
Scripts are organized by analysis type and experimental context; most analyses can be run independently.
