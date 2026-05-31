# Immune-Inflammatory Transcriptomic Dysregulation in Dilated Cardiomyopathy

## Overview
A multi-layer computational analysis of immune and inflammatory gene dysregulation in Dilated Cardiomyopathy (DCM) using publicly available RNA-seq data. This project integrates differential expression, co-expression network analysis, miRNA regulatory network construction, and drug repurposing to identify therapeutic targets.

## Author
**Marjia Islam Tia**
BSc Bioinformatics and Biotechnology, Asian University for Women, Bangladesh

## Datasets
| Accession | Description | Role |
|-----------|-------------|------|
| GSE116250 | RNA-seq, human left ventricle, DCM vs non-failing (n=51) | Primary discovery |
| GSE141910 | RNA-seq, human cardiac tissue, DCM vs non-failing (n=332) | Validation |

## Pipeline
Raw GEO Data → Preprocessing → Differential Expression → Functional Enrichment → WGCNA → miRNA Network → Drug Repurposing → Validation

## Methods Summary

### 1. Data Preprocessing (Script 01)
- Downloaded RPKM expression data from GSE116250
- Filtered low-expression genes (RPKM > 1 in ≥50% samples)
- Normalized using voom transformation
- Quality control via PCA — confirmed disease-driven transcriptomic separation

### 2. Differential Expression Analysis (Script 02)
- Tool: limma-voom (appropriate for RPKM-normalized data)
- Threshold: adj.P.Val < 0.05, |log2FC| > 1.5
- Result: 200 significant DEGs (129 upregulated, 71 downregulated in DCM)
- Key genes: NAMPT, TRBV5-4, NQO1, SERPINA3

### 3. Functional Enrichment Analysis (Script 03)
- GO Biological Process, Molecular Function, Cellular Component
- Top processes: striated muscle development, ECM organization, TGF-β signaling, oxidative stress response

### 4. WGCNA Co-expression Network (Script 04)
- Soft thresholding power: 9 (R² > 0.80)
- 40 modules identified
- Blue module (n=1,733 genes): strongest DCM correlation (r=0.86, p=3.7×10⁻¹⁶)
- Hub genes: PCBP1, PPP2R1A, TKFC, TAOK2, LCAT, TRIM41

### 5. Regulatory Network Analysis (Script 05)
- 2,408 validated miRNA-mRNA interactions via multiMiR
- Key miRNAs: hsa-miR-34a-5p, hsa-let-7 family, hsa-miR-16-5p
- Drug repurposing: TAOK2→AMG28, LCAT→MEDI6012

### 6. Cross-Validation (Script 06)
- Validated in GSE141910 (n=332 patients)
- 11/17 hub genes replicated with consistent upregulation
- Immune pathway conservation confirmed independently

## Key Findings
- DCM shows consistent immune infiltration and inflammatory activation across two independent cohorts
- PPP2R1A and TAOK2 identified as hub inflammatory kinase regulators
- LCAT identified as druggable hub; MEDI6012 as repurposing candidate
- miR-34a-5p and let-7 family as upstream regulators of DCM hub gene network

## Repository Structure
```
├── scripts/          # R analysis scripts (00-06)
|
├── results/
│   ├── figures/      # All generated plots
│   └── tables/       # DEG lists, enrichment results, hub genes
|
├── data/
│   └── processed/    # Normalized expression objects
```
## Requirements
Run scripts/00_install_packages.R before any other script.
R version: 4.4.3

## Data Availability
- GSE116250: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116250
- GSE141910: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141910

## License
MIT License
