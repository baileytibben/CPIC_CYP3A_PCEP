# CYP3A4 In Vitro Activity Analysis

This repository contains R scripts and supporting files used to analyze in vitro functional data for CYP3A4 star alleles across multiple substrates and experimental systems.

## 📌 Project Overview

The goal of this analysis is to evaluate allele-specific activity of CYP3A4 variants using curated in vitro data and to assess variability across substrates and experimental systems. These results support clinical allele function assignment.

---

## 📊 Dataset Description

The dataset was manually curated from published literature reporting in vitro functional characterization of CYP3A4 alleles.

### Data sources
- Peer-reviewed publications reporting enzyme activity for CYP3A4 variants
- Activity values were extracted from tables within manuscripts

### Activity values
- Reported as **percent activity relative to reference (wild-type) CYP3A4*1
- In most cases, values were directly extracted as reported in the publication
- In cases where percent activity was not explicitly reported, it was **calculated from kinetic parameters** using:
  Percent activity = (Vmax/Km_variant) / (Vmax/Km_reference) × 100

### Additional annotations
Each observation may include:
- Substrate
- Experimental system (e.g., mammalian cells, bacterial expression, insect cells)
- PMID / publication source

###  Data Processing Details

- When multiple measurements were reported for the same allele–substrate pair, values were retained separately
- Missing values were handled as NA and excluded from statistical modeling where appropriate

---
## ▶️ Reproducibility

### Requirements
- R (version 4.5.1)
- Required packages listed in `sessionInfo.txt`

### To reproduce analysis:
1. Run `scripts/01_cleaning.R`
2. Run `scripts/02_analysis.R`
3. Run `scripts/03_figures.R`

---

## ⚠️ Notes on Data

- Data were curated manually and may reflect variability in experimental design across studies
- Differences in substrates and expression systems may influence measured activity
- Calculated activity values (from kinetic parameters) assume comparable assay conditions between variant and reference enzymes within each study

---

