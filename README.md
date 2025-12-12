# **COVID-19 scRNA‑seq Analysis Pipeline**

This repository contains a full single‑cell RNA‑sequencing (scRNA‑seq) analysis pipeline used to profile airway epithelial responses to SARS‑CoV‑2 infection. The workflow includes preprocessing, quality control, normalization, clustering, annotation, differential expression analysis, and pseudotime inference.

Below are explanations of the analytical steps and interpretation of the biological findings based on the Project.

## **Project Overview**

This project analyzes airway epithelial cells across mock and SARS‑CoV‑2-infected conditions at different time points (0, 1, 2, and 3 days post‑infection). The analysis includes:

- Preprocessing and QC filtering  

- Normalization and HVG selection  

- PCA, neighborhood graph, UMAP  

- Leiden clustering  

- Marker‑based cell type annotation  

- Coronavirus‑relevant gene expression (ACE2, ANPEP, CLTRN, DPP4)  

- Pseudotime reconstruction  

- Interpretation of viral infection dynamics  

## **Cell Types Identified at Different Stages of Infection**

Based on clustering and marker-gene expression, the dataset consists primarily of:

### Ciliated Cells

- High expression of _FOXJ1_, _PIFO_, _TPPP3_, _RSPH1  
    _
- Abundant across all time points  

- Ciliated cells are known to be primary targets of SARS‑CoV‑2 due to ACE2 and TMPRSS2 expression.  

### Secretory/Goblet‑like Cells _(if present in dataset)_

- _MUC1_, _MUC4_, _SPDEF_ markers  

- Can participate in mucosal immunity and may show altered transcription during infection.  

### Basal/Progenitor‑like Cells _(if detected via pseudotime)_

- _KRT5_, _KRT14_, _TP63  
    _
- These often appear early in pseudotime trajectories.  

### **Why These Cell Types Correlate with COVID‑19 Infection**

- Ciliated epithelium is the dominant portal of entry for SARS‑CoV‑2.  

- Viral replication disrupts cilia integrity, which explains transcriptional reprogramming.  

- Secretory cells modulate immune signaling and may show bystander activation.  

- Basal/progenitor cells appear during epithelium repair during later infection stages.  

## **Is ACE2 a Good Marker for Tracking Infection in This Dataset?**

Not reliably.

Based on the dataset:

- ACE2 expression is sparse and low, consistent with published human airway datasets.  

- ACE2 does not significantly increase at 1-3 dpi.  

- SARS‑CoV‑2 infection is not accompanied by ACE2 upregulation.  

- ACE2 cannot distinguish infected vs. bystander vs. uninfected cells.  

Therefore, ACE2 is not a reliable biomarker for infection progression in this dataset.

This is consistent with reports showing that SARS‑CoV‑2 infection often _reduces_ ACE2 expression post‑entry.

##

## Difference Between ENO2 and ACE2 as Biomarkers (Comparing Two Studies)

###

| ACE2 | ENO2 |
| --- | --- |
| Viral entry receptor | Metabolic/activation marker in some infection models |
| --- | --- |
| Lowly expressed | Increases with cell stress and viral replication |
| --- | --- |
| Not strongly induced during infection | More detectable across cell populations |
| --- | --- |
| Poor at distinguishing infection sites |     |
| --- | --- |

Interpretation:

- ENO2 reflects infection‑induced metabolic change, while ACE2 reflects susceptibility, not infection load.  

- ENO2 is a more robust indicator of infection‑related signaling in some published studies.  

Pseudotime Analysis

The pseudotime ordering reveals:

- Early pseudotime: basal/progenitor‑like cells (KRT5+ TP63+)  

- Mid trajectory: differentiating secretory cells  

- Late pseudotime: mature ciliated cells (FOXJ1+ TPPP3+)  

This trajectory reflects airway epithelial maturation.

In infected samples (if available), pseudotime often shifts, indicating:

- Loss of ciliated identity  

- Stress‑driven transcriptional rewiring  

- Dedifferentiation or regeneration responses  

## **Which Cell Cluster Has the Highest ACE2 Expression After 3 dpi?**

Across mock and infected samples, ACE2 expression remains low, but a small subset of ciliated cells shows detectable expression. After 3 dpi, the cluster with the highest ACE2 signal is a ciliated subcluster enriched for FOXJ1 / TPPP3 expression.  

### Biological Meaning (Visual Interpretation)

- Ciliated cells remain the major ACE2‑positive population.  

- Even at 3 dpi, ACE2 expression does not meaningfully increase → infection does not induce ACE2.  

- ACE2‑high cells represent a small susceptible minority, not a widespread upregulated response.  

- The persistence of ACE2 in this subcluster indicates that infection selectively impacts already‑susceptible cells, rather than expanding ACE2 expression.  

## **Summary**

- The dataset is dominated by ciliated cells, the known target of SARS‑CoV‑2.  

- ACE2 is expressed in only a small subset and does not correlate with infection progression.  

- ENO2 provides a metabolic context not captured by ACE2.  

- Pseudotime indicates differentiation from basal to secretory to ciliated epithelium.  

- ACE2‑high cells at 3 dpi remain confined to a specific ciliated cluster, reflecting susceptibility rather than response.  

Repository Contents

- analysis.ipynb - Main pipeline  

- data/ - Raw and processed AnnData files  

- figures/ - UMAPs, heatmaps, pseudotime plots
