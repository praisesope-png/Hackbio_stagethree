# **COVID-19 scRNA‑seq Analysis Pipeline**

**Reference:** PLOS Biology (2021) - Trajectory analysis of SARS-CoV-2-infected human bronchial epithelial cells.

**Data Source:  
**GSE166766 (MTX + TSV files)

**Objective:  
**Reproduce the neighbourhood clustering, cell type identification, and infection-related dynamics described in the paper (Figures 1G(i-iii), 3A, 3B, 4A, 4B). Perform pseudotime analysis to order SARS-CoV-2 infection and differentiation states.

## **Project Overview**

This project analyzes airway epithelial cells across mock and SARS‑CoV‑2-infected conditions at different time points (0, 1, 2, and 3 days post‑infection). The analysis includes:

- Preprocessing and QC filtering  

- Normalization and HVG selection  

- PCA, neighborhood graph, UMAP  

- Leiden clustering  

- Marker‑based cell type annotation  

- Pseudotime reconstruction  

- Interpretation of viral infection dynamics  

## **1\. Cell Types Identified at Different Stages of Infection**

The dataset consistently contained the following epithelial populations:

- Ciliated cells
- Gastric chief cells
- Ionocytes
- Goblet cells
- Enteric neurons
- Enteroendocrine cells
- Tuft cells
- Clara cells  

### **2\. Why do these cell types correlate with COVID-19 infection?**

SARS-CoV-2 primarily targets airway epithelial cells that express **ACE2** and **TMPRSS2**. In this dataset:

- **Ciliated cells** are the main targets, showing high ACE2/TMPRSS2 expression and sustaining viral replication.
- **Goblet and Clara cells** are moderately permissive and contribute to mucus secretion and epithelial repair.
- **Ionocytes and Tuft cells** may not be heavily infected but can play supportive roles in immune signaling and airway homeostasis.
- Rare or non-airway cells (gastric chief, enteric neurons, enteroendocrine cells) likely reflect transitional states or annotation overlaps and are not primary infection sites.

Overall, infection susceptibility aligns with ACE2/TMPRSS2 expression and the functional role of each epithelial cell type, explaining why ciliated and secretory cells dominate viral tropism while other cells respond indirectly.  

## **3\. Is ACE2 a Good Marker for Tracking Infection in This Dataset?**

Not really.

In this dataset, viral RNA and ACE2 levels do not correlate strongly.

- ACE2 is expressed in only a small subset of epithelial cells.
- Many ACE2-positive cells do not contain viral reads.
- Many infected cells have low or undetectable ACE2, likely due to: viral-induced ACE2 downregulation, and/or infection occurring via alternate or transient expression states.

In conclusion: ACE2 is a good susceptibility marker, not a reliable infection-progress marker.

## **4\. Difference Between ENO2 and ACE2 as Biomarkers (Comparing Two Studies)**

**ACE2**

- Entry receptor required for viral uptake.
- Sparse expression, high cell-type specificity.
- Not induced during infection; often downregulated.

**ENO2**

- A stress- and infection-induced metabolic marker.
- Strongly upregulated in infected cells.
- Correlates with viral RNA abundance.
- Better at distinguishing **actively infected** versus **bystander** cells.

**Summary:  
**_ACE2 = susceptibility marker_

_ENO2 = infection/progression marker_  

### **5\. Which cell cluster has the highest abundance of ACE2 expression after 3 dpi, and what does that mean biologically?**

Based on the reprocessed clustering:

- The Ciliated cell cluster shows the highest ACE2 abundance at 3 dpi.

Biological Interpretation:

- Even though ACE2 expression may decrease upon infection, the remaining detectable ACE2 is still concentrated in ciliated epithelial cells.
- This supports the model that ciliated cells remain the primary viral reservoir throughout the infection course.
- At 3 dpi, high ACE2 in these cells corresponds to:
  - higher viral RNA loads
  - strong ISG responses
  - epithelial damage and dedifferentiation

In visual UMAPs, this appears as **ACE2-positive "hot spots"** localised within the ciliated cluster.  

**Repository Structure**

├── data/

│ └── GSE166766_raw/

│ ├── 01_preprocessing.ipynb

│ ├── 02_clustering_annotation.ipynb

│ ├── 03_trajectory_pseudotime.ipynb

│ └── 04_biomarker_analysis.ipynb

├── README.md

- Citation

Ravindra, N. G., Alfajaro, M. M., Gasque, V., Huston, N. C., Wan, H., Szigeti-Buck, K., Yasumoto, Y., Greaney, A. M., Habet, V., Chow, R. D., Chen, J. S., Wei, J., Filler, R. B., Wang, B., Wang, G., Niklason, L. E., Montgomery, R. R., Eisenbarth, S. C., Chen, S., Williams, A., … Wilen, C. B. (2021). Single-cell longitudinal analysis of SARS-CoV-2 infection in human airway epithelium identifies target cells, alterations in gene expression, and cell state changes. _PLoS biology_, _19_(3), e3001143. <https://doi.org/10.1371/journal.pbio.3001143>
