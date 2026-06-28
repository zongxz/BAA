# Brain Age Acceleration (BAA) - Multi-omic Analysis

## Repository Overview
This repository contains all custom scripts and pipelines for the data analysis presented in the manuscript:

**"Integrative Analysis of Brain Age Acceleration: Health Associations and Multi-Omics Profiles"**

The code encompasses the entire analysis pipeline, from brain age prediction using UK Biobank imaging data to downstream genome-wide association studies (GWAS), mediation analyses, and multi-omic integrative analyses.

## Study Background & Key Findings
Brain-Age Acceleration (BAA) reflects the discrepancy between an individual's predicted brain age (based on neuroimaging) and their chronological age. Its biological underpinnings and health implications are not fully understood.

**Key findings from this study include:**
*   **Health Impacts**: BAA is significantly associated with mortality and diseases affecting multiple organ systems.
*   **Lifestyle Mediation**: BAA mediates the effects of lifestyle factors on health outcomes.
*   **Genetic Architecture**: A GWAS identified multiple novel genetic risk loci for BAA. A polygenic risk score (PRS) for BAA was associated with multi-system diseases.
*   **Cellular & Molecular Mechanisms**: Brain single-cell transcriptomics implicate endothelial cells, pericytes, and astrocytes. Integrative analyses (SMR, TWAS) identified genes like *BTN3A2* with robust pan-brain correlations to BAA.
*   **Circulating Biomarkers**: Proteomic and metabolomic profiling revealed a broad spectrum of biomarkers that mediate the effects of modifiable factors and genetic risk on brain aging.

## Data Source
This study analyzed data from **over 60,000 individuals** from the **UK Biobank**.
*   **Brain MRI**: Used to predict brain age.
*   **Genomic Data**: Used for GWAS and PRS calculation.
*   **Proteomic & Metabolomic Data**: Used for circulating biomarker profiling.
*   **Single-cell Transcriptomic Data**: Used for cell-type specificity analysis.

*Note: Due to data access restrictions, raw data cannot be shared in this repository. Access to UK Biobank data is managed through the UK Biobank Access Management System.*

## Repository Structure & Scripts
The scripts are organized to reproduce the major steps of the analysis pipeline:
