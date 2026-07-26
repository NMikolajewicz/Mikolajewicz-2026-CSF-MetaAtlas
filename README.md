# CSF proteomic meta-atlas

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21608891.svg)](https://doi.org/10.5281/zenodo.21608891)

This repository contains the processed data and R code used to reproduce the
main analyses reported in **“A shared injury–homeostasis gradient in the human
CSF proteome across neurological disease.”** The study combines 35
cerebrospinal-fluid (CSF) proteomic datasets, defines a 747-protein
mass-spectrometry-accessible core CSF proteome, identifies 13 recurrent
proteomic meta-programs, and evaluates C5 and NRCAM as a two-protein estimate
of a shared injury–homeostasis axis.

## Repository contents

```text
data/       Processed input tables used by the numbered scripts
scripts/    R scripts for analyses 01–09
analysis/   R Markdown analysis files; generated tables and figures are written here
README.md   Analysis description and instructions
CITATION.cff  Citation metadata for archived releases
```

All paths are relative to the repository root. The R scripts use base R only;
no package installation or workflow software is required.

## Software environment

The primary analysis requires only base R. Release 1.0.0 was tested from a
clean clone with R 4.5.2 on macOS 15.6.1 (arm64) and was also checked with R
4.6.1. The complete run does not install packages or access the network.

The optional HTML rendering step was tested with `rmarkdown` 2.31, `knitr`
1.51 and Pandoc 3.8.3. These packages are not used by the numbered R scripts.

## Running the analyses

To reproduce the archived release, clone the tagged source and run the complete
set from the repository root:

```sh
git clone --branch v1.0.0 --depth 1 https://github.com/NMikolajewicz/Mikolajewicz-2026-CSF-MetaAtlas.git
cd Mikolajewicz-2026-CSF-MetaAtlas
Rscript scripts/run_all.R
```

The final command runs analyses 01–09 and checks the principal published
values. A successful run ends with `All analyses completed and verified.` The
verification can also be repeated separately with:

```sh
Rscript scripts/verify_results.R
```

The scripts can also be run individually:

```sh
Rscript scripts/01_matrix_correction.R
Rscript scripts/02_assayable_universe.R
Rscript scripts/03_core_proteome.R
Rscript scripts/04_detection_model.R
Rscript scripts/05_core_classes.R
Rscript scripts/06_meta_programs.R
Rscript scripts/07_injury_homeostasis_axis.R
Rscript scripts/08_enrichment.R
Rscript scripts/09_axis_prediction.R
```

Results are written to a matching subdirectory under `analysis/`. Existing
files in those output directories are replaced when a script is rerun. The
complete processed-data analysis normally finishes in a few minutes on a
standard laptop.

The R Markdown files provide short, executable descriptions of each analysis.
With the optional rendering software listed above, all nine can be rendered
with:

```sh
Rscript -e 'for (f in list.files("analysis", pattern = "[.]Rmd$", full.names = TRUE)) rmarkdown::render(f, quiet = TRUE)'
```

## Analyses and manuscript figures

| Script | Analysis | Manuscript panels |
|---|---|---|
| `01_matrix_correction.R` | Ten-state matrix-correction benchmark and composite method ranking | Extended Data Fig. 4 |
| `02_assayable_universe.R` | Rolling three-study local-support estimate of assayability in 21 MS studies | Extended Data Fig. 1a,b |
| `03_core_proteome.R` | Coverage-aware detection meta-analysis, recurrence tiers and core-proteome list | Fig. 2a–e; Extended Data Fig. 1c–f |
| `04_detection_model.R` | Assessment of the source-specific elastic-net detection model | Fig. 2f–i |
| `05_core_classes.R` | Gower ordination and characterization of six Core Classes | Fig. 2j,k; Extended Data Fig. 2 |
| `06_meta_programs.R` | NMF/ICA/SSN filtering funnel, consensus meta-programs and held-out-study transfer | Fig. 3; Extended Data Fig. 5 and 6 |
| `07_injury_homeostasis_axis.R` | Meta-program clustering, PC1 loadings and protein-axis associations | Fig. 4a–c,f; Extended Data Fig. 8–10 |
| `08_enrichment.R` | Core Class and meta-program over-representation results and axis-ranked enrichment | Fig. 2k, Fig. 3e and Fig. 4d,e |
| `09_axis_prediction.R` | Bader, Riviere-Cazaux and all-study C5/NRCAM prediction analyses | Fig. 4m–o |

The scripts generate individual analytical panels and their source tables.
Assembly of the final multi-panel figures, panel lettering and conceptual
schematics was done separately for the manuscript.

## Input data

The files in `data/` are compressed tab-separated tables where compression is
useful. R reads these files directly; they do not need to be unpacked.

| Directory | Contents |
|---|---|
| `01_matrix_correction/` | Nine batch-removal and biology-preservation metrics for ten matrix states on the 917-protein benchmark panel |
| `02_assayable_universe/` | Protein-by-study assayability flags and the depth-ordered rolling windows for 21 studies |
| `03_core_proteome/` | Random-effects detection estimates, confidence intervals, heterogeneity and leave-one-study-out summaries for 12,606 proteins |
| `04_detection_model/` | Cross-validated predictions for 2,762 proteins and the fitted source-head coefficients |
| `05_core_classes/` | Forty mixed protein features and the six manuscript Core Class assignments for 747 proteins |
| `06_meta_programs/` | Component filtering counts, consensus memberships, 524-protein catalogue, leave-one-study-out results and variance summaries |
| `07_injury_homeostasis_axis/` | Meta-program correlations, PC1 loadings, permutation result and protein-axis associations |
| `08_enrichment/` | Core Class, meta-program and axis enrichment statistics and the meta-program pathway network |
| `09_axis_prediction/` | Aggregate model ladders, study-level C5+NRCAM performance and random-effects meta-analysis inputs |

All included data are aggregate or protein-level. There are no participant
identifiers, sample identifiers, diagnoses, clinical dates or individual-level
measurements in the repository. Some protein annotations and pathway names were
derived from third-party resources and remain subject to the terms of their
original providers.

## Reproducibility scope

These scripts reproduce the reported summary statistics, filtering results,
protein lists, meta-analyses and analytical figure panels from processed inputs.
They do not repeat vendor-file processing, peptide identification, the original
sample-level NMF/ICA/SSN decompositions, GSVA scoring, or participant-level model
fitting. Those stages require study data that cannot be distributed openly in
this repository. The distinction is also noted in the relevant R Markdown file.

The main checks reproduced by `scripts/run_all.R` include:

- a limma correction composite score of 0.674 across ten benchmark states;
- core-proteome tier counts of 337, 410, 497, 1,518 and 9,844 proteins;
- detection-model Spearman correlation of 0.540, AUROC of 0.833 and AUPRC of 0.691;
- six Core Classes containing 190, 185, 92, 193, 42 and 45 proteins;
- the 28,795 → 22,271 → 10,805 → 634 → 247 → 13 program-filtering series and 524 unique meta-program proteins;
- an axis PC1 explaining 55.6% of meta-program variation, with 71 injury-associated and 236 homeostasis-associated proteins; and
- all-study random-effects C5+NRCAM R² of 0.812.

## Citation and contact

The archived software release is available at
[https://doi.org/10.5281/zenodo.21608891](https://doi.org/10.5281/zenodo.21608891).

Nicholas Mikolajewicz, Cecile Riviere-Cazaux, Kyle Tuohy, Sruthi Ranganathan,
Rahul Kumar, Alexandra M. Miller, Manmeet S. Ahluwalia, Paul C. Boutros,
Chetan Bettegowda, Terry Burns, Thomas Kislinger and Alireza Mansouri. **A
shared injury–homeostasis gradient in the human CSF proteome across neurological
disease.** 2026.

Please update the article citation with the journal reference and article DOI
after publication.
Questions about the analysis may be directed to Nicholas Mikolajewicz. Study
correspondence: Alireza Mansouri (`amansouri@pennstatehealth.psu.edu`).

## Licence

The R code is released under the MIT Licence. Data tables remain subject to the
terms of the original studies and annotation resources and are provided here
for reproduction of the analyses reported in the paper.
