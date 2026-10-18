# Otolith shape analysis of Atlantic herring populations in the Western Baltic Sea

R code for the Master's thesis by Laura Charlotte Thiel (University of Rostock, 2026),
written in cooperation with the Thünen Institute of Baltic Sea Fisheries.

Herring otoliths from the German Acoustic Autumn Survey (GERAS, 2023–2025, ICES
Subdivisions 21–24) were genetically assigned to five populations (WBSS, BAS, NSAS,
Downs, CBH). Their outlines were described by elliptic Fourier analysis
(12 harmonics), and the populations were separated by linear discriminant analysis
with repeated 4-fold cross-validation. The approach follows Smoliński, Schade & Berg
(2020, *Can. J. Fish. Aquat. Sci.* 77: 674–683).

**This repository contains code only.** The otolith photographs and the metadata
(genetic assignment, age) belong to the Thünen Institute and are not included.

## Requirements

R ≥ 4.6.1 with the packages

```r
install.packages(c("tidyverse", "readxl", "Momocs", "caret", "MASS",
                   "doParallel", "foreach", "jpeg", "imager", "magick",
                   "ragg", "systemfonts", "scales"))
```

`imager` and `magick` need ImageMagick and FFTW (macOS: `brew install imagemagick fftw`).

## Paths

All paths are defined in [`config.R`](config.R); every script loads it. The project
folder is found automatically, or can be set in `~/.Renviron`:

```
SHAPE_ANALYSIS_DIR=/path/to/Shape Analysis
```

`source("config.R"); check_paths()` checks that all folders are found.
All analysis results are written to `Fourier and Classification/Results/`.

## Run order

| Step | Script | What it does |
|---|---|---|
| 1 | `1 outline extraction/ex_oto.R` | Binarises the photographs and extracts the otolith outlines |
| 2 | `2 metadata/age_join.R` | Adds age to the metadata |
| 3 | `3 analysis/main_analysis.R` | Joins outlines and metadata, removes outliers (Mahalanobis D² on PC1–5), elliptic Fourier analysis, LDA with repeated 4-fold cross-validation for all analysis approaches (age groups, one-vs-rest, class imbalance) |
| 4 | `3 analysis/additional_analyses.R` | Spawning-group model |
| 5 | `3 analysis/harmonics_sweep.R` | Classification accuracy and kappa over 2–35 harmonics |

Step 3 must run before steps 4–5 and before all figure scripts.

## Figures and tables of the thesis

| Thesis | Script (`4 figures and tables/`) |
|---|---|
| Figs. 10, 11 – sample composition; Figs. 16, 17 – mean outlines | `fig_composition_and_mean_shapes.R` |
| Fig. 12 – number of harmonics | `fig_number_of_harmonics.R` |
| Fig. 13 – harmonic series of one otolith | `fig_harmonic_series.R` (uses `1 outline extraction/harmonic_series_single_otolith.R`) |
| Fig. 14 – outlier removal; χ² test of the removal rates (Section 3.3) | `fig_outliers.R` |
| Fig. 15 – principal component analysis | `fig_pca.R` |
| PCA numbers in Sections 3.4 and 3.5 (variance explained, group means) | `pca_numbers.R` |
| Table 2 – population composition; retention per population (Section 3.1) | `stock_composition.R` |
| Tables 3, 4 – class imbalance, classification performance | `results_tables.R` (needs `stock_composition.R` first) |
| Table 5 – yes/no accuracy for WBSS and BAS | `table_yes_no.R` |
| Tables A1, A2 – age composition by population and by subdivision | `appendix_tables_age.R` |
| Fig. A1 – examples of otoliths removed as outliers | `appendix_outlier_examples.R` |
| Fig. A2 – PCAs of the age subsets | `appendix_pca_age0_3_and_age1_10.R` |
| Fig. A3 – mean outlines of the age subsets | `appendix_mean_shapes_age0_3_and_age1_10.R` |
| Figs. A4–A6 – confusion matrices | `appendix_confusion_matrices.R` |

The remaining figures of the thesis are schematic illustrations, photographs, maps or
figures from the literature and were not produced with this code.


