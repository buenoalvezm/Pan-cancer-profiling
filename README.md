# Pan-cancer-profiling
 
## Description

This code was generated in the context of the study “****Next generation pan-cancer blood proteome profiling using proximity extension assay”,**** were we performed a comprehensive analysis of the plasma proteome of a pan-cancer cohort representing the major cancer types. 


### The study

In this study, the plasma profiles of 1,463 proteins were measured for 1,477 cancer patients representing 12 common cancer types, including the most prevalent types such as colorectal-, breast-, lung- and prostate cancer. The plasma proteome was measured in minute amounts of blood plasma collected at the time of diagnosis and before treatment using the antibody-based PEA technology combined with next generation sequencing (developed by Olink). The plasma profiles of patients with a specific cancer type were compared to the patients with other cancer diagnosis in order to find cancer-specific signatures that can distinguish each type of cancer from other cancer types. Both differential expression and disease prediction models were used as tools for the identification of specific cancer signatures. 

The results from the study are available as a pre-print:
Mathias Uhlen, María Bueno Álvez, Fredrik Edfors et al. Next generation pan-cancer blood proteome profiling using proximity extension assay, 01 November 2022, PREPRINT (Version 1) available at Research Square [[https://doi.org/10.21203/rs.3.rs-2025767/v1](https://doi.org/10.21203/rs.3.rs-2025767/v1)].

### Content

This repository includes the code to generate the results describe above, as well as a synthetic dataset to test the code:

1. `/data`: contains example data to test the code, as well as additional data files to reproduce the analysis. Note that this is ****not real data**** and should not be used in any research. 
2. `/scripts`: contains all necessary scripts to reproduce the analysis.
3. `/results`: all the plots resulting from the analysis will be stored in this directory. Note that the results are based on synthetic data and should ****not**** be interpreted as ****valid biological resutlts****.
4. `Pan-cancer-profiling.Rproj`: R project file.

## Demo

### Before running
Before running the code, make sure you have R, R studio and the packages indicated in `data/processed/Sessions` installed.

The code has been developed in the following system:

- Processor Intel(R) Core(TM) i9-9880H CPU @ 2.30GHz
- Installed RAM 32,0 GB
- System type 64-bit operating system, x64-based processor

The provided run times apply for computers with similar specifications.

### Instructions to run data

1. Clone the repository (should take ~15 seconds).
2. Open `R studio` and open the `Pan-cancer-profiling.Rproj`.
3. Start by running through the differential expression analysis using the `Differential_expression.rmd` markdown script located in the scripts folder. This will perform differential expression analysis based on the example data located in the data folder.
4. Continue by running the `Disease_prediction.rmd` markdown script located in the scripts folder. This script will find protein signatures for the different groups of patients using prediction models based on the glmnet and random forest algorithms. It will also combine the results from the differential expression analysis to select a panel of upregulated proteins relevant for the prediction of the disease groups.
5. Explore the generated results.

The expected runtimes are:

- `Differential_expression.rmd`: ~3 minutes
- `Disease_prediction.rmd`: ~17 minutes

If you want to run the code using your data, make sure to format it according to the data provide here or adjust the script accordingly.  

### Expected output

All plots resulting from the analysis are stored in `results/YYYY-MM-DD` (a new folder will be created if you re-run the analysis on a different date).

All resulting data files (e.g. differential expression analysis, results from prediction models …) are stored in subfolders of the `data/processed` directory as R objects.

## Citation

To use this code in your own research, please cite our pre-print: 

Mathias Uhlen, María Bueno Álvez, Fredrik Edfors et al. Next generation pan-cancer blood proteome profiling using proximity extension assay, 01 November 2022, PREPRINT (Version 1) available at Research Square [[https://doi.org/10.21203/rs.3.rs-2025767/v1](https://doi.org/10.21203/rs.3.rs-2025767/v1)].
