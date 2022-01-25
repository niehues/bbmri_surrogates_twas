# README

This repository contains the code used to run data analyses for *Metabolic predictors of phenotypic traits can replace and complement measured clinical variables in transcriptome-wide association studies* by Anna Niehues, Daniele Bizzarri, Marcel J.T. Reinders, P. Eline Slagboom, Alain J. van Gool, Erik B. van den Akker, and Peter A.C. 't Hoen, with the BBMRI-NL BIOS and Metabolomics Consortia.

## Installation

To run analyses in the [SURF](https://www.surf.nl/en) Research Cloud of the [BBMRI-NL BIOS consortium](https://www.bbmri.nl/acquisition-use-analyze/bios), create a *BBMRI BIOS flavoured R-Studio version 4.0.3* workspace and clone this repository. Run the RStudio server, install the R package [`renv`](https://cran.r-project.org/package=renv) - `install.packages("renv")`, open the R project *surrogate_vs_reported\surrogate_vs_reported.Rproj*, and restore the project's dependencies - `renv::restore()`. 

## Executing the workflow

This project uses [`drake`](https://cran.r-project.org/package=drake) for workflow management. Run the script *surrogate_vs_reported\run.R* to run the complete pipeline and render markdown reports with the results.

See *surrogate_vs_reported\interactive.R* for example commands to run `drake` workflow in an interactive R session.

## Troubleshooting and help

### Using `renv` package environment

A brief overview of commands is given below. See [this introduction to the renv workflow](https://rstudio.github.io/renv/articles/renv.html#workflow-1) for more detail.

```{r}
renv::init()      # initialize the project package environment, create lockfile
renv::snapshot()  # save current state of library to lockfile `renv.lock`
renv::restore()   # restore environment from lockfile (e.g. in a new VM)
```

### Manual installation of Bioconductor libraries

This might be necessary if packages are not automatically restored from *surrogate_vs_reported\renv.lock*.

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("preprocessCore", update = FALSE)
BiocManager::install("bacon", update = FALSE)
BiocManager::install("limma", update = FALSE)
BiocManager::install("edgeR", update = FALSE)
BiocManager::install("fgsea", update = FALSE)
BiocManager::install("org.Hs.eg.db", update = FALSE)
BiocManager::install("reactome.db", update = FALSE)
renv::snapshot()
```

### Checking data availability via [`BBMRIomics`](https://bbmri-nl.github.io/BBMRIomics/index.html)

```{r}
library(BBMRIomics)
data(package = "BBMRIomics")
# example data load
# BBMRIomics::bbmri.data(rnaSeqData_ReadCounts_NTR_Freeze2_unrelated)
```




