# recountmethylationManuscriptSupplement

[![DOI](https://zenodo.org/badge/302553109.svg)](https://zenodo.org/badge/latestdoi/302553109)

This repo contains supplemental methods, data, code, and scripts for the manuscript "Human methylome variation across Infinium 450K data on the Gene Expression Omnibus."

## Install as an R package

This repo is designed to work like a typical R packaage. You can install and load the repo from an active R session with:

```
pkg.name <- "recountmethylationManuscriptSupplement"
repo.path <- file.path("metamaden", pkg.name)

devtools::install_github(repo.path, build_vignettes = TRUE)

library(recountmethylationManuscriptSupplement)
```

Documentation, including Supplemental Information and vignettes that reproduce the manuscript figures, can be found using:

```
doc.path <- system.file("docfiles", package = pkg.name)
list.files(doc.path)
```

## Resources of interest

The following resources include data, code, and scripts to reproduce findings in the main manuscript.

* `recountmethylation` -- Bioconductor package to access, query, and analyze full database 
compilations and sample metadata. https://bioconductor.org/packages/release/bioc/html/recountmethylation.html

* recountmethylationManuscriptSupplement -- GitHub repo for manuscript supplemental files, scripts, and data. https://github.com/metamaden/recountmethylationManuscriptSupplement

* figshare project -- Supplemental files, metadata, and DOIs for the manuscript https://figshare.com/account/home#/projects/90758

*  recount.bio/data -- Location of large supplemental data files. https://recount.bio/data/recountmethylation_manuscript_supplement/

* `recount-methylation-server` -- Server software for identifying, downloading, and managing IDATs and SOFT files. https://github.com/metamaden/recount-methylation-server

* `rmpipipeline` -- GitHub repo with resources for processing GEO IDATs and SOFT files.
https://github.com/metamaden/rmpipeline
