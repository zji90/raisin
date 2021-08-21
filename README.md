RAISIN: Regression Analysis in Single-cell RNA-Seq with multiple samples
====

## Overview
RAISIN is a software tool to perform Regression Analysis In SINgle-cell RNA-seq datasets with multiple samples. RAISIN performs DE analysis using a flexible mixed effects regression framework that accounts for both sample-level and cell- level variances (Fig. 1a). The classical linear mixed effects model (LMM) do not consider small sample size or small cell number in rare cell populations, which are common in scRNA-seq studies and can lead to poor variance estimation and reduced statistical power. Fitting mixed effects models to large datasets consisting of many samples and millions of cells is also computationally challenging. To address these issues, RAISIN combines the mixed model with a hierarchical model to regularize variances, and a new model fitting algorithm is developed to efficiently handle large datasets. RAISIN can be applied to identify differential cell type proportions, differential gene expression between cell types, and differential gene expression between sample groups within each cell type.

## RAISIN Installation

RAISIN software can be installed via Github.
Users should have R installed on their computer before installing RAISIN. R can be downloaded here: http://www.r-project.org/.
To install the latest version of RAISIN package via Github, run following commands in R:
```{r }
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("zji90/raisin")
```
## Contact the Author
Author: Zhicheng Ji, Wenpin Hou, Hongkai Ji

Report bugs and provide suggestions by sending email to:

Maintainer: Zhicheng Ji (zhicheng.ji@duke.edu)

Or open a new issue on this Github page
