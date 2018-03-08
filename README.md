---
title: "patelGBMSC -- CONQUER quantification of single-cell RNA-seq in glioblastoma"
author: "Vincent J. Carey, stvjc at channing.harvard.edu"
date: "`r format(Sys.time(), '%B %d, %Y')`"
vignette: >
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{patelGBMSC -- a single-cell RNA-seq dataset in glioblastoma}
  %\VignetteEncoding{UTF-8}
output:
  BiocStyle::html_document:
    highlight: pygments
    number_sections: yes
    theme: united
    toc: yes
---

```{r setup,echo=FALSE,results="hide"}
suppressPackageStartupMessages({
suppressMessages({
library(patelGBMSC)
})
})
```

# Introduction

[Patel et al. 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4123637/)
describe a single-cell RNA-seq study of several glioblastoma
samples.  The data were reprocessed with the 
[CONQUER](http://imlspenticton.uzh.ch:3838/conquer/)
pipeline (see the [QC report](http://imlspenticton.uzh.ch/robinson_lab/conquer/report-multiqc/GSE57872_multiqc_report.html)).

The rds file distributed by CONQUER is large as it includes
multiple gene-level and transcript-level quantifications.
As of Oct 30 2017, the CONQUER distribution does not include
sample-level information beyond the GSM identifier.  This
package includes a smaller image of the data (the `count_lstpm`
quantifications, that are estimated counts created using the
salmon algorithm, rescaled to account for library size).
The data image is over 200MB, so the `r Biocpkg("BiocFileCache")` 
discipline is used to perform a one-time download, insertion
and bookkeeping in cache; the `loadPatel` function takes
care of the download and retrieval from cache as appropriate.

# Quick view of the data

We'll randomly sample 5000 genes to reduce runtime
in this vignette.  We filter down to the 430 patient samples
that passed quality control.

```{r getdat}
library(patelGBMSC)
patelGeneCount = loadPatel()
#
# use metadata on sample QC to exclude failed samples
#
qdrop = grep("excluded", patelGeneCount$description) # QC issues
patelGeneCount = patelGeneCount[,-qdrop]
#
# drop gliospheres
#
ispat = grep("MGH", patelGeneCount$characteristics_ch1)
patelGeneCount = patelGeneCount[,ispat]
patelERCCCount = patelGeneCount  # save for ERCC check later
#
# drop ERCC spikeins
#
patelGeneCount = patelGeneCount[-grep("ERCC", rownames(patelGeneCount)),] 
#
# keep ERCC spikeins
#
patelERCCCount = patelERCCCount[grep("ERCC", rownames(patelERCCCount)),]
#
# derive patient code
#
patelGeneCount$sampcode = factor(gsub("patient id: ", "", patelGeneCount$characteristics_ch1))
tcol = as.numeric(tfac <- factor(patelGeneCount$sampcode))
patelERCCCount$sampcode = factor(gsub("patient id: ", "", patelERCCCount$characteristics_ch1))
etcol = as.numeric(tfac <- factor(patelERCCCount$sampcode))
#
# sample 5000 genes for t-SNE
#
set.seed(1234)
samp = assay(patelGeneCount[sample(1:nrow(patelGeneCount), size=5000),])
library(Rtsne)
RTL = Rtsne(t(log(samp+1)))
myd = data.frame(ts1=RTL$Y[,1], ts2=RTL$Y[,2], 
        code = patelGeneCount$sampcode, tcol=tcol)
library(ggplot2)
ggplot(myd, aes(x=ts1, y=ts2, group=code, colour=code)) + geom_point() +
  ggtitle("t-SNE for 5000 randomly chosen genes in five GBM scRNA samples")
```

# Some views of ERCC spikeins

## t-SNE
```{r liker}
ercc = assay(patelERCCCount)
set.seed(1234)
ERTL = Rtsne(t(log(ercc+1)))
ed = data.frame(ets1=ERTL$Y[,1], ets2=ERTL$Y[,2],
    code = patelERCCCount$sampcode, tcol=etcol)
ggplot(ed, aes(x=ets1, y=ets2, group=code, colour=code)) + geom_point() +
  ggtitle("t-SNE for ERCC spikeins in five GBM scRNA samples")
```

## PCA

```{r lkr2}
pcs = prcomp(t(log(ercc+1)))
plot(pcs$x[,1], pcs$x[,2], pch=19, col=etcol)
boxplot(split(pcs$x[,1], patelERCCCount$sampcode))
boxplot(split(pcs$x[,2], patelERCCCount$sampcode))
```
