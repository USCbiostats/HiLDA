# HiLDA: "*Hi*erarchical *L*atent *D*irichlet *A*llocation" 
[![](https://img.shields.io/badge/release%20version-1.7.5-green.svg)](https://www.bioconductor.org/packages/HiLDA) 
[![](https://img.shields.io/badge/download-2676/total-blue.svg)](https://bioconductor.org/packages/stats/bioc/HiLDA) 
[![](https://img.shields.io/badge/doi-10.7717/peerj.7557-yellow.svg)](https://doi.org/10.7717/peerj.7557) [![](https://raw.githubusercontent.com/USCbiostats/badges/master/tommy-image-badge.svg)](https://image.usc.edu)

## Introduction

The R package `HiLDA` is developed under the Bayesian framework to allow 
statistically testing whether there is a change in the mutation burdens of 
mutation signatures between two groups. The mutation signature is defined based 
on the independent model proposed by Shiraishi's et al. 

## Paper

- **Yang Z**, Pandey P, Shibata D, Conti DV, Marjoram P, Siegmund KD. 2019. HiLDA: a statistical approach to investigate differences in mutational signatures. PeerJ 7:e7557 https://doi.org/10.7717/peerj.7557

## Installation 

Now you can download the pacakge from Bioconductor https://bioconductor.org/packages/devel/bioc/html/HiLDA.html

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("HiLDA")
```

## Documentation

To view documentation for the version of this package installed in your system, start R and enter:

```
browseVignettes("HiLDA")
```

## An introduction to HiLDA

Tutorials: 
https://bioconductor.org/packages/devel/bioc/vignettes/HiLDA/inst/doc/HiLDA.html

R Scripts:
https://bioconductor.org/packages/devel/bioc/vignettes/HiLDA/inst/doc/HiLDA.R

## Funding
This work was supported by NCI grant numbers 5P30 CA014089 and P01 CA196569. 
