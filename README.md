# infuse <img src="man/figures/logo.png" align="right" height="139" />

The goal of `infuse` is to provide fast, mathematically consistent inference framework for Win statistics (net treatment benefit NTB, win ratio WR and a generalisation of the NNT) using influence functions.

## Installation
```r
# install.packages("devtools")
devtools::install_github("AlbertoAlvarezIglesias2019/infuse")
```

## Documentation
The package is designed around a three-stage lifecycle: sow(), harvest(), and reap().

For a full tutorial with real-world examples (Continuous, Binary, and Survival data), see the **Quick Start Guide**:

```r
vignette("infuse_start_guide", package = "infuse")
```
or just

```r
browseVignettes("infuse")
```
