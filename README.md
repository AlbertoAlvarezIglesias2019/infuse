# infuse <img src="man/figures/logo.png" align="right" height="139" />

The goal of `infuse` is to provide a fast, mathematically consistent alternative framework for generalised pairwise comparisons. By leveraging empirical survival functions and influence functions, `infuse` provides rapid inference for win statistics, including:

* Net Treatment Benefit (NTB)

* Win Ratio (WR)

* Generalized Number Needed to Treat (GNNT)

## Installation
You can install the development version of infuse from GitHub.
**Note**: Building the vignettes may take a few minutes as it executes the coding examples.

```r
# install.packages("devtools")
devtools::install_github("AlbertoAlvarezIglesias2019/infuse", build_vignettes = TRUE)
```

## Documentation
The package is designed around a three-stage lifecycle: sow(), harvest(), and reap().

For a detailed tutorial covering real-world examples (Continuous, Binary, and Survival data), please refer to the **Quick Start Guide**:

```r
library(infuse)
vignette("infuse_start_guide", package = "infuse")
```
Alternatively, you can browse all available documentation by running:
```r
library(infuse)
browseVignettes("infuse")
```
