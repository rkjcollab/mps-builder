# mps-builder
A tool for building a methylation profile score (MPS) from DNA methylation data.

**TO NOTE: this tool is currently under development.**


## Development Instructions

### R Environmen Setup

In order to work on development of this tool, setup the R environment:

1. Enter R:

```bash
R

```

2. In R, restore the R environment:

```r
renv::restore()

```

If you add or update an R package, install it as usual and then snapshot the
R environment:

```r
renv::install("package_name")
usethis::use_package("package_name")  # optionally declare as a dependency
renv::snapshot()

```


### Pre-commit Setup

One time before making any commits, install precommit hooks:

```bash
pre-commit install

```
