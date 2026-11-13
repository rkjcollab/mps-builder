# mps-builder
A tool for building a methylation profile score (MPS) from DNA methylation data.

**TO NOTE: this tool is currently under development.**


## Development Instructions

### R Environment Setup

In order to work on development of this tool, setup the R environment:

1. Enter R:

```bash
R

```

2. In R, restore the R environment:

```r
renv::restore()

```

If you add or update an R package, install it as usual and then snapshot the R
environment:

```r
renv::install("package_name")  # installs to project library
usethis::use_package("package_name")  # adds to Imports in DESCRIPTION
renv::snapshot()

```

### Pre-commit Setup

One time before making any commits, install precommit hooks:

```bash
pre-commit install

```

The hooks run general checks, and styler and lintr for R files at each commit.
To run on all files:

```bash
pre-commit run --all-files

```


### Development Workflow

Load the package from source:

```r
devtools::load_all()

```

To use a package function inside another function, add it to the `@importFrom`
lines in `R/mpsbuilder-package.R` or call it as `package::function()`.


### Documentation

Function documentation is written as roxygen2 `#'` comments above each function
in `R/`. Add `#' @export` to make a function available to package users, or
`#' @keywords internal` for internal helpers.

roxygen2 does not run automatically. After adding, renaming, or deleting
functions, or editing `#'` comments, regenerate `NAMESPACE` and `man/`:

```r
devtools::document()

```

Commit the updated `NAMESPACE` and `man/` files with your change. CI fails if
they are out of date.


### Running Tests

Unit tests use small simulated data files in `tests/testthat/testdata/`. Run
them from R:

```r
devtools::test()

```

To run a single test file:

```r
devtools::test(filter = "filter_cpgs")

```

To regenerate the test data files, from the repository root:

```bash
Rscript data-raw/make_testdata.R

```

Tests for known issues in the original code are marked with `skip()` and the
reason; remove the `skip()` once the issue is fixed.


### Package Check

Before opening a pull request, run R CMD check. It rebuilds the documentation,
runs the tests, and checks the package structure:

```r
devtools::check()

```

### Continuous Integration

GitHub Actions runs on every pull request and on pushes to `main`:

- **R-CMD-check**: confirms `NAMESPACE` and `man/` are up to date with roxygen2,
  then runs R CMD check (including the tests) using the `renv.lock` environment.
- **Lint**: runs the pre-commit hooks on all files. Fails on any lint.
