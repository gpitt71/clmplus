[![R-CMD-check](https://github.com/gpitt71/clmplus/actions/workflows/r-checkrelease.yml/badge.svg)](https://github.com/gpitt71/clmplus/actions/workflows/r-checkrelease.yml)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/grand-total/clmplus?color=blue)](https://r-pkg.org/pkg/clmplus)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10987181.svg)](https://doi.org/10.5281/zenodo.10987181)


# clmplus

Repository GitHub that contains the code for the package `clmplus`. 

About `clmplus`:

* It adds flexibilty to the chain-ladder model: we provide a new and more versatile framework for claims reserving given the cumulative payments data.

* It achieves the ambitious objective of showing the contact point between models often used in life insurance and non-life insurance models for claims reserving.

* `clmplus` relies on the powerful `StMoMo` package. Practitioners can either use the default models we programmed for them or set their own hazard model.
Examples of model configurations we support: 

|      Model      | Life insurance                  |Non-life insurance                |
| :-------------: |:-------------------------------:|---------------------------------:|
| a               | age                             |development (chain-ladder model)  |
| ac              | age-cohort                      |development-accident              |
| ap              | age-period                      |development-calendar              |
| apc             | age-period-cohort               |development-calendar-accident     |

## Installation 

The developer version of `clmplus` can be installed from GitHub.

```
library(devtools)
devtools::install_github("gpitt71/clmplus")

```
The current version of `clmplus` can be installed from CRAN.

```
install.packages('clmplus')

```





