[![R-CMD-check](https://github.com/gpitt71/clmplus/actions/workflows/r-checkrelease.yml/badge.svg)](https://github.com/gpitt71/clmplus/actions/workflows/r-checkrelease.yml)
[![R-hub](https://github.com/gpitt71/clmplus/actions/workflows/rhub.yaml/badge.svg)](https://github.com/gpitt71/clmplus/actions/workflows/rhub.yaml)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/grand-total/clmplus?color=blue)](https://r-pkg.org/pkg/clmplus)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10987181.svg)](https://doi.org/10.5281/zenodo.10987181)


# clmplus

`clmplus` is an `R` package for implementing the age-period-cohort models for the claim development presented in the manuscript 'Replicating and extending chain-ladder via an age-period-cohort structure on the claim development in a run-off triangle' <doi:10.48550/arXiv.2301.03858>. 

## Our models

`clmplus` relies on the powerful `StMoMo` package. Users can either rely on our default models or set their own configuration for the claim development.


|      Model      | Lexis dimension                 |Claims reserving                  |
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

## Get Started

In this brief example, we work with the `sifa.mtpl` data from the `clmplus` package. Further examples can be found in the [package vignettes](https://github.com/gpitt71/clmplus/tree/main/vignettes). The data set of cumulative claim payments is transformed into an `AggregateDataPP` object that pre-processes the data for claim development modelling.

```
library(clmplus)

data ("sifa.mtpl")
dataset = sifa.mtpl
datapp = AggregateDataPP(cumulative.payments.triangle = dataset, eta= 1/2)
```

Our models can be fit with the `clmplus` function.

```
a.model.fit=clmplus(datapp,
                 hazard.model = "a") # age-model replicates the chain ladder
                 
ac.model.fit=clmplus(datapp,
                 hazard.model = "ac")

ap.model.fit=clmplus(datapp,
                 hazard.model = "ap")

apc.model.fit=clmplus(datapp,
                  hazard.model = "apc")

```

The `plot` function can be be used to explore the scaled deviance residuals of fitted models. Below, an example for the age-period-cohort (`apc`) model for the claim development. 

```
plot(apc.model.fit)
```

Predictions are performed with the `predict` function.

```
a.model=predict(a.model.fit)
                 
# clmplus reserve (age model)
sum(a.model$reserve)
#226875.5


ac.model=predict(ac.model.fit,
                 gk.fc.model = 'a',
                 gk.order = c(1,1,0))
                 
# clmplus reserve (age-cohort model)
sum(ac.model$reserve)
#205305.7

ap.model= predict(ap.model.fit,
                 ckj.fc.model = 'a',
                 ckj.order = c(0,1,0))

# clmplus reserve (age-period model)
sum(ap.model$reserve)
#215602.8
          
                 
apc.model= predict(apc.model.fit,
                  gk.fc.model = 'a',
                  ckj.fc.model = 'a',
                  gk.order = c(1,1,0),
                  ckj.order = c(0,1,0))
# clmplus reserve (age-period-cohort model)
sum(apc.model$reserve)
#213821.6
```

The fitted effect (and extrapolated) effects can be inspected with the `plot` function. We continue below the example with the `apc` model.

```
plot(apc.model)
```




