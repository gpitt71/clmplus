# Fit Chain Ladder Plus to reverse time triangles.

Default method to fit Chain Ladder plus models.

## Usage

``` r
# Default S3 method
clmplus(
  AggregateDataPP,
  hazard.model = NULL,
  link = c("log", "logit"),
  staticAgeFun = TRUE,
  periodAgeFun = "NP",
  cohortAgeFun = NULL,
  effect_log_scale = TRUE,
  constFun = function(ax, bx, kt, b0x, gc, wxt, ages) list(ax = ax, bx = bx, kt = kt, b0x
    = b0x, gc = gc),
  ...
)
```

## Arguments

- AggregateDataPP:

  `AggregateDataPP` object, reverse time triangle to be fitted.

- hazard.model:

  `character`, hazard model supported from our package. The model can be
  chosen from:

  - 'a': Age model, this is equivalent to the Mack chain-ladder.

  - 'ac': Age and cohort effects.

  - 'ap': Age and cohort effects.

  - 'apc': Age cohort and period effects.

- link:

  `character`, defines the link function and random component associated
  with the mortality model. `"log"` would assume that deaths follow a
  Poisson distribution and use a log link while `"logit"` would assume
  that deaths follow a Binomial distribution and a logit link. To be
  disregarded unless the practitioner specifies his own hazard model in
  StMoMo.

- staticAgeFun:

  `logical`, indicates if a static age function \\\alpha_x\\ is to be
  included. To be disregarded unless the practitioner specifies his own
  hazard model in StMoMo.

- periodAgeFun:

  `list`, a list of length \\N\\ with the definitions of the period age
  modulating parameters \\\beta_x^{(i)}\\. Each entry can take values:
  `"NP"` for non-parametric age terms, `"1"` for \\\beta_x^{(i)}=1\\ or
  a predefined parametric function of age (see details). Set this to
  `NULL` if there are no period terms in the model. To be disregarded
  unless the practitioner specifies his own hazard model in StMoMo.

- cohortAgeFun:

  `character` or `function`, defines the cohort age modulating parameter
  \\\beta_x^{(0)}\\. It can take values: `"NP"` for non-parametric age
  terms, `"1"` for \\\beta_x^{(0)}=1\\, a predefined parametric function
  of age (see details) or `NULL` if there is no cohort effect. To be
  disregarded unless the practitioner specifies his own hazard model in
  StMoMo.

- effect_log_scale:

  `logical`, whether effects should be on the logarithmic scale. By
  default, `TRUE`.

- constFun:

  `function`, it defines the identifiability constraints of the model.
  It must be a function of the form
  `constFun <- function(ax, bx, kt, b0x, gc, wxt, ages)` taking a set of
  fitted model parameters and returning a list
  `list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)` of the model
  parameters with the identifiability constraints applied. If omitted no
  identifiability constraints are applied to the model. To be
  disregarded unless the practitioner specifies his own hazard model in
  StMoMo.

- ...:

  parameters to be passed to clmplus.

## Value

No return value, called to pass method `clmplus.AggregateDataPP`. See
`clmplus.AggregateDataPP` documentation.

## References

Pittarello, Gabriele, Munir Hiabu, and Andrés M. Villegas. "Replicating
and extending chain ladder via an age-period-cohort structure on the
claim development in a run-off triangle." arXiv preprint
arXiv:2301.03858 (2023).

Hiabu, Munir. “On the relationship between classical chain ladder and
granular reserving.” Scandinavian Actuarial Journal 2017 (2017): 708 -
729.
