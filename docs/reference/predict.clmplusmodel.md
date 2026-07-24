# Predict the Reserve using Chain Ladder Plus Models

Predict the lower triangle with a `clmplus` model.

## Usage

``` r
# S3 method for class 'clmplusmodel'
predict(
  object,
  gk.fc.model = "a",
  ckj.fc.model = "a",
  gk.order = c(1, 1, 0),
  ckj.order = c(0, 1, 0),
  forecasting_horizon = NULL,
  constrained_development_factors = FALSE,
  ...
)
```

## Arguments

- object:

  `clmplusmodel`, Model to predict from.

- gk.fc.model:

  `character`, model to forecast the cohort component for the last
  accident period. It can be either arima ('a') or linear model ('l').
  Disregarded for models that do not have a cohort effect.

- ckj.fc.model:

  `character`, model to forecast the calendar period effect. It can be
  either arima ('a') or linear model ('l'). Disregarded for models that
  do not have a period effect.

- gk.order:

  `integer`, order of the arima model with drift for the accident year
  effect extrapolation. Default to (1,1,0).

- ckj.order:

  `integer`, order of the arima model with drift for the calendar year
  effect extrapolation. Default to (0,1,0).

- forecasting_horizon:

  `integer`, between 1 and the triangle width. Calendar periods ahead
  for the predictions. Default predictions are to run-off.

- constrained_development_factors:

  `logical`, if `TRUE` the predict function will set negative
  development factors to 1.

- ...:

  Extra arguments to be passed to the predict function.

## Value

Returns the following output:

- reserve:

  `numeric` The reserve for each accident period.

- ultimate_cost:

  `numeric` The ultimate cost for each accident period.

- full_triangle:

  `matrix array` The complete run-off triangle of cumulative payments,
  it includes the (input) upper triangle and the predicted (output)
  lower triangle.

- lower_triangle:

  `matrix array` The predicted lower triangle of cumulative payments.

- development_factors_predicted:

  `matrix array` The predicted lower triangle of the extrapolated
  development factors.

- apc_output:

  `list` The following output from the age-period-cohort representation:
  `model.fit` (`fitStMoMo`) age-period-cohort model fit. `alphaij`
  (`matrix array`) predicted claim development. `lower_triangle_apc`
  (`matrix array`) predicted lower triangle of cumulative payments in
  age-period-cohort form. `development_factors_apc` (`matrix array`)
  development factors in age-period-cohort representation.

## References

Pittarello, Gabriele, Munir Hiabu, and Andrés M. Villegas. "Replicating
and extending chain ladder via an age-period-cohort structure on the
claim development in a run-off triangle." arXiv preprint
arXiv:2301.03858 (2023).
