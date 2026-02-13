# Pre-process Run-Off Triangles

Pre-process Run-Off Triangles.

## Usage

``` r
AggregateDataPP(
  cumulative.payments.triangle,
  entries.weights = NULL,
  eta = 1/2
)
```

## Arguments

- cumulative.payments.triangle:

  `triangle matrix` or `matrix array` object, input triangle of
  cumulative payments.

- entries.weights:

  `triangle matrix` or `matrix array` model entries weights.

- eta:

  `numeric`, individual claims exposure in the cell, also known as lost
  exposure. It must be in the interval (0,1\].

## Value

An object of class `AggregateDataPP`. Lists the following elements:

- cumulative.payments.triangle:

  `triangle matrix` object, input triangle of cumulative payments.

- occurrance:

  `matrix array` object, the occurrence derived from the input triangle.

- exposure:

  `matrix array` object, the exposure derived from the input triangle,
  under the `eta` claims arrival assumption.

- incremental.payments.triangle:

  `triangle matrix` object, incremental payments derived from the input.

- fit.w:

  `matrix array` object, the weights used during estimation.

- J:

  `integer`, Run-off triangle dimension.

- diagonal:

  `numeric`, cumulative payments last diagonal.

- eta:

  `numeric`, Expected time-to-event in the cell. I.e., lost exposure.

## References

Pittarello, G., Hiabu, M., & Villegas, A. M. (2023). Replicating and
extending chain-ladder via an age-period-cohort structure on the claim
development in a run-off triangle. arXiv preprint arXiv:2301.03858.

## Examples

``` r
data(sifa.mtpl)
sifa.mtpl.rtt <- AggregateDataPP(cumulative.payments.triangle=sifa.mtpl)
```
