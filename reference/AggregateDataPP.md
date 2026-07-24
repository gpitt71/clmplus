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

  A square numeric matrix with at least two rows. Rows are accident
  periods, columns are development periods, and observed upper-triangle
  cells satisfy \`row + column \<= J + 1\`. Values are non-negative
  cumulative paid amounts in the source data's monetary units,
  non-decreasing across each row; unavailable cells may be \`NA\`.
  Recoveries (negative incremental payments) are not supported.

- entries.weights:

  Optional non-negative numeric \`J\` by \`J\` matrix of fitting weights
  in the same accident/development layout. \`NULL\` gives observed cells
  weight one. The first development period and missing cells are always
  zero-weighted after conversion to calendar coordinates.

- eta:

  One finite numeric value in \`(0, 1\]\`, default \`0.5\`, describing
  expected within-cell payment timing (lost exposure). It is used to
  derive exposure and to convert fitted hazards to development factors.

## Value

An \`AggregateDataPP\` list with:

- cumulative.payments.triangle:

  The input \`J\` by \`J\` cumulative paid triangle, unchanged.

- occurrance:

  A \`J\` by \`J\` matrix of incremental paid amounts in
  development-period by calendar-period coordinates. The misspelling is
  retained as a stable public field name.

- exposure:

  A \`J\` by \`J\` numeric matrix in the same calendar coordinates,
  calculated as cumulative payments minus \`(1 - eta) \* occurrence\`.

- incremental.payments.triangle:

  A \`J\` by \`J\` accident/development matrix of incremental paid
  amounts.

- fit.w:

  The \`J\` by \`J\` fitting-weight matrix in development/calendar
  coordinates.

- J:

  The integer triangle dimension.

- diagonal:

  A length-\`J\` numeric vector containing the latest observed
  cumulative diagonal in calendar representation.

- eta:

  The supplied within-cell timing scalar.

## References

Pittarello, G., Hiabu, M., & Villegas, A. M. (2023). Replicating and
extending chain-ladder via an age-period-cohort structure on the claim
development in a run-off triangle. arXiv preprint arXiv:2301.03858.

## Examples

``` r
data(sifa.mtpl)
sifa.mtpl.rtt <- AggregateDataPP(cumulative.payments.triangle=sifa.mtpl)
```
