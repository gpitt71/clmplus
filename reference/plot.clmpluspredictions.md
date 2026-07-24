# Plot the hazard model fitted and forecasted parameters

This function allows to define the behavior of the triangle payments.

## Usage

``` r
# S3 method for class 'clmpluspredictions'
plot(x, cy.type = "fe", ...)
```

## Arguments

- x:

  A \`clmpluspredictions\` object returned by
  \[predict.clmplusmodel()\].

- cy.type:

  Either \`"fe"\` (the default) to include extrapolated calendar effects
  or \`"f"\` to show only fitted calendar effects.

- ...:

  Reserved; currently ignored.

## Value

A \`gtable\` containing one ggplot panel for each effect included in the
fitted model. The table is returned visibly after being drawn.

## References

Pittarello, G., Hiabu, M., & Villegas, A. M. (2023). Replicating and
extending chain-ladder via an age-period-cohort structure on the claim
development in a run-off triangle. arXiv preprint arXiv:2301.03858.

## Examples

``` r
data(sifa.mtpl)
sifa.mtpl.rtt <- AggregateDataPP(cumulative.payments.triangle=sifa.mtpl)
clm.fit<-clmplus(sifa.mtpl.rtt, 'a')
clm <- predict(clm.fit)
plot(clm)

```
