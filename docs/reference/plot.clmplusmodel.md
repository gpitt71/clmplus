# Plot the hazard model residuals

This function allows to plot the hazard model residuals on the triangle
payments.

## Usage

``` r
# S3 method for class 'clmplusmodel'
plot(x, heat.lim = c(-2.5, 2.5), ...)
```

## Arguments

- x:

  A fitted \`clmplusmodel\` object.

- heat.lim:

  A length-two numeric vector giving the lower and upper fill scale
  limits for scaled deviance residuals.

- ...:

  Reserved; currently ignored.

## Value

A \`ggplot\` object showing scaled deviance residuals in accident-year
by development-year triangle form.

## References

Pittarello, Gabriele, Munir Hiabu, and Andrés M. Villegas. "Replicating
and extending chain ladder via an age-period-cohort structure on the
claim development in a run-off triangle." arXiv preprint
arXiv:2301.03858 (2023).

## Examples

``` r
data(sifa.mtpl)
sifa.mtpl.rtt <- AggregateDataPP(cumulative.payments.triangle=sifa.mtpl)
clm.fit<-clmplus(sifa.mtpl.rtt, 'a')
#> Warning: StMoMo: 66 missing values which have been zero weighted
plot(clm.fit)

```
