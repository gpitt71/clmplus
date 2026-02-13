# Plot the hazard model fitted and forecasted parameters

This function allows to define the behavior of the triangle payments.

## Usage

``` r
# S3 method for class 'clmpluspredictions'
plot(x, cy.type = "fe", ...)
```

## Arguments

- x:

  `clmpluspredictions`, Model effects (fitted and extrapolated) to be
  plotted.

- cy.type:

  `character`, whether to show fitted period effect with or without
  extrapolatio Default is "fe", standing for fitted and extrapolated.
  Alternative is to specify "f" for fitted effect.

- ...:

  Arguments to be passed to plot.

## Value

No return value, plots coefficients of the hazard models.

## References

Pittarello, G., Hiabu, M., & Villegas, A. M. (2023). Replicating and
extending chain-ladder via an age-period-cohort structure on the claim
development in a run-off triangle. arXiv preprint arXiv:2301.03858.

## Examples

``` r
data(sifa.mtpl)
sifa.mtpl.rtt <- AggregateDataPP(cumulative.payments.triangle=sifa.mtpl)
clm.fit<-clmplus(sifa.mtpl.rtt, 'a')
#> Warning: StMoMo: 66 missing values which have been zero weighted
#> StMoMo: The following ages have been zero weigthed: 1 
#> StMoMo: The following years have been zero weigthed: 1 
#> StMoMo: The following cohorts have been zero weigthed: -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 11 
#> StMoMo: Start fitting with gnm
#> Warning: non-integer x = 860.063891
#> Warning: non-integer x = 458.156066
#> Warning: non-integer x = 559.340054
#> Warning: non-integer x = 281.869679
#> Warning: non-integer x = 456.270650
#> Warning: non-integer x = 727.142070
#> Warning: non-integer x = 49353.589600
#> Warning: non-integer x = 50605.820000
#> Warning: non-integer x = 29251.193380
#> Warning: non-integer x = 36106.295210
#> Warning: non-integer x = 40125.390780
#> Warning: non-integer x = 44498.942240
#> Warning: non-integer x = 45490.189570
#> Warning: non-integer x = 48040.321560
#> Warning: non-integer x = 49991.357650
#> Warning: non-integer x = 49694.294990
#> Warning: non-integer x = 20880.666980
#> Warning: non-integer x = 18304.246140
#> Warning: non-integer x = 18603.713030
#> Warning: non-integer x = 12463.983790
#> Warning: non-integer x = 13441.333110
#> Warning: non-integer x = 12951.176400
#> Warning: non-integer x = 15370.423630
#> Warning: non-integer x = 15339.426370
#> Warning: non-integer x = 17842.901760
#> Warning: non-integer x = 19570.203660
#> Warning: non-integer x = 10047.067150
#> Warning: non-integer x = 8201.559214
#> Warning: non-integer x = 8833.487780
#> Warning: non-integer x = 5144.085289
#> Warning: non-integer x = 5868.234632
#> Warning: non-integer x = 6033.755464
#> Warning: non-integer x = 5593.853546
#> Warning: non-integer x = 5478.076757
#> Warning: non-integer x = 7035.199194
#> Warning: non-integer x = 3933.606063
#> Warning: non-integer x = 5749.890057
#> Warning: non-integer x = 4714.168250
#> Warning: non-integer x = 2726.947823
#> Warning: non-integer x = 2881.972323
#> Warning: non-integer x = 3009.655856
#> Warning: non-integer x = 2615.516340
#> Warning: non-integer x = 2540.639658
#> Warning: non-integer x = 2906.368627
#> Warning: non-integer x = 2725.682936
#> Warning: non-integer x = 3312.691467
#> Warning: non-integer x = 2359.362462
#> Warning: non-integer x = 2421.987896
#> Warning: non-integer x = 1264.396267
#> Warning: non-integer x = 1984.311792
#> Warning: non-integer x = 2136.804699
#> Warning: non-integer x = 1294.338023
#> Warning: non-integer x = 2266.599873
#> Warning: non-integer x = 1333.856105
#> Warning: non-integer x = 918.004039
#> Warning: non-integer x = 1249.952517
#> Warning: non-integer x = 1134.706143
#> Warning: non-integer x = 1184.355428
#> Warning: non-integer x = 1124.336355
#> Warning: non-integer x = 1238.429974
#> Warning: non-integer x = 1075.978340
#> Warning: non-integer x = 733.741025
#> Warning: non-integer x = 904.371165
#> Warning: non-integer x = 872.633331
#> Warning: non-integer x = 941.440475
#> StMoMo: Finish fitting with gnm
clm <- predict(clm.fit)
plot(clm)
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_line()`).

```
