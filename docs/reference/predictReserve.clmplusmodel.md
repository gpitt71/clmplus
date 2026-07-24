# Convert clmplus predictions to the ReSurv reserving interface

Convert clmplus predictions to the ReSurv reserving interface

## Usage

``` r
# S3 method for class 'clmplusmodel'
predictReserve(object, granularity = NULL, lower_triangle_only = TRUE, ...)

# S3 method for class 'clmpluspredictions'
predictReserve(object, granularity = NULL, lower_triangle_only = TRUE, ...)
```

## Arguments

- object:

  A \`clmplusmodel\` or \`clmpluspredictions\` object.

- granularity:

  Ignored; retained for ReSurv compatibility.

- lower_triangle_only:

  Return only forecast lower-triangle cells.

- ...:

  Arguments passed to \[predict.clmplusmodel()\].

## Value

A data frame containing \`AP\`, \`DP\`, \`CP\`, and \`IBNR\`.
