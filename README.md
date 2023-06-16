
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SplineHazardRegression

<!-- badges: start -->
<!-- badges: end -->

There is no R software available for the direct estimation of hazards
with uncertainty estimation (confidence intervals).

The primary goal of SplineHazardRegression is to make available the
methods for flexible estimation of hazards using (cubic) b-splines
published in Philip S. Rosenberg. “Hazard Function Estimation Using
B-Splines” In Biometrics, Vol. 51, No. 3 (Sep., 1995), pp. 874-887
<https://doi.org/10.2307/2532989> The input data is time-to-event data
(e.g. time to death), possibly right-censored and with late entries
(both meaning that patients are followed-up for unequal times).

The package also allows the flexible estimation of the cumulative hazard
and cumulative survival functions, as well as the compuatation of
aggregate measures for those (average, median, interquartile range,
etc).

Different methods for the automatic selection of knots and for variance
estimation are implemented.

## Installation

You can install the development version of SplineHazardRegression from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fespuny/SplineHazardRegression")
```

## Example

This is a basic example which shows you how to simulate time-to-event
data and fit a hazard function, deriving then cumulative hazard and
survival estimates.

### Data Simulation

``` r
library(SplineHazardRegression)
## simulation parameters
   knots = c(0, 1, 3, 6, 10, NA, NA)
   betac = 1 * c(0.05, 0.05, 0.05, 0.05, 0.40, 0.1, 0.05)
   HParm = data.frame(knots, betac) # 'A Simple B-Spline'
   cll = c(0, 5)
   cup = c(5, 10)
   cih = c(0.0125, 0.025)
   CParm = data.frame(cll, cup, cih) # 'Light Censoring'
## calculate simulation true hazard and censoring distributions
   INPUTS = etsim_inputs( HParam=HParm, CParam=CParm, SampleSize = 101 )
## simulate time-to-event data using true distribution
   SimDat = etsim(INPUTS)
## histogram of the fully observed hazard data (gray) and censored observations (light blue)
hist( SimDat$time[ which(SimDat$status==1)], main="", xlab="t", breaks="Freedman-Diaconis", xlim=c(0,10) )
hist( SimDat$time[ which(SimDat$status==0)], main="", xlab="t", breaks="Freedman-Diaconis", xlim=c(0,10),add=TRUE, col=rgb(173,216,230,max=255,alpha=100) )
```

<img src="man/figures/README-simulation-1.png" width="100%" />

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
## Fit a cubic B-spline regression model estimating the needed knots
   timeout = seq( 0, 10, length.out = 101 )
   Result = hspcore(yd=SimDat, ORDER=4, knots=c(0,1,3,6,10), time=timeout )
#> [1] "K= 3 DOF= 7 knots= 0 0.1 0.3 0.6 1"
#> [1] "L-BFGS-B Convergence (0=yes,1=maxiter)? 52"
#> [1] "L-BFGS-B Warning: ERROR: ABNORMAL_TERMINATION_IN_LNSRCH"
#> [1] "PORT Convergence (0=yes)? 1"
#> [1] "PORT Warning: false convergence (8)"
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
