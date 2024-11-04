# SplineHazardRegression

**Ferran Espuny Pujol,  Philip S. Rosenberg**

Survival hazards are used in regression (Cox proportional hazards
models) but rarely reported or visualized in descriptive analyses, where
Kaplan-Meier cumulative survival and at best Nelson-Aalen cumulative
hazard are commonly used. It turns out that describing survival hazards
can be very informative, allowing to spot “instantaneous” patterns that
are not visible in cumulative measures of hazard and survival, which in
fact can be derived from hazard functions.

The primary goal of SplineHazardRegression is to make available the
methods for flexible estimation of hazards using (cubic) b-splines
published in Philip S. Rosenberg. “Hazard Function Estimation Using
B-Splines” In Biometrics, Vol. 51, No. 3 (Sep., 1995), pp. 874-887
<https://doi.org/10.2307/2532989> The input data is time-to-event data
(e.g. time to death), possibly right-censored and with late entries
(both meaning that patients are followed-up for unequal times).

The package also allows the flexible estimation of the cumulative hazard
and cumulative survival functions, and (in progress) the computation of
aggregate measures for those (average, median, interquartile range,
etc).

The automatic selection of knots and bootstrap variance estimation have
been implemented.

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
   INPUTS = etsim_inputs( HParam=HParm, CParam=CParm, SampleSize = 301 )
## simulate time-to-event data using true distribution
   SimDat = etsim(INPUTS)
## histogram of the fully observed hazard data (gray) and censored observations (light blue)
hist( SimDat$time[ which(SimDat$status==1)], main="", xlab="t", breaks="Freedman-Diaconis", xlim=c(0,10) )
hist( SimDat$time[ which(SimDat$status==0)], main="", xlab="t", breaks="Freedman-Diaconis", xlim=c(0,10),add=TRUE, col=rgb(173,216,230,max=255,alpha=100) )
```

<img src="man/figures/README-simulation-1.png" width="100%" />

### Hazard regression

``` r
## Fit a cubic B-spline regression searching for the best number of knots
   timeout = seq( 0, 10, length.out = 501 )
   Result = hspcore(yd=SimDat, ORDER=4, Exterior.knots = c(0,10), Interior.knots=NULL, SelectBestKnots = TRUE, time=timeout, Bootstrap = 200, verbose=FALSE )  
#> [1] "Automatic search for K the number of interior knots of the B-spline hazard function"
#> [1] "K=0, m2loglik=287.04, DOF=4, AICc=295.17, knots= 0 10"
#> [1] "K=1, m2loglik=282.82, DOF=5, AICc=293.03, knots= 0 4.5 10"
#> [1] "K=2, m2loglik=277.18, DOF=6, AICc=289.47, knots= 0 3.2 5.5 10"
#> [1] "K=3, m2loglik=278.18, DOF=7, AICc=292.56, knots= 0 2.1 4.5 6 10"
#> [1] "K=4, m2loglik=275.48, DOF=8, AICc=291.98, knots= 0 1.6 4 5 6.3 10"
#> [1] "K=5, m2loglik=274.06, DOF=9, AICc=292.68, knots= 0 1.4 3.2 4.5 5.5 6.6 10"
#> [1] "K=6, m2loglik=273.47, DOF=10, AICc=294.23, knots= 0 1.2 2.7 4.1 4.8 5.8 6.8 10"
#> [1] "K=7, m2loglik=272.24, DOF=11, AICc=295.15, knots= 0 0.9 2.1 3.7 4.5 5.2 6 7.1 10"
#> [1] "K=8, m2loglik=273.11, DOF=12, AICc=298.19, knots= 0 0.8 1.8 3.2 4.2 4.8 5.5 6.2 7.3 10"
#> [1] "K=9, m2loglik=271.62, DOF=13, AICc=298.89, knots= 0 0.7 1.6 2.9 4 4.5 5 5.7 6.3 7.4 10"
#> [1] "K=10, m2loglik=271.2, DOF=14, AICc=300.67, knots= 0 0.6 1.4 2.3 3.7 4.3 4.7 5.3 5.9 6.5 7.6 10"
#> [1] "K=11, m2loglik=271.12, DOF=15, AICc=302.81, knots= 0 0.6 1.4 2.1 3.2 4 4.5 4.9 5.5 6 6.6 7.7 10"
#> [1] "K=12, m2loglik=270.63, DOF=16, AICc=304.54, knots= 0 0.5 1.2 1.9 3 3.9 4.3 4.7 5.2 5.6 6.1 6.7 7.8 10"
#> [1] "K=13, m2loglik=271.15, DOF=17, AICc=307.31, knots= 0 0.4 1.2 1.7 2.7 3.6 4.1 4.5 4.8 5.4 5.8 6.3 6.8 7.9 10"
#> [1] "K=14, m2loglik=270.39, DOF=18, AICc=308.81, knots= 0 0.4 1 1.6 2.3 3.2 4 4.3 4.7 5 5.5 5.9 6.3 7 7.9 10"
#> [1] "SEARCH RESULT: We use 2 interior B-spline knots"
#> [1] "K= 2 DOF= 6 knots= 0 3.2 5.5 10"
#> [1] "Variance estimation using bootstrap"
```

<img src="man/figures/README-result-figures-side-1.png" width="50%" /><img src="man/figures/README-result-figures-side-2.png" width="50%" /><img src="man/figures/README-result-figures-side-3.png" width="50%" /><img src="man/figures/README-result-figures-side-4.png" width="50%" />

## PACKAGE DEVELOPMENT NOTE

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. In that case,
don’t forget to commit and push the resulting figure files, so they
display on GitHub and CRAN.
