pkgname <- "SplineHazardRegression"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "SplineHazardRegression-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('SplineHazardRegression')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("etsim")
### * etsim

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: etsim
### Title: Simulate time-to-event data
### Aliases: etsim

### ** Examples

    #Generate input for a b-spline hazard and piecewise exponential survival censoring
    knots = c(0, 1, 3, 6, 10, NA, NA)
    betac = 1 * c(0.05, 0.05, 0.05, 0.05, 0.40, 0.1, 0.05)
    HParm = data.frame(knots, betac) # 'A Simple B-Spline'
    cll = c(0, 5)
    cup = c(5, 10)
    cih = c(0.0125, 0.025)
    CParm = data.frame(cll, cup, cih) # 'Light Censoring'
    INPUTS = etsim_inputs( HParam=HParm, CParam=CParm)

    #Then, we can generate the time-to event data
    SimDat = etsim(INPUTS)
    table( SimDat$status )




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("etsim", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("etsim_inputs")
### * etsim_inputs

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: etsim_inputs
### Title: Generate hazard and censoring distributions
### Aliases: etsim_inputs

### ** Examples

    #Generate input for a b-spline hazard simulation with piecewise exponential survival censoring
    knots = c(0, 1, 3, 6, 10, NaN, NaN )
    betac = 1 * c(0.05, 0.05, 0.05, 0.05, 0.40, 0.1, 0.05)
    HParm = data.frame(knots, betac) # 'A Simple B-Spline'
    cll = c(0, 5)
    cup = c(5, 10)
    cih = c(0.0125, 0.025)
    CParm = data.frame(cll, cup, cih) # 'Light Censoring'
    INPUTS = etsim_inputs( HParam=HParm, CParam=CParm )




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("etsim_inputs", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("generate_bspline_basis")
### * generate_bspline_basis

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: generate_bspline_basis
### Title: Generate B-Spline basis
### Aliases: generate_bspline_basis

### ** Examples

# A basis of cubic B-splines with no interior points
B = generate_bspline_basis( time = 0:10, Interior.knots=c(), Boundary.knots=c(0,10) )
B
# Note that B has 4=3+1 functions (as many as the order of cubic B-splines)
# we can plot the basis using matplot
matplot( 0:10, B, type="l")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("generate_bspline_basis", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
