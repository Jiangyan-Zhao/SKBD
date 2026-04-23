
# SKBD

`SKBD` implements the shared keyboard design (SKBD) for model-assisted
phase I dose-finding trials. The package provides tools for

- constructing decision tables for the shared keyboard design,
- simulating operating characteristics under the standard SKBD,
- extending SKBD to time-to-event settings with late-onset toxicity, and
- evaluating adaptive dose insertion under a shared-keyboard framework.

The package is intended for method development, operating-characteristic
evaluation, and practical exploration of shared keyboard decision rules
in phase I oncology trials.

## Installation

`SKBD` is not currently on CRAN. You can install the development version
from GitHub:

``` r
install.packages("remotes")
remotes::install_github("Jiangyan-Zhao/SKBD")
```

Then load the package:

``` r
library(SKBD)
```

## Main functions

The current public interface includes five main functions:

- `get_boundary_SKBD()` generates pre-tabulated decision boundaries for
  the shared keyboard design.
- `get_OC_SKBD()` simulates operating characteristics for the standard
  SKBD.
- `get_OC_TITE_SKBD()` simulates operating characteristics for the
  time-to-event shared keyboard design.
- `get_OC_Insert_SKBD()` simulates operating characteristics when
  adaptive dose insertion is allowed.
- `PUA()` generates random monotone dose-toxicity scenarios using a
  pseudo-uniform algorithm.

## Decision tables for SKBD

A typical starting point is to generate a decision table at a current
dose level, conditional on the observed data across doses.

``` r
y <- c(0, 1, 2, 2, 0)
n <- c(3, 6, 9, 3, 0)

out_boundary <- get_boundary_SKBD(
  target_prob = 0.30,
  d = 3,
  y = y,
  n = n,
  table_type = "continue"
)

out_boundary$boundary_tab
#>                                                         
#> Number of patients treated  3  6  9 12 15 18 21 24 27 30
#> Escalate if # of DLT <=    NA NA NA  2  2  3  4  4  5  6
#> de-escalate if # of DLT >= NA NA NA  4  5  6  7  8  9 10
#> Eliminate if # of DLT >=   NA NA NA NA NA 10 11 13 14 15
```

This returns a keyboard-style table with escalation, de-escalation, and
elimination boundaries based on the SKBD pseudo-posterior at the current
dose.

## Simulating operating characteristics under SKBD

The standard SKBD can be evaluated under a prespecified toxicity
scenario as follows:

``` r
out_skbd <- get_OC_SKBD(
  target_prob = 0.30,
  tox_prob = c(0.05, 0.12, 0.30, 0.45, 0.60),
  n_cohort = 10,
  cohort_size = 3,
  n_trial = 1000
)

out_skbd$PCS
#> [1] 61.8
out_skbd$PCA
#> [1] 37.49
out_skbd$ROD60
#> [1] 1.8
```

The returned object includes accuracy and safety summaries such as
percent correct selection (`PCS`), percent correct allocation (`PCA`),
and overdosing risk.

## Time-to-event SKBD

Late-onset toxicity can be handled through the time-to-event shared
keyboard design:

``` r
out_tite <- get_OC_TITE_SKBD(
  target_prob = 0.20,
  tox_prob = c(0.05, 0.12, 0.20, 0.35, 0.50),
  n_cohort = 10,
  cohort_size = 3,
  tau = 3,
  accrual = 2,
  dist_DLT = "weibull",
  dist_enter = "exp",
  n_trial = 1000
)

out_tite$PCS
#> [1] 44.3
out_tite$duration_mean
#> [1] 23.56665
```

This version accounts for pending toxicity outcomes through weighted
follow-up within the DLT assessment window.

## Adaptive dose insertion

The package also supports simulations under an insertion-enabled shared
keyboard design:

``` r
out_insert <- get_OC_Insert_SKBD(
  target_prob = 0.30,
  tox_prob = c(0.14, 0.45, 0.63, 0.74, 0.80),
  dose_set = c(5, 15, 25, 35, 45),
  n_trial = 1000
)

out_insert$insertion
```

The insertion summary reports how often new doses are inserted, where
inserted doses tend to be selected, and how patients are allocated to
inserted doses.

## Random monotone scenarios

For random-scenario simulation studies, `PUA()` can be used to generate
monotone dose-toxicity curves with a well-defined target dose.

``` r
scen <- PUA(
  dose_set = 1:5,
  target_prob = 0.30,
  n_scenarios = 5
)

scen
#>            [,1]       [,2]      [,3]      [,4]      [,5]
#> [1,] 0.01944263 0.02871687 0.1536818 0.1940659 0.2761874
#> [2,] 0.04828428 0.25536821 0.3181795 0.5105885 0.5796248
#> [3,] 0.19502259 0.27180464 0.5551730 0.6691403 0.7765835
#> [4,] 0.17182210 0.33929376 0.4044292 0.6745865 0.8366659
#> [5,] 0.09726789 0.10416956 0.1283349 0.1508885 0.2799666
```

## Development notes

`SKBD` is under active development. The package currently focuses on
core simulation and decision-table tools for the shared keyboard
framework. Additional examples, validation materials, and extended
documentation can be added as the package evolves.

## Interactive Shiny app

You can launch an interactive Shiny interface to explore SKBD decision
boundaries and run basic operating-characteristic simulations:

```r
run_SKBD_shiny()
```

If needed, install Shiny first:

```r
install.packages("shiny")
```
