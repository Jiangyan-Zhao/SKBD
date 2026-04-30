
# SKBD

<!-- badges: start -->

[![R-CMD-check](https://github.com/Jiangyan-Zhao/SKBD/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Jiangyan-Zhao/SKBD/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`SKBD` implements the **Shared Keyboard Design (SKBD)** for
model-assisted phase I dose-finding trials.

The package provides practical tools for:

- generating SKBD decision tables,
- evaluating trial operating characteristics through simulation,
- extending SKBD to time-to-event settings with delayed toxicity
  outcomes, and
- assessing adaptive dose insertion under a shared-borrowing framework.

`SKBD` is designed for method development, simulation-based validation,
and practical trial planning in early-phase oncology studies.

## Installation

`SKBD` is not currently on CRAN. Install the development version from
GitHub:

``` r
# install.packages("pak")
pak::pak("Jiangyan-Zhao/SKBD")
```

Then load the package:

``` r
library(SKBD)
```

## Main functions

The current public interface includes five core functions:

- `get_boundary_SKBD()` generates pre-tabulated decision boundaries for
  the shared keyboard design.
- `get_OC_SKBD()` simulates operating characteristics for the standard
  SKBD.
- `get_OC_TITE_SKBD()` simulates operating characteristics for
  time-to-event SKBD.
- `get_OC_Insert_SKBD()` simulates operating characteristics when
  adaptive dose insertion is enabled.
- `PUA()` generates monotone dose-toxicity scenarios using a
  pseudo-uniform algorithm.

## Shiny app

A Shiny app is included for interactive use of SKBD decision tables and
simulation outputs.

### What you can do in the app

- Configure trial design settings, including dose levels, target
  toxicity probability, borrowing strength, sample size, and
  overdose-control parameters.
- Input current trial observations and generate SKBD decision tables
  interactively.
- Run multi-scenario operating-characteristic simulations using manual
  input or uploaded CSV templates.

### Launch the Shiny app

After installing **SKBD**, load the package and launch the app by
running:

``` r
library(SKBD)

run_SKBD_shiny()
```

By default, `run_SKBD_shiny()` follows the standard launching behavior
of `shiny::runApp()`. In RStudio, the app is typically opened in the
RStudio Viewer or Shiny window.

To open the app in the system default browser, use:

``` r
run_SKBD_shiny(launch.browser = TRUE)
```

To start the app without automatically opening a browser, use:

``` r
run_SKBD_shiny(launch.browser = FALSE)
```

The app provides two main modules:

- **Trial Setting**: interactively generates dose-escalation and
  de-escalation boundary tables.
- **Simulation**: runs batch simulations under user-specified scenarios
  and summarizes the operating characteristics.

## Decision tables for SKBD

A typical starting point is to generate a decision table at the current
dose level, conditional on the observed data across all dose levels.

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

The output is a keyboard-style table with escalation, de-escalation, and
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

The returned object includes accuracy and safety summaries, such as the
percentage of correct selection (`PCS`), percentage of correct
allocation (`PCA`), and overdose risk.

## Time-to-event SKBD

Delayed toxicity outcomes can be handled using time-to-event SKBD:

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

This extension accounts for pending toxicity outcomes by incorporating
weighted follow-up information within the DLT assessment window.

## Adaptive dose insertion

`SKBD` also supports simulations for insertion-enabled shared keyboard
designs:

``` r
out_insert <- get_OC_Insert_SKBD(
  target_prob = 0.30,
  tox_prob = c(0.14, 0.45, 0.63, 0.74, 0.80),
  dose_set = c(5, 15, 25, 35, 45),
  n_trial = 1000
)

out_insert$insertion
#> $sel_pct
#> [1] 80.7
#> 
#> $pts_pct
#> [1] 49.38
#> 
#> $dose_mean
#> [1] 10.64145
#> 
#> $dose_sd
#> [1] 3.651881
#> 
#> $trial_pct
#> [1] 95.3
#> 
#> $cohort_mean
#> [1] 3.373277
#> 
#> $n_median
#> [1] 3
```

The insertion summary reports how often new doses are inserted, where
inserted doses tend to be located, how often they are selected as the
final MTD, and how patients are allocated to inserted doses.

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

`SKBD` is under active development. The current version focuses on
decision-table utilities, simulation engines, time-to-event extensions,
and adaptive dose insertion for the SKBD framework. Additional examples,
validation materials, and extended tutorials will be added as the
package evolves.
