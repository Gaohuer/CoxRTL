# CoxRTL

Deployment-oriented recalibrated transfer learning for survival prediction.

The package provides a reproducible simulation and evaluation workflow for
transfer-learning methods in survival analysis under covariate shift. It can

- simulate target, deployment and external survival datasets under covariate
  shift (`gen_surv_data()`),
- estimate density-ratio weights for covariate-shift correction
  (`calculate_weights()`),
- fit either of the two models per `fit_coxrtl()` call: a penalized Cox model
  trained on the target data only (`transfer = FALSE`) or a two-step
  transfer-learning model (`transfer = TRUE`) that borrows the external data.
  Both are recalibrated toward the deployment population when the density-ratio weights are supplied and are left uncalibrated otherwise,

The methodology is described in Pan et al. (2026) <doi:10.1038/s42256-026-01285-x>.

## Installation

```r
# install.packages("remotes")
remotes::install_github("Gaohuer/CoxRTL")
```

## Example

A fully worked example that mirrors the simulation study in the companion
paper is available in the package vignette:

```r
vignette("CoxRTL", package = "CoxRTL")
```

## Citation

Pan, L., Zhao, G., Yu, Y. et al. Transfer learning with deployment-covariate recalibration for survival prediction under covariate shift. Nat Mach Intell 8, 1312–1326 (2026). https://doi.org/10.1038/s42256-026-01285-x