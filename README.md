# DTAtoolkit

`DTAtoolkit` provides streamlined R workflows for meta-analysis of diagnostic
test accuracy from study-level 2 x 2 tables. It supports unified sensitivity
and specificity forest plots, SROC curves, subgroup comparisons,
likelihood-ratio matrices, Fagan nomograms, influential-study analysis,
bootstrap inference, and publication-bias assessment.

The main diagnostic-accuracy analyses are based on the bivariate Reitsma model
for jointly modelling logit sensitivity and specificity
([Reitsma et al., 2005](https://doi.org/10.1016/j.jclinepi.2005.02.022)). SROC
AUC uncertainty can be estimated by bootstrap methods informed by the
`dmetatools` approach
([Freeman et al., 2021](https://doi.org/10.1080/23737484.2021.1894408)), and
heterogeneity summaries include the bivariate I-squared approach described by
[Holling et al.](https://doi.org/10.1080/03610918.2018.1489553). Selected
univariate diagnostic and influence procedures use methods provided by
[`metafor`](https://doi.org/10.18637/jss.v036.i03). Publication-bias assessment
includes the generalized Egger regression framework described by
[Noma (2020)](https://doi.org/10.1111/biom.13343) and
[Noma (2022)](https://doi.org/10.48550/arXiv.2209.07270).

The Fagan plot and likelihood-ratio matrix presentation were inspired by the
[MIDAS Stata module](https://ideas.repec.org/c/boc/bocode/s456880.html). Pooled
likelihood-ratio results in `DTAtoolkit` are derived from the fitted bivariate
model rather than by separately pooling diagnostic likelihood ratios, in line
with the methodological concerns discussed by
[Zwinderman and Bossuyt (2008)](https://doi.org/10.1002/sim.2992).

This package is based on a personal project available at
[Payam-PJD/My-mada-functions](https://github.com/Payam-PJD/My-mada-functions).

## Installation

To install a downloaded source tarball, run the following from the directory
containing the file:

```r
install.packages(
  "DTAtoolkit_0.1.0.tar.gz",
  repos = NULL,
  type = "source",
  dependencies = TRUE
)
```

To install the current GitHub version:

```r
install.packages("remotes")
remotes::install_github("Payam-PJD/DTAtoolkit", dependencies = TRUE)
```

Load the package with:

```r
library(DTAtoolkit)
```


## Example data

`dta_example` is a simulated dataset containing study identifiers, a subgroup
variable, and true-positive, true-negative, false-positive, and false-negative
counts.

```r
data(dta_example)
head(dta_example)
```

## Examples

Create a combined sensitivity and specificity forest plot:

```r
dta_example |>
  dta_forest(
    study = study,
    tp = tp, tn = tn, fp = fp, fn = fn,
    type = "combined"
  )
```

Add subgroup estimates and the bivariate subgroup interaction test:

```r
dta_example |>
  dta_forest(
    subgroup = subgroup,
    study = study,
    tp = tp, tn = tn, fp = fp, fn = fn,
    type = "combined",
    pairwise = "never"
  )
```

Draw subgroup SROC curves:

```r
dta_example |>
  dta_sroc(
    subgroup = subgroup,
    study = study,
    tp = tp, tn = tn, fp = fp, fn = fn,
    pairwise = "never"
  )
```

Request bootstrap confidence intervals for the SROC AUC:

```r
dta_example |>
  dta_sroc(
    subgroup = subgroup,
    study = study,
    tp = tp, tn = tn, fp = fp, fn = fn,
    auc_ci = TRUE,
    n_boots = 2000,
    seed = 2026,
    pairwise = "never"
  )
```

Create a likelihood-ratio matrix:

```r
dta_example |>
  dta_lr_matrix(
    subgroup = subgroup,
    study = study,
    tp = tp, tn = tn, fp = fp, fn = fn
  )
```

Run within-subgroup influential-study analysis and remove flagged rows:

```r
reduced_data <- dta_example |>
  dta_influentials(
    subgroup = subgroup,
    study = study,
    tp = tp, tn = tn, fp = fp, fn = fn
  ) |>
  remove_inf()
```

## To do

1. Add multiple-cutoff support.
2. Add a manual and vignette.
3. Optimize nomogram plots.
4. Optimize heterogeneity reporting.
