# DTAtoolkit 0.1.0

Initial public beta with bivariate diagnostic-accuracy meta-analysis,
forest and SROC plots, subgroup comparison, influential-study analysis,
likelihood-ratio matrices, Fagan nomograms, bootstrap inference, and
publication-bias workflows.

- Combined forest plots now recover their natural content widths from the
  captured `meta::forest()` viewport when `meta::forest_dims()` is unavailable.
  This removes the unintended fixed 70/30 fallback, packs panels without a
  large internal gap, and iteratively fits narrow graphics devices.
- `left_panel_width_adjustment` is the preferred name for manual
  percentage-point changes to the automatic panel split. The former
  `sensitivity_width_adjustment` name remains available for compatibility.
