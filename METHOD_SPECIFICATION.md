# Corrective-Weighted Stacked DiD

## Target and Stack Construction

This implements the saturated WFH event-study regression, as in equations
(1) and the corrective weighting results of Wing, Freedman and Hollingsworth,
NBER Working Paper 32054 (2024), and the default saturated event-study path
of hollina/stacked at commit 3116918a5773c6415acc9bce6699c5c15b3a1792.

Integer calendar periods and absorbing treatment are required. Never-treated
cohorts are missing or positive infinity (Stata: missing). Zero is a real date.
The full inclusive event window is (-pre, post), with pre >= 1 and post >= 0.
The reference is -1. No arbitrary window is selected. A harness post-only
window (0,H) must be explicitly mapped to (-1,H). Nonzero post starts are
not implemented. All post effects from 0 through H receive equal weight.

For cohort a, include a-pre through a+post, treated cohort g=a, and clean
controls g>a+post or never treated. control_group can restrict controls to
never_treated or not_yet_treated. Trim cohorts whose full window is outside
observed calendar support; drop cohorts with no controls. Report used and
dropped cohorts. Fail on empty event-by-treatment cells. Panel mode requires
unique unit-time rows, a time-invariant cohort and cluster per unit, and every
included unit observed throughout each stack. Unbalanced panels fail rather
than silently receive panel weights.

## Q Weights

Let w be the positive supplied observation weight (one by default). Let T_a
and C_a be sums of w over treated and control rows in stack a, including all
event times. Let T and C sum these quantities across stacks, counting repeated
control copies. Target cohort shares are pi_a = T_a/T.

For balanced panels, treated Q = w; control Q = w*pi_a/(C_a/C).
Supplied panel weights must be constant within unit across each stack.
With a common window, pi_a is the treated-unit share (weighted if requested).

For repeated cross sections, let T_ae and C_ae sum w in each stack/event cell,
and T_e and C_e sum across stacks at event e. Treated Q = w*pi_a/(T_ae/T_e);
control Q = w*pi_a/(C_ae/C_e). This is the frozen official RCS extension.
It fixes the relative cohort composition in each event/treatment cell, without
asserting that sampled people are panel units. Here pi_a is a treated sampled
observation-weight share over the full window, not automatically an
administrative-unit population share. Survey-design variance beyond one-way
clustering is not implemented.

These weights equalize cohort shares in treated and control event cells.
Without covariates, weighted event-time means therefore combine within-cohort
DID contrasts with identical pi_a on both sides; naive stacking generally
uses different treated and control shares and does not recover this target.

## Regression and Aggregation

WLS regressors are an intercept, treated status within stack, nonreference
event indicators, and their treated interactions. No stack fixed effects
are added. Scalar estimate is mean of the post interaction coefficients.
Event covariance and the full cross-event covariance contribution to the
scalar standard error are retained. Baseline event -1 is exactly zero.

## Covariates and Interpretation

Requested numeric covariates enter additively at their contemporaneous row
values. Time-varying values are not frozen, differenced, or lagged; time-
invariant covariates also work if full rank. This matches the official
stackreg covariates option but is an additional regression-adjustment
specification, not a covariate-balancing or propensity-score estimator.
Arbitrary post-treatment covariates may induce bias. Causal interpretation
requires defensible covariate exogeneity and a suitable common additive
outcome model in addition to parallel trends/no anticipation. The adjusted
output is labelled a regression-adjusted Q-weighted post-event contrast,
not unconditionally asserted to equal the unadjusted trimmed ATT.

## Linearity

For fixed design and weights, theta = L'(X'QX)^(-1)X'QY. The stack duplicates
outcomes through a fixed selection matrix S, so theta = L'(X'QX)^(-1)X'QS y
is exactly linear in original outcomes. Covariates are not fit using y.
Native and jackknife variances and studentized statistics are not linear.
Assignments and deletions change S, Q, X and L, and require full rebuilding.

## Explicit Limitations

No interacted FE, propensity-score balancing, population-weight variants,
multiple-way clustering, missing-data deletion, nonabsorbing treatment,
automatic window selection, or large-data acceleration. Unknown options fail.
The package is not a complete port of every author-package extension.

