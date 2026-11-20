{smcl}
{title:faststacked - corrective-weighted saturated stacked DiD}
{p 4 4 2}
Development reference. See VALIDATION.md for tested scope and limitations.
{p 4 4 2}
faststacked outcome, time(year) gvar(cohort) unit(id) cluster(state)
window(-2 1) [datamode(panel|rcs|auto) covariates(x1 x2)
weightvar(w) controlgroup(both|never_treated|not_yet_treated)
inference(point|native|jackknife|ri_beta|ri_t_native||)
assignments(matrix_name)]
{p 4 4 2}
Missing cohort means never treated. Zero is a real date. Window contains the
entire event study, including baseline -1. Integer time only. Panel must be
balanced in each stack. RCS needs cohort input and original cluster IDs.
Additive covariates enter contemporaneously without transformation and require
additional identifying assumptions; this is not covariate balancing.
{p 4 4 2}
Returns r(estimate), r(native_se), r(native_p), r(native_df), r(df), r(n_stacked),
r(n_clusters), r(event_study), r(b), r(V), r(cohorts), r(target_label).
Jackknife also returns r(jk_se), r(jk_p), r(jk_df), r(jackknife_estimates).
RI returns method_p scalars and r(randomized_statistics). Supply unique,
cohort-count-preserving alternative assignments, columns in sorted original
cluster-ID order, missing for never treated. RI methods may be combined as a
space-separated list. Point must be used alone. Automatic assignment generation,
treatment-only input, interacted fixed effects, population-weight
extensions, and missing-data deletion are unsupported and fail explicitly.
