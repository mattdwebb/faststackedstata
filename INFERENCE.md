# Inference

## Point Only

inference='point' skips sandwich calculations and succeeds when the point
design is full rank even if there is only one original cluster. Native fields
are missing/None. Python/R deletion fits and Stata jackknife use this path.
Stata still uses native regress internally but does not construct clustered V.

## Native

CR1 uses the score summed over each ORIGINAL cluster across all stack copies.
B = (X'QX)^(-1), score_c = sum_{i in c} q_i x_i residual_i.
V = B sum_c(score_c score_c') B * G/(G-1)*(Nstack-1)/(Nstack-K).
Nstack counts copied rows; K is the full design rank. No singleton adjustment.
At least two used clusters and positive residual degrees of freedom are
required. The scalar variance is L'VL, not the mean of marginal variances.
native_df and df are G-1, and native_p uses a two-sided Student-t reference.
The frozen author's Stata _stacked_reg calls clustered regress followed by
lincom for the scalar average; this is Student t with G-1 df. The underlying
R fixest::feols coefficient table also uses t, and degrees_freedom(type='t')
was verified as 11 for the 12-cluster fixture (fixest 0.14.2). In contrast,
the separate R .compute_avg_post_att helper constructs CIs with qnorm and
does not return a p-value. Canonical native_p here follows author Stata scalar
lincom and the underlying R regression t-reference. It does NOT reproduce
the separate R helper's normal CI, and no such CI parity is claimed. Explicit
oracle tests validate the df, coefficient p-values, normal helper critical
value, and Stata scalar p-values. Zero native SE yields missing p.

## Jackknife

Delete every original input cluster and ALL its copies, rebuilding stacks,
eligibility, Q, cohort shares and regression each time. Variance is
(G-1)/G * sum((theta_minus_c - mean(theta_minus))^2), df=G-1.
G here counts input clusters, including any not used in observed stacks.
Every deletion must yield an estimable scalar aggregate. If a treated cohort
disappears, rebuild the eligible stacks, corrective weights and aggregate over
surviving supported effects. Losing a cohort alone is not a failure. Zero
variance is a failure.

## Supplied-Assignment RI

Python accepts a DataFrame with cluster-ID columns, one row per assignment.
R accepts a numeric matrix/data.frame with cluster-ID column names.
Stata accepts assignments(matrix_name), columns ordered by sorted original
cluster ID (numeric ascending or lexical string order); missing denotes never.

## November 2026 support revision

The Stata delete-cluster calculation now recomputes the estimator on surviving support. Disappearance of an adoption cohort by itself does not invalidate a deletion. Every required deletion must yield an estimate for a jackknife standard error.
