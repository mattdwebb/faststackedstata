# Validation (2026-09-29)

Status: tested reference implementation, PARTIAL against the full suite brief.
Frozen external author: hollina/stacked commit
3116918a5773c6415acc9bce6699c5c15b3a1792 (MIT).

Stata oracle: eight corresponding scenarios against frozen author Stata code;
maximum point difference 2.89e-15 and CR1 SE difference 6.22e-15. StataNow 19.5
SE, native regress backend (not reghdfe). Standalone smoke covers JK, all three
RI methods, homogeneous-effect recovery and one-cluster point-only fitting.

Python: 70 pytest cases pass, including eight frozen external-oracle scenarios,
hand-derived Q weights, full deletion rebuilding, supplied-assignment inference,
strict failures, outcome linearity, two-worker spawn isolation and 42 R pt
probability cases (both signs). NumPy 2.5.3, pandas 3.0.6, Python 3.14.
R two-worker process isolation and independent package install/build pass.
R CMD check --no-manual --no-vignettes finishes Status: OK, with zero
warnings/notes, using LC_ALL=C and the optional oracle dependency library.
Across all three native languages and eight scenarios, max point disagreement
is 1.16e-14, native SE 1.31e-14, and native p-value 1.30e-17.
The reproducible development check is faststackedpy/tests/check_cross_language.py.

The R author scalar helper uses normal CI critical values; native Stata lincom
uses Student t. fixest::degrees_freedom(fit,type='t') was verified as G-1 (11)
on the clustered fixture. Our canonical native_p follows author Stata scalar
lincom and the underlying clustered fixest t-reference, NOT the separate
R avg_post_att helper's normal CI. No parity claim for that helper's CI.

No full acceptance certification: automatic enumeration/sampling is not yet
provided by these packages, broad design and large-N benchmarks remain pending.
Tests are deterministic synthetic data, not a Monte Carlo. Python and R
Windows/Linux CI passed on the final implementation/audit commits. Stata's
licensed self-hosted workflow is configured; Stata checks ran locally.
CI evidence: Python run 36523668906; R run 36523678300.
Raw Stata startup/license logs are ignored
and never committed.
