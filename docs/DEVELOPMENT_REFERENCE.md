> Historical development notes retained for provenance. For current behavior, installation and examples, use the [main user guide](../README.md). Old audit statements below may describe superseded behavior.

# faststackedstata

Native Stata/Mata corrective-weighted WFH saturated event study. MIT, 2026 Matt Webb.
No R/Python runtime calls and no user-written runtime dependencies.

From this directory, `net install faststacked, from(.) replace`, or use
`adopath ++ "path/to/faststackedstata"`. Stata 16+ required.

```stata
faststacked y, time(time) gvar(g) unit(id) cluster(id) window(-2 1) inference(native)
display r(estimate)
```

`do tests/validate.do` is a standalone synthetic runnable example and smoke test.
`do tests/oracle.do` compares frozen author Stata results when STACKED_ORACLE
points to hollina/stacked commit 3116918a5773c6415acc9bce6699c5c15b3a1792.
Tests use native insheet, not Java-dependent import delimited.