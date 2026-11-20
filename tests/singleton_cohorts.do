version 16
clear all
set more off
adopath ++ "."
forvalues cv=0/1 {
    insheet using "tests/singleton_reference.csv", comma names clear double
    keep if adjusted==`cv'
    mkmat estimate, matrix(expected)
    scalar expected_se=jk_se[1]
    insheet using "tests/singleton_cohorts.csv", comma names clear double
    local controls
    if `cv' local controls covariates(x)
    faststacked y, time(time) cluster(state) gvar(g) datamode(rcs) window(-2 1) weightvar(w) `controls' inference(jackknife)
    matrix deletes=r(jackknife_estimates)
    assert abs(r(estimate)-expected[1,1])<1e-9
    assert abs(r(jk_se)-expected_se)<1e-9
    assert r(jk_df)==7
    forvalues i=1/8 {
        assert abs(deletes[`i',1]-expected[`i'+1,1])<1e-9
        preserve
        drop if state==`i'
        faststacked y, time(time) cluster(state) gvar(g) datamode(rcs) window(-2 1) weightvar(w) `controls' inference(point)
        assert abs(r(estimate)-deletes[`i',1])<1e-9
        restore
    }
}
display "SINGLETON_COHORT_STATA_PYTHON_PARITY_PASS"
