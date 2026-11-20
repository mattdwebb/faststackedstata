version 16
clear all
set more off
set type double
capture set processors 1
adopath ++ "."
local oracle : environment STACKED_ORACLE
if "`oracle'" == "" {
    display as error "STACKED_ORACLE environment variable must point to frozen hollina/stacked"
    exit 198
}
adopath ++ "`oracle'/Stata/ado"
tempfile results
tempname handle
postfile `handle' str8 mode byte weighted adjusted double estimate double se double native_p double estimate_error double se_error double p_error using "`results'"
foreach mode in panel rcs {
    foreach w in FALSE TRUE {
        foreach x in FALSE TRUE {
            insheet using "tests/fixtures/`mode'_`w'_`x'.csv", clear comma names double
            capture confirm numeric variable g
            if _rc destring g, replace force
            local wo ""
            local xo ""
            if "`w'" == "TRUE" local wo "weightvar(w)"
            if "`x'" == "TRUE" local xo "covariates(x x0)"
            faststacked y, time(time) gvar(g) unit(id) window(-2 1) datamode(`mode') `wo' `xo'
            scalar ours_b = r(estimate)
            scalar ours_s = r(native_se)
            scalar ours_p = r(native_p)
            stacked build, time(time) unit(id) adopt(g) kpre(2) kpost(1) datatype(`mode') `wo'
            stacked reg y, cluster(id) `xo'
            scalar err_b = abs(ours_b-r(att))
            scalar err_s = abs(ours_s-r(att_se))
            quietly lincom (_sdtet_0+_sdtet_1)/2
            scalar err_p = abs(ours_p-r(p))
            assert err_b < 1e-9
            assert err_s < 1e-9
            assert err_p < 1e-12
            post `handle' ("`mode'") ("`w'"=="TRUE") ("`x'"=="TRUE") (ours_b) (ours_s) (ours_p) (err_b) (err_s) (err_p)
        }
    }
}
postclose `handle'
use "`results'", clear
format estimate se native_p estimate_error se_error p_error %24.17g
outsheet using "tests/oracle_results.csv", replace comma
summarize estimate_error se_error
display "STACKED_STATA_ORACLE_PASS: 8 scenarios x point and clustered SE"
