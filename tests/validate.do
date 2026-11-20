version 16
clear all
set more off
capture set processors 1
adopath ++ "."
set obs 84
gen id = floor((_n-1)/7)
gen time = mod(_n-1,7)+1
gen g = cond(id<4,3,cond(id<8,5,.))
gen double x = sin(id*2.3+time*.8)
gen x0 = mod(id,2)
gen double w = 1+id/20
gen double y = id*.3+time*.2+(time>=g)*(1+id*.2)+.3*x+.2*cos(id*time)
faststacked y, time(time) gvar(g) unit(id) window(-2 1) inference(jackknife)
assert r(native_se)>0
assert r(jk_se)>0
assert r(df)==11
assert abs(r(jk_se)-.302848859548824)<1e-9
matrix AA = J(3,12,.)
forvalues b=1/3 {
    forvalues c=1/12 {
        local original = mod(`c'-1-`b'+12,12)
        matrix AA[`b',`c'] = cond(`original'<4,3,cond(`original'<8,5,.))
    }
}
faststacked y, time(time) gvar(g) unit(id) window(-2 1) inference(ri_beta ri_t_native ) assignments(AA)
assert r(ri_beta_p)==.25
assert r(ri_t_native_p)==.25
matrix DUP = AA \ AA
capture faststacked y, time(time) gvar(g) unit(id) window(-2 1) inference(ri_beta) assignments(DUP)
assert _rc==459
assert _N==84
gen constant=1
capture faststacked y, time(time) gvar(g) unit(id) window(-2 1) covariates(constant)
assert _rc==498
assert _N==84
preserve
drop if id>0 & id<4
faststacked y, time(time) gvar(g) unit(id) window(-2 1) inference(jackknife)
assert r(jk_se)>0 & r(jk_se)<.
restore
replace y = id*.4+time*.3+2*(time>=g)
faststacked y, time(time) gvar(g) unit(id) window(-2 1)
assert abs(r(estimate)-2)<1e-9
gen allcluster=1
faststacked y, time(time) gvar(g) unit(id) cluster(allcluster) window(-2 1) inference(point)
assert abs(r(estimate)-2)<1e-9
assert missing(r(native_se))
assert missing(r(native_df))
capture faststacked y, time(time) gvar(g) unit(id) window(-2 1) inference(ri_beta_p)
assert _rc==198
display "Native Stata smoke assertions passed"
