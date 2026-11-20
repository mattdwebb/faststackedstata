*! faststacked 0.1.0: corrective-weighted saturated reference implementation
program define faststacked, rclass
    version 16
    syntax varname(numeric), TIME(varname numeric) GVAR(varname numeric) ///
        WINDOW(numlist integer min=2 max=2) [UNIT(varname) CLuster(varname) ///
        COVariates(varlist numeric) WEIGHTvar(varname numeric) ///
        DATAMode(string) CONTROLGroup(string) INference(string) ASSIGNments(name)]
    if "`inference'" == "" local inference "native"
    local wantjk = 0
    local wantri = 0
    foreach method of local inference {
        if !inlist("`method'", "point", "native", "jackknife", "ri_beta", "ri_t_native") {
            di as error "UNSUPPORTED inference method: `method'"
            exit 198
        }
        if inlist("`method'", "jackknife") local wantjk = 1
        if inlist("`method'", "ri_beta", "ri_t_native") local wantri = 1
    }
    if strpos(" `inference' "," point ") & "`inference'" != "point" exit 198
    if `wantri' & "`assignments'" == "" {
        di as error "RI requires assignments(matrix_name), columns in sorted original-cluster order"
        exit 198
    }
    if !`wantri' & "`assignments'" != "" exit 198
    if "`cluster'" == "" local cluster "`unit'"
    if "`cluster'" == "" {
        di as error "An original cluster identifier is required"
        exit 198
    }
    local opts "time(`time') gvar(`gvar') window(`window') cluster(`cluster')"
    if "`unit'" != "" local opts "`opts' unit(`unit')"
    if "`covariates'" != "" local opts "`opts' covariates(`covariates')"
    if "`weightvar'" != "" local opts "`opts' weightvar(`weightvar')"
    if "`datamode'" != "" local opts "`opts' datamode(`datamode')"
    if "`controlgroup'" != "" local opts "`opts' controlgroup(`controlgroup')"
    local pointopt ""
    if "`inference'" == "point" local pointopt "point"
    quietly _faststacked_point `varlist', `opts' `pointopt'
    tempname est se p df ns nc bb vv es jk jkse jkp
    scalar `est' = r(estimate)
    scalar `se' = r(native_se)
    scalar `p' = r(native_p)
    scalar `df' = r(df)
    scalar `ns' = r(n_stacked)
    scalar `nc' = r(n_clusters)
    matrix `bb' = r(b)
    matrix `vv' = r(V)
    matrix `es' = r(event_study)
    local cohorts "`r(cohorts)'"
    local target "`r(target_label)'"
    if `wantjk' {
        preserve
        tempfile base
        tempvar cid
        quietly egen long `cid' = group(`cluster')
        quietly save "`base'"
        quietly levelsof `cid', local(ids)
        local G : word count `ids'
        matrix `jk' = J(`G', 1, .)
        local j = 0
        foreach c of local ids {
            quietly use "`base'", clear
            quietly drop if `cid' == `c'
            quietly _faststacked_point `varlist', `opts' point
            local ++j
            matrix `jk'[`j',1] = r(estimate)
        }
        restore
        mata: st_numscalar("`jkse'", sqrt((rows(st_matrix("`jk'"))-1)/rows(st_matrix("`jk'"))*sum((st_matrix("`jk'"):-mean(st_matrix("`jk'"))):^2)))
        if scalar(`jkse') <= 1e-14 {
            di as error "IDENTIFICATION_FAILURE: zero jackknife variance"
            exit 498
        }
        scalar `jkp' = 2*ttail(`G'-1,abs(scalar(`est')/scalar(`jkse')))
    }
    if `wantri' {
        confirm matrix `assignments'
        preserve
        tempvar rid
        quietly egen long `rid' = group(`cluster')
        bysort `rid': assert `gvar' == `gvar'[1]
        tempfile original
        quietly save "`original'"
        quietly by `rid': keep if _n == 1
        quietly sort `rid'
        local J = _N
        tempname observed valid stats raw tse
        quietly mkmat `gvar', matrix(`observed')
        matrix `observed' = `observed''
        local B = rowsof(`assignments')
        if colsof(`assignments') != `J' | `B' < 1 {
            di as error "Assignment matrix dimensions do not match original clusters"
            exit 503
        }
        matrix `stats' = J(`B',2,.)
        matrix colnames `stats' = ri_beta ri_t_native
        forvalues b = 1/`B' {
            mata: st_numscalar("`valid'", all(sort(st_matrix("`assignments'")[`b',.]',1):==sort(st_matrix("`observed'")',1)) & !all(st_matrix("`assignments'")[`b',.]:==st_matrix("`observed'")))
            if !scalar(`valid') {
                di as error "Assignments must preserve cohort counts and exclude observed assignment"
                exit 459
            }
            if `b' > 1 {
                forvalues h = 1/`=`b'-1' {
                    mata: st_numscalar("`valid'", !all(st_matrix("`assignments'")[`b',.]:==st_matrix("`assignments'")[`h',.]))
                    if !scalar(`valid') {
                        di as error "Duplicate supplied assignments are not allowed"
                        exit 459
                    }
                }
            }
            quietly use "`original'", clear
            forvalues j = 1/`J' {
                quietly replace `gvar' = `assignments'[`b',`j'] if `rid' == `j'
            }
            local submethod "native"
            quietly faststacked `varlist', `opts' inference(`submethod')
            if "`r(cohorts)'" != "`cohorts'" {
                di as error "IDENTIFICATION_FAILURE: assignment changes cohort support"
                exit 498
            }
            scalar `raw' = r(estimate)
            scalar `tse' = r(native_se)
            matrix `stats'[`b',1] = scalar(`raw')
            if strpos(" `inference' "," ri_t_native ") {
                if scalar(`tse') <= 1e-14 | missing(scalar(`tse')) exit 498
                matrix `stats'[`b',2] = scalar(`raw')/scalar(`tse')
            }
        }
        restore
        tempname statistic rp
        foreach method of local inference {
            local col = 0
            if "`method'" == "ri_beta" {
                local col = 1
                scalar `statistic' = scalar(`est')
            }
            if "`method'" == "ri_t_native" {
                if scalar(`se') <= 1e-14 | missing(scalar(`se')) exit 498
                local col = 2
                scalar `statistic' = scalar(`est')/scalar(`se')
            }
            if `col' {
                mata: st_numscalar("`rp'", (1+sum(abs(st_matrix("`stats'")[.,`col']):>=abs(st_numscalar("`statistic'")):-64*2.220446049250313e-16:*rowmax((J(rows(st_matrix("`stats'")),1,1),abs(st_matrix("`stats'")[.,`col']),J(rows(st_matrix("`stats'")),1,abs(st_numscalar("`statistic'")))))))/(rows(st_matrix("`stats'"))+1))
                return scalar `method'_p = scalar(`rp')
            }
        }
        return matrix randomized_statistics = `stats'
        return scalar ri_reps = `B'
    }
    if `wantjk' {
        return scalar jk_se = scalar(`jkse')
        return scalar jk_p = scalar(`jkp')
        return scalar jk_df = `G'-1
        return matrix jackknife_estimates = `jk'
    }
    return scalar estimate = scalar(`est')
    return scalar native_se = scalar(`se')
    return scalar native_p = scalar(`p')
    return scalar df = scalar(`df')
    return scalar native_df = scalar(`df')
    return scalar n_stacked = scalar(`ns')
    return scalar n_clusters = scalar(`nc')
    return matrix b = `bb'
    return matrix V = `vv'
    return matrix event_study = `es'
    return local cohorts "`cohorts'"
    return local target_label "`target'"
    return local native_reference "Student t, G-1: author Stata lincom; R scalar helper CI uses normal"
    return local status "PASS"
    di as text "Q-weighted post-event estimate = " as result %12.7g scalar(`est') ///
        as text "  CR1 SE = " as result %12.7g scalar(`se')
end

program define _faststacked_point, rclass
    version 16
    syntax varname(numeric), TIME(varname numeric) GVAR(varname numeric) ///
        WINDOW(numlist integer min=2 max=2) CLuster(varname) [UNIT(varname) ///
        COVariates(varlist numeric) WEIGHTvar(varname numeric) DATAMode(string) CONTROLGroup(string) POINT]
    local lo : word 1 of `window'
    local hi : word 2 of `window'
    if `lo' > -1 | `hi' < 0 {
        di as error "window(-pre post) must contain baseline -1 and event 0"
        exit 198
    }
    if inlist("`datamode'", "", "auto") {
        if "`unit'" == "" local datamode "rcs"
        else local datamode "panel"
    }
    if !inlist("`datamode'", "panel", "rcs") exit 198
    if "`datamode'" == "panel" & "`unit'" == "" exit 198
    if "`controlgroup'" == "" local controlgroup "both"
    if !inlist("`controlgroup'", "both", "never_treated", "not_yet_treated") exit 198
    foreach v of varlist `varlist' `time' `cluster' `unit' `covariates' `weightvar' {
        quietly count if missing(`v')
        if r(N) {
            di as error "Missing required values are not silently dropped"
            exit 459
        }
    }
    assert `time' == floor(`time')
    assert `gvar' == floor(`gvar') if !missing(`gvar')
    preserve
    tempvar row w a ev z q tr co ta ca tt ct tse cse ec et cc clean
    quietly gen long `row' = _n
    quietly gen double `w' = 1
    if "`weightvar'" != "" quietly replace `w' = `weightvar'
    assert `w' > 0
    if "`datamode'" == "panel" {
        isid `unit' `time'
        bysort `unit': assert `gvar' == `gvar'[1]
        bysort `unit': assert `cluster' == `cluster'[1]
    }
    quietly summarize `time', meanonly
    local mint = r(min)
    local maxt = r(max)
    quietly levelsof `gvar' if !missing(`gvar'), local(allcohorts)
    tempfile original stacked
    quietly save "`original'"
    local first = 1
    local cohorts ""
    foreach focal of local allcohorts {
        if `focal'+`lo' < `mint' | `focal'+`hi' > `maxt' continue
        quietly use "`original'", clear
        quietly gen byte `clean' = missing(`gvar') | `gvar' > `focal'+`hi'
        if "`controlgroup'" == "never_treated" quietly replace `clean' = missing(`gvar')
        if "`controlgroup'" == "not_yet_treated" quietly replace `clean' = !missing(`gvar') & `gvar' > `focal'+`hi'
        quietly keep if (`gvar' == `focal' | `clean') & inrange(`time',`focal'+`lo',`focal'+`hi')
        quietly gen double `a' = `focal'
        quietly gen int `ev' = `time'-`focal'
        quietly gen byte `z' = `gvar' == `focal'
        quietly summarize `z', meanonly
        if r(min) == r(max) | r(N) == 0 continue
        local supported = 1
        forvalues e = `lo'/`hi' {
            forvalues d = 0/1 {
                quietly count if `ev' == `e' & `z' == `d'
                if r(N) == 0 local supported = 0
            }
        }
        if !`supported' continue
        if "`datamode'" == "panel" {
            bysort `unit': assert _N == `hi'-`lo'+1
            bysort `unit': assert `w' == `w'[1]
        }
        local cohorts "`cohorts' `focal'"
        if !`first' quietly append using "`stacked'"
        quietly save "`stacked'", replace
        local first = 0
    }
    if `first' {
        di as error "IDENTIFICATION_FAILURE: no identified cohorts in window"
        exit 498
    }
    quietly use "`stacked'", clear
    quietly gen double `tr' = `w'*`z'
    quietly gen double `co' = `w'*(1-`z')
    quietly egen double `ta' = total(`tr'), by(`a')
    quietly egen double `ca' = total(`co'), by(`a')
    quietly egen double `tt' = total(`tr')
    quietly egen double `ct' = total(`co')
    quietly gen double `q' = `w'
    if "`datamode'" == "panel" {
        quietly replace `q' = `w'*(`ta'/`tt')/(`ca'/`ct') if `z' == 0
    }
    else {
        quietly egen double `tse' = total(`tr'), by(`a' `ev')
        quietly egen double `cse' = total(`co'), by(`a' `ev')
        quietly egen double `et' = total(`tr'), by(`ev')
        quietly egen double `ec' = total(`co'), by(`ev')
        quietly replace `q' = `w'*(`ta'/`tt')/(`tse'/`et') if `z' == 1
        quietly replace `q' = `w'*(`ta'/`tt')/(`cse'/`ec') if `z' == 0
    }
    local rhs "`z'"
    local expr ""
    local interactions ""
    forvalues e = `lo'/`hi' {
        if `e' == -1 continue
        tempvar dd dz
        quietly gen byte `dd' = `ev' == `e'
        quietly gen byte `dz' = `dd'*`z'
        local rhs "`rhs' `dd' `dz'"
        local interactions "`interactions' `dz'"
        if `e' >= 0 {
            if "`expr'" == "" local expr "`dz'"
            else local expr "`expr' + `dz'"
        }
    }
    quietly egen long `cc' = group(`cluster')
    local k : word count `rhs' `covariates'
    tempname designrank
    mata: st_numscalar("`designrank'", rank((J(st_nobs(),1,1),st_data(.,tokens("`rhs' `covariates'")))))
    if scalar(`designrank') != `k'+1 {
        di as error "IDENTIFICATION_FAILURE: singular regression design"
        exit 498
    }
    quietly summarize `cc', meanonly
    local G = r(max)
    if "`point'" == "" & `G' < 2 {
        di as error "IDENTIFICATION_FAILURE: singular regression or too few clusters"
        exit 498
    }
    if "`point'" == "" quietly regress `varlist' `rhs' `covariates' [aw=`q'], vce(cluster `cc')
    else quietly regress `varlist' `rhs' `covariates' [aw=`q']
    tempname bb vv es estimate se p
    matrix `bb' = e(b)
    matrix `vv' = e(V)
    local N = e(N)
    matrix `es' = J(`hi'-`lo'+1,3,0)
    matrix colnames `es' = event_time estimate native_se
    local j = 0
    forvalues e = `lo'/`hi' {
        local pos = `e'-`lo'+1
        matrix `es'[`pos',1] = `e'
        if `e' == -1 continue
        local ++j
        local term : word `j' of `interactions'
        matrix `es'[`pos',2] = _b[`term']
        if "`point'" == "" matrix `es'[`pos',3] = _se[`term']
        else matrix `es'[`pos',3] = .
    }
    quietly lincom (`expr')/(`hi'+1)
    scalar `estimate' = r(estimate)
    scalar `se' = r(se)
    scalar `p' = r(p)
    local df = `G'-1
    if "`point'" != "" {
        scalar `se' = .
        scalar `p' = .
        local df = .
        matrix `vv' = J(colsof(`bb'), colsof(`bb'), .)
    }
    restore
    return scalar estimate = scalar(`estimate')
    return scalar native_se = scalar(`se')
    return scalar native_p = scalar(`p')
    return scalar df = `df'
    return scalar n_clusters = `G'
    return scalar n_stacked = `N'
    return matrix b = `bb'
    return matrix V = `vv'
    return matrix event_study = `es'
    local cohorts = strtrim("`cohorts'")
    return local cohorts "`cohorts'"
    if "`covariates'" == "" return local target_label "trimmed aggregate ATT: equal mean over post event times"
    else return local target_label "regression-adjusted Q-weighted post-event contrast"
end
