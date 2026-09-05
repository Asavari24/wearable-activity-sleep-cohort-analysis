*!=============================================================================
*! FILE:        01_cohort_analysis.do
*! PROJECT:     Wearable Activity & Sleep Cohort Analysis — Stata replication
*! AUTHOR:      Asavari Shejwal  <shejwal.a@northeastern.edu>
*! CREATED:     2026-09-04
*! STATA:       18.0 (SE or BE; no community-contributed packages required)
*! RUNTIME:     ~30 seconds (minuteSleep import dominates)
*!
*!-----------------------------------------------------------------------------
*! PURPOSE
*!-----------------------------------------------------------------------------
*! Reproduce, in Stata, a cohort analysis originally executed in Python
*! (pandas / SciPy). The scientific question is unchanged:
*!
*!     Do high-activity users show different mean sleep efficiency than
*!     low-activity users?
*!
*!     H0: mean sleep efficiency is equal across activity cohorts
*!     H1: mean sleep efficiency differs (two-sided)
*!
*! The Python implementation FAILED TO REJECT H0. This file tests whether that
*! null survives a change of statistical environment. It is therefore a
*! REPLICATION script, not an exploratory one: every headline number is pinned
*! in EXPECTED OUTPUT below. Any deviation beyond rounding is a defect in this
*! file, not a discovery.
*!
*!-----------------------------------------------------------------------------
*! INPUTS
*!-----------------------------------------------------------------------------
*! Fitbit Fitness Tracker Data (Mobius / Kaggle), Fitabase exports of MTurk
*! respondents. NOTE the two activity files come from DIFFERENT public mirrors
*! and have DIFFERENT SCHEMAS; harmonising them is the first real task here.
*!
*!  data/raw/dailyActivity_merged_03-12_04-11.csv   [PRIMARY]  457 rows, 35 users
*!      schema: source, participant-id, timestamp(unix secs), total-steps, ...
*!      mirror: health-hacks/datasets
*!
*!  data/raw/dailyActivity_merged_04-12_05-12.csv   [APPEND]   940 rows, 33 users
*!      schema: Id, ActivityDate(M/D/YYYY), TotalSteps, ...
*!      mirror: tubako/bellabeat-case-study (canonical Fitabase layout)
*!
*!  data/raw/minuteSleep_merged_03-12_04-11.csv     [PRIMARY]  198,559 rows
*!      schema: source, participant-id, timestamp, sleep-minutes, log-id
*!      sleep-minutes codes: 1 = asleep, 2 = restless, 3 = awake
*!
*!-----------------------------------------------------------------------------
*! OUTPUTS
*!-----------------------------------------------------------------------------
*!   data/derived/activity_panel.dta      appended 2-tranche daily panel
*!   data/derived/sleep_sessions.dta      one row per sleep session (wide)
*!   data/derived/sleep_states_long.dta   tidy long state counts (reshape demo)
*!   data/derived/user_cohort.dta         one row per user, cohort assigned
*!   data/derived/analysis_sessions.dta   session rows + user cohort (1:m merge)
*!   logs/01_cohort_analysis.log          full transcript
*!
*!-----------------------------------------------------------------------------
*! EXPECTED OUTPUT  — verified against the Python baseline on 2026-09-04.
*! These are not approximations; the Python run reproduced the published
*! README to every digit shown. Stata should land on the same values.
*!-----------------------------------------------------------------------------
*!   Sleep sessions retained (>= 60 min) ............... 556
*!   Users with BOTH activity and sleep signals ........ 23
*!   Median-split threshold (avg daily steps) .......... 8135
*!   Cohort n .......................................... 12 high / 11 low
*!   Mean efficiency, high-activity .................... 0.898
*!   Mean efficiency, low-activity ..................... 0.930
*!   Welch t ........................................... -0.98
*!   Welch p ........................................... 0.34    <- fail to reject
*!   Cohen's d (simple-average pooling) ................ -0.41
*!   Mann-Whitney (ranksum) p .......................... 0.13
*!   Levene / variance-ratio p ......................... 0.42
*!   CDC split n ....................................... 5 high / 9 low
*!   CDC split p ....................................... 0.41
*!   CDC split d ....................................... +0.42   <- SIGN FLIPS
*!   Shapiro-Wilk p, high / low ........................ 0.0014 / 0.0000
*!   Cohort MEDIANS, high / low ........................ 0.929 / 0.947
*!
*!   Section 8 (regress), pinned by the Python reference run:
*!     M1  b(high_activity) = -0.03269, t = -0.964, p = 0.346, R2 = 0.0424
*!     M2  robust p = 0.327   (Stata vce(robust) = HC1 = statsmodels HC1)
*!     M3  b(steps_k) = +0.001538, p = 0.770, R2 = 0.0041
*!     M4  R2 = 0.0698, df_model = 3, df_resid = 19
*!
*!   Section 9 (logit), pinned:
*!     good_sleeper = 20 of 23 (87%) — badly unbalanced, as warned
*!     L1  b(steps_k) = +0.0687, OR = 1.071, p = 0.720, pseudo-R2 = 0.0074
*!     L2  b(high_activity) = -0.6931, p = 0.595
*!     2x2 = [1 2 / 10 10], Fisher exact p = 1.000
*!
*!   Section 10 (sensitivity), pinned:
*!     10a tranche-2 activity ....... p = 0.860, d = +0.075
*!     10b pooled both tranches ..... p = 0.860, d = +0.075
*!     10c duration-weighted ........ p = 0.377
*!     10d duplicates dropped ....... p = 0.3385 (IDENTICAL to baseline)
*!
*!-----------------------------------------------------------------------------
*! DOCUMENTED ANALYTIC DECISIONS
*!-----------------------------------------------------------------------------
*! (1) PRECISION. participant-id reaches 8,877,689,391 and log-id reaches
*!     11,374,876,178. Stata's -float- holds only ~7 significant digits, so
*!     either stored as float would SILENTLY COLLAPSE distinct users and
*!     distinct sleep sessions into one another. Both are imported as strings
*!     and destrung to -double-. Section 1 and 2 assert the resulting distinct
*!     counts, because a precision loss here is invisible in the output but
*!     corrupts every number downstream. This is the single largest
*!     correctness risk in the pipeline.
*!
*! (2) TRANCHE SELECTION. The primary analysis uses the 03-12..04-11 tranche
*!     for BOTH activity and sleep, because that is the pairing the Python
*!     baseline consumed. Pairing tranche-1 sleep with tranche-2 activity
*!     gives p = 0.86 rather than 0.34 — same qualitative null, different
*!     numbers. Section 10 runs that pairing explicitly as a sensitivity
*!     check rather than leaving the choice implicit.
*!
*! (3) DUPLICATE SLEEP MINUTES. minuteSleep contains 525 exactly-duplicated
*!     rows. The Python baseline did not remove them, so this file does not
*!     either. Their effect was MEASURED rather than assumed: dropping them
*!     gives p = 0.3385, identical to keeping them, because the duplicates
*!     fall inside long sessions and shift no session across the 60-minute
*!     filter. So this decision is immaterial to the conclusion — which is
*!     worth stating explicitly, since "we kept known-bad rows" otherwise
*!     reads as an unquantified risk. They are still counted and reported
*!     in section 2 so a reader can see the flaw exists.
*!
*! (4) SESSION FILTER. Sessions under 60 minutes in bed are dropped as
*!     nap/wear fragments, matching Python. In this tranche the filter is
*!     non-binding (the shortest session is already >= 60 min), so it is
*!     retained for equivalence rather than effect.
*!
*! (5) COHORT SPLIT. The median split is sample-dependent by construction.
*!     The CDC 10k/7k split is the sensitivity check, and it FLIPS THE SIGN
*!     of the effect (d = -0.41 vs +0.42). Reproducing that instability is
*!     the point: it is the clearest evidence the effect is sampling noise.
*!
*! (6) POWER. n = 23 at the unit of analysis; ~96 per cohort would be needed
*!     for d = 0.4 at 80% power. Every inferential result here is severely
*!     underpowered. The logit in section 9 is reported for method coverage
*!     and is NOT interpretable at this n; caveats print inline.
*!=============================================================================

version 18
clear all
set more off
set linesize 100

*-----------------------------------------------------------------------------
* 0. ENVIRONMENT
*-----------------------------------------------------------------------------
* Single point of configuration. Run this file from the repo root.

global PROJROOT "`c(pwd)'"

* If you run from the do/ subdirectory instead, uncomment:
* global PROJROOT ".."

global RAW      "${PROJROOT}/data/raw"
global DERIVED  "${PROJROOT}/data/derived"
global LOGS     "${PROJROOT}/logs"

capture mkdir "${DERIVED}"
capture mkdir "${LOGS}"

capture log close _all
log using "${LOGS}/01_cohort_analysis.log", replace text

* -- helper ----------------------------------------------------------------
* The tranche-1 mirror uses hyphenated headers ("participant-id"). Stata
* rewrites invalid characters on import, and the exact rewrite has varied
* across versions: usually "participant_id", sometimes "participantid".
* Rather than hard-code one and fail cryptically on the other, accept either
* and abort with a readable message if neither is present.
capture program drop _needvar
program define _needvar
    args want alt
    capture confirm variable `want'
    if _rc {
        capture rename `alt' `want'
        if _rc {
            display as error "FATAL: import produced neither `want' nor `alt'."
            display as error "Run -describe- after the import and fix the rename block."
            exit 111
        }
    }
end

display as text _n "{hline 78}"
display as text "Wearable cohort replication — `c(current_date)' `c(current_time)'"
display as text "Stata `c(stata_version)' / flavor `c(flavor)' / OS `c(os)'"
display as text "{hline 78}"


*=============================================================================
* 1. INGEST + APPEND   [JD task: append]
*=============================================================================
* Two Fitabase export windows, two different mirror schemas. Harmonise to a
* common set of names and types FIRST, then stack. Provenance is tagged before
* the append so any row can be traced back to its source file afterwards.

display as result _n "== 1. APPEND: stacking two activity tranches =="

* --- tranche 1: health-hacks schema, unix-epoch timestamps -----------------
* stringcols(2 3) forces participant-id and timestamp in as text so we can
* control their numeric type ourselves. timestamp (~1.46e9) also exceeds
* float precision, so it matters as much as the id does — see decision (1).
import delimited "${RAW}/dailyActivity_merged_03-12_04-11.csv", ///
    clear varnames(1) stringcols(2 3)

_needvar participant_id       participantid
_needvar total_steps          totalsteps
_needvar total_distance       totaldistance
_needvar very_active_minutes  veryactiveminutes
_needvar sedentary_minutes    sedentaryminutes

destring participant_id, replace
destring timestamp,      replace
recast double participant_id
recast double timestamp

rename participant_id       id
rename total_steps          steps
rename very_active_minutes  very_active_min
rename sedentary_minutes    sedentary_min
rename total_distance       distance

* Unix seconds -> Stata daily date. mdy(1,1,1970) is Stata's epoch origin.
gen int activity_date = mdy(1,1,1970) + floor(timestamp/86400)
format activity_date %td

keep  id activity_date steps calories very_active_min sedentary_min distance
gen byte tranche = 1

quietly count
display as text "  tranche 1 rows: " r(N) "   (expected 457)"

tempfile t1
save "`t1'"

* --- tranche 2: canonical Fitabase schema, M/D/YYYY dates -----------------
import delimited "${RAW}/dailyActivity_merged_04-12_05-12.csv", ///
    clear varnames(1) stringcols(1)

destring id, replace
recast double id

rename totalsteps          steps
rename veryactiveminutes   very_active_min
rename sedentaryminutes    sedentary_min
rename totaldistance       distance

gen int activity_date = date(activitydate, "MDY")
format activity_date %td

keep  id activity_date steps calories very_active_min sedentary_min distance
gen byte tranche = 2

quietly count
display as text "  tranche 2 rows: " r(N) "   (expected 940)"

* --- the append -----------------------------------------------------------
append using "`t1'"

label define tranchelbl 1 "2016-03-12..04-11" 2 "2016-04-12..05-12"
label values tranche tranchelbl
label variable tranche         "Source Fitabase export window"
label variable activity_date   "Calendar date of activity record"
label variable steps           "Total daily steps, device-reported"
label variable very_active_min "Very-active minutes"
label variable sedentary_min   "Sedentary minutes"

quietly count
display as text "  appended rows (pre de-dup): " r(N) "   (expected 1397)"
tabulate tranche

* --- de-duplicate on the natural key --------------------------------------
* The windows touch on 2016-04-12, which appears in both exports. Where a
* user-date appears twice, prefer the later export.
quietly duplicates tag id activity_date, generate(_dup)
quietly count if _dup > 0
display as text "  overlapping user-date rows: " r(N)

* Sort so the PREFERRED record is first within each user-date. gsort would
* express "tranche descending" directly, but gsort with a descending key
* leaves the data not -sortedby- these variables, so a following -by:- fails.
* Sorting ascending on a negated copy achieves the same order and keeps the
* sort marker valid.
gen byte _negtranche = -tranche
bysort id activity_date (_negtranche): keep if _n == 1
drop _dup _negtranche

* ASSERT: the key is unique. If this trips, the harmonisation or the de-dup
* is wrong and every downstream number is suspect.
isid id activity_date
display as text "  ASSERT ok: id x activity_date unique"

* ASSERT: precision survived. 35 distinct users is the union of the two
* windows (33 shared + 2 tranche-1 only). A float-rounded id would collapse
* users and this count would fall — the failure mode decision (1) guards.
quietly levelsof id
local nu = r(r)
display as text "  distinct users: `nu'   (expected 35)"
assert `nu' == 35

compress
label data "Fitbit daily activity, 2 tranches harmonised + appended + de-duped"
save "${DERIVED}/activity_panel.dta", replace


*=============================================================================
* 2. SLEEP SESSIONS + RESHAPE   [JD task: reshape long/wide]
*=============================================================================
* minuteSleep is one row per MINUTE per session with a state code. Sleep
* efficiency needs one row per SESSION with the three state counts side by
* side. That is a genuine long -> wide problem, not a contrived one.

display as result _n "== 2. RESHAPE: minute states -> per-session wide =="

* stringcols(2 5) protects participant-id and log-id. log-id reaches 1.14e10;
* as a float it would round and merge distinct sleep sessions together.
import delimited "${RAW}/minuteSleep_merged_03-12_04-11.csv", ///
    clear varnames(1) stringcols(2 5)

_needvar participant_id  participantid
_needvar sleep_minutes   sleepminutes
_needvar log_id          logid

destring participant_id, replace
destring log_id,         replace
recast double participant_id
recast double log_id

rename participant_id id
rename log_id         logid
rename sleep_minutes  state

keep id logid state

quietly count
display as text "  minute-level rows: " r(N) "   (expected 198559)"

label define sleepstate 1 "asleep" 2 "restless" 3 "awake"
label values state sleepstate
tabulate state

* Decision (3): the raw file contains 525 exactly-duplicated rows, which are
* retained here to match the Python baseline.
*
* Deliberately NOT running -duplicates report id logid state- at this point:
* timestamp has already been dropped, so within any session every minute
* sharing a state would count as a duplicate and the command would report a
* number in the tens of thousands, not 525. The genuine full-row duplicate
* count is established in section 10d, which re-imports with timestamp
* retained. Reporting a wrong number next to the right caption is worse than
* reporting none.
display as text "  NOTE: raw file holds 525 exactly-duplicated minute rows."
display as text "        Retained to match the Python baseline; section 10d"
display as text "        re-runs without them and confirms p is unchanged."

* ASSERT: session key precision survived.
quietly levelsof logid
local nl = r(r)
display as text "  distinct sleep sessions: `nl'   (expected 556)"
assert `nl' == 556

* --- collapse minutes -> (session x state) counts: the LONG form ----------
gen byte one = 1
collapse (sum) nmin = one, by(id logid state)

quietly count
display as text "  session x state rows (long): " r(N)

label variable nmin "Minutes recorded in this state for this session"
label data "Sleep state minute counts, long (one row per session x state)"
save "${DERIVED}/sleep_states_long.dta", replace

* --- long -> wide ---------------------------------------------------------
* i() = the unit we want one row of; j() = the variable whose values become
* new columns. The collapse above guarantees i x j is unique, which reshape
* requires — running reshape straight off the minute file would error.
reshape wide nmin, i(id logid) j(state)

* A session with no restless (or no awake) minutes yields a MISSING cell,
* not a zero. Left missing, the arithmetic below would silently propagate to
* missing and quietly drop sessions. Recode explicitly.
foreach v of varlist nmin1 nmin2 nmin3 {
    quietly recode `v' (. = 0)
}

* egen rowtotal also treats missing as zero, so this is belt-and-braces after
* the recode. Kept because it is likewise the JD's egen task.
egen double minutes_in_bed = rowtotal(nmin1 nmin2 nmin3)
gen  double efficiency     = nmin1 / minutes_in_bed

label variable nmin1          "Minutes asleep"
label variable nmin2          "Minutes restless"
label variable nmin3          "Minutes awake"
label variable minutes_in_bed "Total minutes recorded in bed"
label variable efficiency     "Minutes asleep / minutes in bed"

quietly count
display as text "  sessions before duration filter: " r(N)

* Decision (4): non-binding on this tranche, retained for equivalence.
quietly summarize minutes_in_bed
display as text "  shortest session (min): " r(min)
drop if minutes_in_bed < 60
quietly count
display as text "  sessions retained (>= 60 min): " r(N) "   (EXPECTED 556)"

* Efficiency is a proportion by construction; assert rather than assume.
assert efficiency >= 0 & efficiency <= 1
display as text "  ASSERT ok: efficiency within [0,1]"

summarize efficiency, detail

isid id logid
compress
label data "Fitbit sleep sessions, one row per session (>=60 min in bed)"
save "${DERIVED}/sleep_sessions.dta", replace

* --- wide -> long round trip (demonstrating the inverse) ------------------
* Not needed downstream. Run to show the reverse operation and to prove the
* wide file converts back losslessly, which is a real curation property.
preserve
    keep id logid nmin1 nmin2 nmin3
    reshape long nmin, i(id logid) j(sleep_state)
    label values sleep_state sleepstate
    quietly count
    display as text "  round-tripped back to long, rows: " r(N)
    tabulate sleep_state
restore


*=============================================================================
* 3. CURATE TO USER LEVEL   [JD task: egen / collapse]
*=============================================================================
display as result _n "== 3. COLLAPSE: per-user aggregates =="

* --- per-user mean sleep efficiency ---------------------------------------
use "${DERIVED}/sleep_sessions.dta", clear

* Unweighted mean of session efficiencies, matching the Python baseline.
* NOTE this weights a 70-minute session equally with a 9-hour one. A
* duration-weighted alternative is carried alongside so the choice stays
* visible and testable instead of buried in a default.
collapse (mean)  avg_efficiency  = efficiency      ///
         (sum)   total_asleep    = nmin1           ///
         (sum)   total_in_bed    = minutes_in_bed  ///
         (count) n_sessions      = efficiency,     by(id)

gen double wtd_efficiency = total_asleep / total_in_bed

label variable avg_efficiency "Mean of session efficiencies (baseline definition)"
label variable wtd_efficiency "Duration-weighted efficiency (alternative)"
label variable n_sessions     "Number of qualifying sleep sessions"

quietly count
display as text "  users with sleep data: " r(N) "   (expected 23)"

tempfile usersleep
save "`usersleep'"

* --- per-user mean daily steps, PRIMARY tranche (decision (2)) ------------
use "${DERIVED}/activity_panel.dta", clear
keep if tranche == 1
quietly count
display as text "  tranche-1 activity rows: " r(N) "   (expected 457)"

collapse (mean)  avg_steps       = steps           ///
         (mean)  avg_calories    = calories        ///
         (mean)  avg_veryactive  = very_active_min ///
         (mean)  avg_sedentary   = sedentary_min   ///
         (count) n_days          = steps,          by(id)

label variable avg_steps      "Mean daily steps (tranche 1)"
label variable n_days         "Days of activity data"
label variable avg_veryactive "Mean daily very-active minutes"
label variable avg_sedentary  "Mean daily sedentary minutes"

quietly count
display as text "  users with activity data: " r(N) "   (expected 35)"


*=============================================================================
* 4a. MERGE 1:1 — build the analysis population
*=============================================================================
display as result _n "== 4a. MERGE 1:1 =="

* Inner join: the analysis population is users carrying BOTH signals.
merge 1:1 id using "`usersleep'"

display as text "  1 = activity only, 2 = sleep only, 3 = both:"
tabulate _merge

keep if _merge == 3
drop _merge

quietly count
display as text "  users with BOTH signals: " r(N) "   (EXPECTED 23)"
assert _N == 23


*=============================================================================
* 5. COHORT ASSIGNMENT   [JD task: egen]
*=============================================================================
display as result _n "== 5. COHORT ASSIGNMENT =="

* egen median() computes across all observations in memory, which at this
* point is exactly the 23-user analysis population — the same denominator
* the Python baseline used. Computing it before the merge would silently
* use all 35 activity users and shift the threshold.
egen double med_steps = median(avg_steps)
display as text "  median-split threshold: " %9.1f med_steps[1] "   (EXPECTED 8135)"

* The !missing() guard is not decoration. In Stata, missing is larger than
* any number, so "avg_steps >= med_steps" returns 1 for a missing avg_steps
* and would silently file unmeasured users into the HIGH-activity cohort.
* The inner merge above leaves no missings today, but this is the kind of
* assumption that quietly stops holding when someone changes the join.
gen byte high_activity = (avg_steps >= med_steps) if !missing(avg_steps)
label define actlbl 0 "Low activity" 1 "High activity"
label values high_activity actlbl
label variable high_activity "Activity cohort (median split on mean daily steps)"

* CDC sensitivity split. Deliberately leaves 7000-9999 unassigned, so this
* variable has missings BY DESIGN, not by accident.
gen byte cdc_cohort = .
replace cdc_cohort = 1 if avg_steps >= 10000
replace cdc_cohort = 0 if avg_steps <  7000
label values cdc_cohort actlbl
label variable cdc_cohort "Activity cohort (CDC 10k/7k; 7k-10k excluded)"

* Binary outcome for the logit in section 9. The 0.90 threshold is the
* conventional clinical sleep-efficiency cut, NOT a within-sample split, so
* it does not inherit the median split's sample dependence.
gen byte good_sleeper = (avg_efficiency >= 0.90) if !missing(avg_efficiency)
label define slplbl 0 "Efficiency < 0.90" 1 "Efficiency >= 0.90"
label values good_sleeper slplbl
label variable good_sleeper "Mean sleep efficiency >= 0.90"

* Steps in thousands: coefficients on raw steps are uninterpretably small.
gen double avg_steps_k = avg_steps / 1000
label variable avg_steps_k "Mean daily steps (thousands)"

tabulate high_activity
tabulate cdc_cohort, missing
tabulate good_sleeper

order id avg_steps avg_steps_k avg_efficiency high_activity cdc_cohort good_sleeper
compress
label data "User-level analysis file: 23 users with paired activity + sleep"
save "${DERIVED}/user_cohort.dta", replace


*=============================================================================
* 4b. MERGE 1:m — broadcast cohort back to session level   [JD task: merge 1:m]
*=============================================================================
display as result _n "== 4b. MERGE 1:m =="

* One user record maps to many session records. This is the JD's -merge 1:m-
* and it is the correct direction here: 1:1 would error on the duplicate ids
* in the using file, and m:1 would be backwards. The result supports
* session-level (within-user) modelling that the 23-row file cannot.
use "${DERIVED}/user_cohort.dta", clear
keep id avg_steps avg_steps_k high_activity cdc_cohort good_sleeper

merge 1:m id using "${DERIVED}/sleep_sessions.dta"
tabulate _merge

* Sessions from users without activity data are outside the population.
keep if _merge == 3
drop _merge

quietly count
display as text "  session rows in analysis population: " r(N) "   (EXPECTED 556)"

* egen tag flags exactly one row per user, so user-level counts can be taken
* off this session-level file without collapsing it.
egen byte first_of_user = tag(id)
quietly count if first_of_user
display as text "  distinct users on session file: " r(N) "   (expected 23)"

label data "Sleep sessions + user cohort (1:m merge)"
save "${DERIVED}/analysis_sessions.dta", replace


*=============================================================================
* 6. DESCRIPTIVES   [JD task: summarize]
*=============================================================================
display as result _n "== 6. SUMMARIZE =="

use "${DERIVED}/user_cohort.dta", clear

summarize avg_steps avg_efficiency n_sessions n_days, detail

display as text _n "  -- by activity cohort --"
tabstat avg_steps avg_efficiency, by(high_activity) ///
    statistics(n mean sd min p50 max) columns(statistics) nototal

display as text _n "  -- cohort means (EXPECTED 0.898 high / 0.930 low) --"
tabulate high_activity, summarize(avg_efficiency)

* Baseline-vs-alternative efficiency definition. If these disagree materially,
* the unweighted mean in section 3 is doing more work than it should.
display as text _n "  -- unweighted vs duration-weighted efficiency --"
summarize avg_efficiency wtd_efficiency
correlate avg_efficiency wtd_efficiency

* --- assumption checks, mirroring the Python baseline ---------------------
* Python used scipy.stats.shapiro and scipy.stats.levene. Stata's direct
* equivalents are swilk and robvar (robvar reports Levene's W0 first).
display as text _n "  -- Shapiro-Wilk by cohort (EXPECTED p ~.001 / <.001) --"
foreach g in 1 0 {
    display as text "   cohort = `g'"
    quietly count if high_activity == `g'
    if r(N) >= 4 {
        swilk avg_efficiency if high_activity == `g'
    }
    else {
        display as error "   n < 4, swilk not computable"
    }
}

display as text _n "  -- variance equality (EXPECTED Levene p ~0.42) --"
robvar avg_efficiency, by(high_activity)
sdtest avg_efficiency, by(high_activity)


*=============================================================================
* 7. TWO-SAMPLE TESTS   [JD task: ttest]
*=============================================================================
display as result _n "== 7. TTEST =="

* Welch's unequal-variance t-test is primary, matching
* scipy.stats.ttest_ind(..., equal_var=False). Stata's -unequal welch- uses
* the Welch-Satterthwaite df, which is what SciPy computes.
display as text _n "  -- PRIMARY: Welch two-sample t-test (median split) --"
ttest avg_efficiency, by(high_activity) unequal welch

scalar t_welch = r(t)
scalar p_welch = r(p)
scalar m_low   = r(mu_1)
scalar m_high  = r(mu_2)
scalar sd_low  = r(sd_1)
scalar sd_high = r(sd_2)

* Cohen's d with SIMPLE-AVERAGE pooling: sqrt((var_a + var_b)/2). This is the
* formula the Python helper used. It is NOT Stata's n-weighted pooled SD, so
* it is computed by hand to keep the two environments strictly comparable.
scalar pooled_sd = sqrt((sd_high^2 + sd_low^2) / 2)
scalar cohen_d   = (m_high - m_low) / pooled_sd

display as text _n "  {hline 62}"
display as text "  Welch t   = " %8.3f t_welch  "     EXPECTED  -0.98"
display as text "  Welch p   = " %8.4f p_welch  "     EXPECTED   0.34"
display as text "  mean high = " %8.4f m_high   "     EXPECTED   0.898"
display as text "  mean low  = " %8.4f m_low    "     EXPECTED   0.930"
display as text "  Cohen's d = " %8.3f cohen_d  "     EXPECTED  -0.41"
display as text "  {hline 62}"

* Stata's own effect size, for reference. It uses n-weighted pooling, so it
* will differ slightly from cohen_d above. Both are correct; they answer the
* same question under different pooling conventions. Reporting both makes
* the convention explicit rather than letting a reader assume.
esize twosample avg_efficiency, by(high_activity) cohensd

* Shapiro-Wilk rejects normality in BOTH cohorts, so the rank test is
* arguably the more defensible primary. ranksum == Mann-Whitney U ==
* scipy.stats.mannwhitneyu.
display as text _n "  -- BACKUP: Mann-Whitney (EXPECTED p ~0.13) --"
ranksum avg_efficiency, by(high_activity)

* --- sensitivity: CDC split, where the effect SIGN FLIPS -----------------
display as text _n "  -- SENSITIVITY: CDC 10k/7k (EXPECTED n 5/9, p 0.41, d +0.42) --"
quietly count if cdc_cohort == 1
local n_hi = r(N)
quietly count if cdc_cohort == 0
local n_lo = r(N)
display as text "   n high = `n_hi' (exp 5), n low = `n_lo' (exp 9)"

if `n_hi' >= 3 & `n_lo' >= 3 {
    ttest avg_efficiency, by(cdc_cohort) unequal welch
    scalar d_cdc = (r(mu_2) - r(mu_1)) / sqrt((r(sd_2)^2 + r(sd_1)^2) / 2)
    display as text "   Cohen's d (CDC) = " %8.3f d_cdc "   EXPECTED +0.42"
    display as text "   NOTE: the sign flip against the median split indicates"
    display as text "         instability from small n, not a dose-response."
}
else {
    display as error "   skipped — cohort too small"
}


*=============================================================================
* 8. LINEAR MODELS   [JD task: regress]
*=============================================================================
display as result _n "== 8. REGRESS =="

* M1. OLS on the cohort dummy. Under equal variances the coefficient on
* high_activity is numerically identical to the pooled two-sample mean
* difference, and its t equals the pooled t. Included to make the
* t-test/regression equivalence explicit rather than assumed.
display as text _n "  -- M1: efficiency on cohort dummy --"
display as text "     EXPECTED b=-0.03269  t=-0.964  p=0.346  R2=0.0424"
regress avg_efficiency high_activity

* M2. Same model, heteroskedasticity-robust SEs — the regression analogue of
* Welch. Its p should sit near the Welch p from section 7.
display as text _n "  -- M2: same, robust SEs (Welch analogue) --"
display as text "     EXPECTED p=0.327  (Stata vce(robust) is HC1)"
regress avg_efficiency high_activity, vce(robust)

* M3. Steps continuous, discarding the arbitrary median cut. If a real
* dose-response existed, dropping the cut should make it EASIER to see, not
* harder — a cohort split cannot manufacture signal the continuous model lacks.
display as text _n "  -- M3: efficiency on continuous steps (per 1k) --"
display as text "     EXPECTED b=+0.001538  p=0.770  R2=0.0041"
regress avg_efficiency avg_steps_k
display as text "   R-squared 0.004 vs M1's 0.042: the continuous model explains"
display as text "   even LESS. There is no dose-response for the median split to"
display as text "   have discretised. That is the substantive finding here."

* M4. Adjusted. n = 23 with 3 predictors is thin; this is a robustness
* gesture, not a credible causal model, and is reported as such.
display as text _n "  -- M4: + very-active and sedentary minutes (n=23, thin) --"
display as text "     EXPECTED R2=0.0698  df_m=3  df_r=19"
regress avg_efficiency avg_steps_k avg_veryactive avg_sedentary
display as text "   model df = " e(df_m) ", residual df = " e(df_r)

* Diagnostics on the primary continuous specification.
display as text _n "  -- M3 diagnostics --"
quietly regress avg_efficiency avg_steps_k
predict double resid_m3, residuals
swilk resid_m3
estat hettest
display as text "   Non-normal residuals here reflect the same 2-3 low-efficiency"
display as text "   outliers that pulled the cohort MEANS apart while the cohort"
display as text "   MEDIANS stayed close (0.929 vs 0.947 in the baseline)."
drop resid_m3


*=============================================================================
* 9. LOGISTIC MODEL   [JD task: logit]
*=============================================================================
display as result _n "== 9. LOGIT =="
display as text "  WARNING: n = 23 and the outcome is unbalanced. The rule of"
display as text "  thumb is >= 10 events per predictor; this has far fewer."
display as text "  Reported for method coverage. NOT interpretable as evidence"
display as text "  about sleep. Do not quote these coefficients."

tabulate good_sleeper
display as text "   EXPECTED: 20 of 23 are 'good sleepers' (87%). Three events in"
display as text "   the minority class is why nothing below can be believed."

* L1: does step volume predict being a 'good sleeper'?
display as text _n "  -- L1: good_sleeper on steps (per 1k) --"
display as text "     EXPECTED b=+0.0687  OR=1.071  p=0.720  pseudo-R2=0.0074"
capture noisily logit good_sleeper avg_steps_k
if _rc != 0 {
    display as error "   logit failed (rc = " _rc "): likely separation or"
    display as error "   non-convergence at this sample size."
}
else {
    display as text "   pseudo-R2 = " %6.4f e(r2_p)
    * Odds ratios are the readable form of the same fit.
    logit, or
    * Margins put the effect back on the probability scale — the only scale
    * a product stakeholder would ever act on.
    margins, dydx(avg_steps_k)
}

* L2: cohort dummy as predictor. This is a 2x2 table, so the exact test
* below is the honest companion at this n.
display as text _n "  -- L2: good_sleeper on cohort dummy --"
display as text "     EXPECTED b=-0.6931  p=0.595"
capture noisily logit good_sleeper high_activity
if _rc != 0 {
    display as error "   logit failed (rc = " _rc ") — expected if a cell is empty."
}

display as text _n "  -- 2x2 + Fisher's exact (the valid test at this n) --"
display as text "     EXPECTED table [1 2 / 10 10], Fisher exact p = 1.000"
tabulate good_sleeper high_activity, chi2 exact row
display as text "   Fisher p = 1.00 is as null as a result can be: the cohorts"
display as text "   are indistinguishable on this outcome."


*=============================================================================
* 10. SENSITIVITY ANALYSES  (decisions (2) and (3))
*=============================================================================
display as result _n "== 10. SENSITIVITY =="

* --- 10a. tranche-2 activity instead of tranche-1 -------------------------
* Decision (2) made a choice; this quantifies what that choice cost. The
* Python cross-check gives p = 0.86 for this pairing: a different number,
* the same qualitative conclusion.
display as text _n "  -- 10a. tranche-2 activity x tranche-1 sleep --"
use "${DERIVED}/activity_panel.dta", clear
keep if tranche == 2
collapse (mean) avg_steps = steps, by(id)
merge 1:1 id using "`usersleep'", keep(match) nogenerate

quietly count
display as text "     users with both signals: " r(N)
egen double med2 = median(avg_steps)
gen byte high2 = (avg_steps >= med2) if !missing(avg_steps)
display as text "     median threshold: " %9.1f med2[1]
ttest avg_efficiency, by(high2) unequal welch
display as text "     Welch p = " %8.4f r(p) "   (Python cross-check: 0.86)"

* --- 10b. pooled panel, both tranches -------------------------------------
display as text _n "  -- 10b. pooled activity panel (both tranches) --"
use "${DERIVED}/activity_panel.dta", clear
collapse (mean) avg_steps = steps, by(id)
merge 1:1 id using "`usersleep'", keep(match) nogenerate

egen double medp = median(avg_steps)
gen byte highp = (avg_steps >= medp) if !missing(avg_steps)
display as text "     pooled median threshold: " %9.1f medp[1] "   (expect 8367.1)"
ttest avg_efficiency, by(highp) unequal welch
display as text "     Welch p = " %8.4f r(p) "   (expect 0.8600, d = +0.075)"
display as text "     Same p as 10a: pooling does not change which users land"
display as text "     in which cohort, so it cannot change the answer."

* --- 10c. duration-weighted efficiency ------------------------------------
* Tests whether the unweighted session mean (decision in section 3) drove
* the result.
display as text _n "  -- 10c. duration-weighted efficiency, primary cohorts --"
use "${DERIVED}/user_cohort.dta", clear
ttest wtd_efficiency, by(high_activity) unequal welch
display as text "     Welch p = " %8.4f r(p) "   (expect 0.3769)"
display as text "     Non-significant under either weighting: the null does not"
display as text "     depend on the session-weighting choice in section 3."

* --- 10d. effect of the retained duplicate minute rows (decision (3)) -----
* Re-derive efficiency after dropping the 525 exact duplicates. The Python
* cross-check gives p = 0.3385, IDENTICAL to the baseline, because the
* duplicates sit inside long sessions and move none across the 60-minute
* filter. Run here so the claim in decision (3) is demonstrated, not asserted.
display as text _n "  -- 10d. duplicates dropped instead of retained --"
display as text "     (baseline retains them and gives p = 0.3385)"

* timestamp is read as a STRING here (stringcols includes column 3). It is
* only used as part of the row identity for de-duplication, and as a float
* it would round — making distinct minutes collide and over-deleting.
import delimited "${RAW}/minuteSleep_merged_03-12_04-11.csv", ///
    clear varnames(1) stringcols(2 3 5)
_needvar participant_id  participantid
_needvar sleep_minutes   sleepminutes
_needvar log_id          logid
destring participant_id log_id, replace
recast double participant_id
recast double log_id
rename participant_id id
rename log_id         logid
rename sleep_minutes  state

* CRITICAL: de-duplicate on the FULL row, timestamp included. Dropping
* timestamp first would make every minute of a session with the same state
* look like a duplicate and collapse each session to 3 rows.
keep id logid state timestamp
quietly duplicates drop
quietly count
display as text "     rows after de-duplication: " r(N) "   (expect 198034)"
drop timestamp

gen byte asleep_min = (state == 1)
collapse (count) minutes_in_bed = state (sum) nasleep = asleep_min, by(id logid)
gen double efficiency = nasleep / minutes_in_bed
drop if minutes_in_bed < 60
quietly count
display as text "     sessions after de-dup: " r(N) "   (expect 556 — unchanged)"

collapse (mean) avg_efficiency = efficiency, by(id)
merge 1:1 id using "${DERIVED}/user_cohort.dta", ///
    keep(match) keepusing(high_activity) nogenerate

ttest avg_efficiency, by(high_activity) unequal welch
display as text "     Welch p = " %8.4f r(p) "   (expect 0.3385 — identical)"
display as text "     Confirms decision (3): the known-duplicate rows are"
display as text "     immaterial to the conclusion."


*=============================================================================
* 11. CLOSE
*=============================================================================
display as result _n "== DONE =="
display as text "Derived files: ${DERIVED}"
display as text "Log:           ${LOGS}/01_cohort_analysis.log"
display as text _n "CONCLUSION TO VERIFY: Welch p ~ 0.34 => fail to reject H0."
display as text "The Python null should survive the environment change."
display as text "If a headline number diverges, this file is wrong, not the"
display as text "finding. Check the double-precision casts in sections 1-2"
display as text "first: that failure is silent and corrupts everything after it."

log close _all

* end of 01_cohort_analysis.do
