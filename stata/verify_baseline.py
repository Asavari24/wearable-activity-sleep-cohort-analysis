"""
Python reference run for the Stata replication.

Two jobs:
  1. Confirm the mirrored CSVs reproduce the published README of
     Asavari24/wearable-activity-sleep-cohort-analysis exactly.
  2. Pin the numbers for every section of do/01_cohort_analysis.do, including
     the regression and logit sections that the original write-up did not
     report, so the Stata run has a target for all of them.

Run from the python/ directory:  ../.venv/bin/python verify_baseline.py

IMPORTANT: the baseline pairs tranche-1 ACTIVITY with tranche-1 SLEEP
(2016-03-12..04-11). Pairing tranche-1 sleep with tranche-2 activity gives
p = 0.86 — same null, different numbers. See decision (2) in the .do file.
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

RAW = "../data/raw"

PUBLISHED = {
    "n_users": 23, "n_sessions": 556, "median_steps": 8135,
    "n_high": 12, "n_low": 11,
    "mean_high": 0.898, "mean_low": 0.930,
    "t": -0.98, "p": 0.34, "d": -0.41,
    "p_mw": 0.13, "levene_p": 0.42,
    "cdc_n_high": 5, "cdc_n_low": 9, "cdc_p": 0.41, "cdc_d": 0.42,
}


def cohens_d(a, b):
    """Simple-average pooling, matching the original helper (not n-weighted)."""
    return (a.mean() - b.mean()) / np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)


# ---------------------------------------------------------------- load ----
sleep = (
    pd.read_csv(f"{RAW}/minuteSleep_merged_03-12_04-11.csv")
    .rename(columns={"participant-id": "Id", "sleep-minutes": "value",
                     "log-id": "logId"})
)
act1 = (
    pd.read_csv(f"{RAW}/dailyActivity_merged_03-12_04-11.csv")
    .rename(columns={"participant-id": "Id", "total-steps": "TotalSteps",
                     "very-active-minutes": "VeryActiveMinutes",
                     "sedentary-minutes": "SedentaryMinutes"})
)
act2 = pd.read_csv(f"{RAW}/dailyActivity_merged_04-12_05-12.csv")

print("=== INPUTS ===")
print(f"sleep t1 : {len(sleep):>7} rows | {sleep.Id.nunique()} users "
      f"| {sleep.duplicated().sum()} exact dup rows (RETAINED, see decision 3)")
print(f"act   t1 : {len(act1):>7} rows | {act1.Id.nunique()} users")
print(f"act   t2 : {len(act2):>7} rows | {act2.Id.nunique()} users")
print(f"union of activity users: {len(set(act1.Id) | set(act2.Id))}  (Stata asserts 35)")

# ------------------------------------------------------------ sessions ----
sess = (
    sleep.groupby(["Id", "logId"])
    .agg(in_bed=("value", "count"), asleep=("value", lambda v: (v == 1).sum()))
    .reset_index()
)
sess["efficiency"] = sess["asleep"] / sess["in_bed"]
sess = sess[sess["in_bed"] >= 60]

user = pd.DataFrame({
    "avg_efficiency": sess.groupby("Id")["efficiency"].mean(),
    "wtd_efficiency": sess.groupby("Id")["asleep"].sum() / sess.groupby("Id")["in_bed"].sum(),
    "n_sessions": sess.groupby("Id")["efficiency"].count(),
})
act_u = act1.groupby("Id").agg(
    avg_steps=("TotalSteps", "mean"),
    avg_veryactive=("VeryActiveMinutes", "mean"),
    avg_sedentary=("SedentaryMinutes", "mean"),
)
u = act_u.join(user, how="inner").dropna(subset=["avg_efficiency"])

med = u["avg_steps"].median()
u["high_activity"] = (u["avg_steps"] >= med).astype(int)
u["avg_steps_k"] = u["avg_steps"] / 1000
u["good_sleeper"] = (u["avg_efficiency"] >= 0.90).astype(int)

hi = u.loc[u.high_activity == 1, "avg_efficiency"]
lo = u.loc[u.high_activity == 0, "avg_efficiency"]
t, p = stats.ttest_ind(hi, lo, equal_var=False)
uu, p_mw = stats.mannwhitneyu(hi, lo, alternative="two-sided")
lev = stats.levene(hi, lo)

hc = u.loc[u.avg_steps >= 10000, "avg_efficiency"]
lc = u.loc[u.avg_steps < 7000, "avg_efficiency"]
tc, pc = stats.ttest_ind(hc, lc, equal_var=False)

got = {
    "n_users": len(u), "n_sessions": len(sess), "median_steps": round(med),
    "n_high": len(hi), "n_low": len(lo),
    "mean_high": round(hi.mean(), 3), "mean_low": round(lo.mean(), 3),
    "t": round(t, 2), "p": round(p, 2), "d": round(cohens_d(hi, lo), 2),
    "p_mw": round(p_mw, 2), "levene_p": round(lev.pvalue, 2),
    "cdc_n_high": len(hc), "cdc_n_low": len(lc),
    "cdc_p": round(pc, 2), "cdc_d": round(cohens_d(hc, lc), 2),
}

print("\n=== REPLICATION vs PUBLISHED README ===")
print(f"{'metric':<14}{'published':>11}{'got':>9}   status")
ok = True
for k, exp in PUBLISHED.items():
    g = got[k]
    match = abs(g - exp) <= (0.015 if isinstance(exp, float) else 0)
    ok &= match
    print(f"{k:<14}{exp:>11}{g:>9}   {'ok' if match else 'MISMATCH'}")
print("\n>>> " + ("MIRRORS FAITHFUL — safe to replicate in Stata" if ok
                  else "MIRRORS DIVERGE — do NOT use for replication"))

sw_hi, sw_lo = stats.shapiro(hi), stats.shapiro(lo)
print(f"\nShapiro-Wilk  high p={sw_hi.pvalue:.4f}  low p={sw_lo.pvalue:.4f}")
print(f"cohort medians: high={hi.median():.3f} low={lo.median():.3f} "
      f"(close, while means differ -> outlier-driven)")

# ------------------------------- pin section 8 (regress) targets ----------
print("\n=== TARGETS FOR .do SECTION 8 (regress) ===")
X1 = sm.add_constant(u[["high_activity"]].astype(float))
m1 = sm.OLS(u["avg_efficiency"], X1).fit()
m2 = sm.OLS(u["avg_efficiency"], X1).fit(cov_type="HC1")
X3 = sm.add_constant(u[["avg_steps_k"]])
m3 = sm.OLS(u["avg_efficiency"], X3).fit()
X4 = sm.add_constant(u[["avg_steps_k", "avg_veryactive", "avg_sedentary"]])
m4 = sm.OLS(u["avg_efficiency"], X4).fit()

print(f"M1 b(high_activity)={m1.params['high_activity']:+.5f} "
      f"t={m1.tvalues['high_activity']:.3f} p={m1.pvalues['high_activity']:.4f} R2={m1.rsquared:.4f}")
print(f"M2 robust(HC1)      p={m2.pvalues['high_activity']:.4f}  "
      f"(Stata vce(robust) uses HC1 -> should match)")
print(f"M3 b(steps_k)={m3.params['avg_steps_k']:+.6f} "
      f"p={m3.pvalues['avg_steps_k']:.4f} R2={m3.rsquared:.4f}")
print(f"M4 R2={m4.rsquared:.4f}  df_model={int(m4.df_model)} df_resid={int(m4.df_resid)}")

# ------------------------------- pin section 9 (logit) targets ------------
print("\n=== TARGETS FOR .do SECTION 9 (logit) ===")
print(f"good_sleeper: {int(u.good_sleeper.sum())} of {len(u)} "
      f"({u.good_sleeper.mean():.1%}) -> events per predictor is low")
try:
    L1 = sm.Logit(u["good_sleeper"], sm.add_constant(u[["avg_steps_k"]])).fit(disp=0)
    print(f"L1 b(steps_k)={L1.params['avg_steps_k']:+.5f} "
          f"OR={np.exp(L1.params['avg_steps_k']):.4f} "
          f"p={L1.pvalues['avg_steps_k']:.4f} pseudoR2={L1.prsquared:.4f}")
except Exception as e:
    print("L1 failed:", e)
try:
    L2 = sm.Logit(u["good_sleeper"],
                  sm.add_constant(u[["high_activity"]].astype(float))).fit(disp=0)
    print(f"L2 b(high)={L2.params['high_activity']:+.5f} p={L2.pvalues['high_activity']:.4f}")
except Exception as e:
    print("L2 failed:", e)

ct = pd.crosstab(u.good_sleeper, u.high_activity)
print("2x2 table (rows good_sleeper, cols high_activity):\n", ct.to_string())
print(f"Fisher exact p = {stats.fisher_exact(ct.values)[1]:.4f}")

# ------------------------------- pin section 10 (sensitivity) -------------
print("\n=== TARGETS FOR .do SECTION 10 (sensitivity) ===")


def welch_by_steps(steps_series, label, expect=None):
    j = pd.DataFrame({"avg_steps": steps_series}).join(user["avg_efficiency"],
                                                      how="inner").dropna()
    m = j.avg_steps.median()
    a = j.loc[j.avg_steps >= m, "avg_efficiency"]
    b = j.loc[j.avg_steps < m, "avg_efficiency"]
    tt, pp = stats.ttest_ind(a, b, equal_var=False)
    extra = f"  (expect {expect})" if expect else ""
    print(f"{label:<34} n={len(j):>2} median={m:>7.1f} p={pp:.4f} "
          f"d={cohens_d(a, b):+.3f}{extra}")


welch_by_steps(act2.groupby("Id")["TotalSteps"].mean(), "10a tranche-2 activity")
pooled = (pd.concat([act1[["Id", "TotalSteps"]], act2[["Id", "TotalSteps"]]])
          .groupby("Id")["TotalSteps"].mean())
welch_by_steps(pooled, "10b pooled both tranches")

tw, pw = stats.ttest_ind(u.loc[u.high_activity == 1, "wtd_efficiency"],
                         u.loc[u.high_activity == 0, "wtd_efficiency"],
                         equal_var=False)
print(f"{'10c duration-weighted efficiency':<34} p={pw:.4f}")

# effect of the retained duplicate rows
ded = sleep.drop_duplicates()
s2 = (ded.groupby(["Id", "logId"])
      .agg(in_bed=("value", "count"), asleep=("value", lambda v: (v == 1).sum()))
      .reset_index())
s2["efficiency"] = s2["asleep"] / s2["in_bed"]
s2 = s2[s2["in_bed"] >= 60]
e2 = s2.groupby("Id")["efficiency"].mean()
j = act_u.join(e2.rename("eff"), how="inner").dropna(subset=["eff"])
m = j.avg_steps.median()
t2, p2 = stats.ttest_ind(j.loc[j.avg_steps >= m, "eff"],
                         j.loc[j.avg_steps < m, "eff"], equal_var=False)
print(f"{'10d if duplicates WERE dropped':<34} n={len(j):>2} "
      f"sessions={len(s2)} p={p2:.4f}  (baseline retains them: p={p:.4f})")
