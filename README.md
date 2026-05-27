# Wearable Activity & Sleep Cohort Analysis

A cohort-based statistical analysis testing whether daily activity level predicts sleep efficiency on consumer wearable telemetry. Built on the Fitbit Fitness Tracker Data (Kaggle, Mobius), 30 users / ~2 months of paired activity + minute-level sleep data.

## Question

Do high-activity users have higher sleep efficiency than low-activity users?

In product language: *if a wearable platform wanted to surface activity-based sleep recommendations to its members, is there a measurable behavioral signal to anchor the feature on?*

## Hypothesis

- **H₀**: There is no difference in mean sleep efficiency between high-activity and low-activity user cohorts.
- **H₁**: Mean sleep efficiency differs between the two cohorts (two-sided).

## Data

- **Source**: Fitbit Fitness Tracker Data (Mobius, Kaggle), MTurk respondents, 3/12/2016–5/12/2016
- **Files used**: `dailyActivity_merged.csv`, `minuteSleep_merged.csv`
- **Users with both signals (after merge)**: 23
- **Sleep sessions ≥ 60 minutes**: 556

**Limitations baked into the dataset**:
- Small, self-selected MTurk sample
- Mix of Fitbit device models (per dataset description)
- No control for age, baseline fitness, or health conditions
- Two-month observation window only

## Methodology

### 1. Sleep efficiency from minute-level data

Sleep efficiency was computed per sleep session as:

```
sleep_efficiency = minutes_asleep / total_minutes_in_bed
```

Where `value == 1` indicates asleep (vs 2 = restless, 3 = awake) in the minute-level file. Sessions shorter than 60 minutes were filtered as likely fragments. Per-user efficiency = mean of that user's session efficiencies.

### 2. Cohort assignment — two splits

Two cohort definitions were used to test the robustness of any observed effect:

- **Median split**: high (≥ 8,135 avg daily steps) vs low (< 8,135). Guarantees balanced groups (n = 12 vs 11).
- **CDC threshold split** (sensitivity check): high (≥ 10,000) vs low (< 7,000). More interpretable but small groups (n = 5 vs 9).

### 3. Statistical tests

- **Assumption checks**: Shapiro-Wilk (normality), Levene (equal variance)
- **Primary**: Welch's two-sample t-test (robust to unequal variance)
- **Non-parametric backup**: Mann-Whitney U (used when normality assumption fails)
- **Effect size**: Cohen's d
- **Power analysis**: required n per cohort to detect observed d at 80% power, α = 0.05

## Results

### Median split (n = 12 high, n = 11 low)

| Metric | Value |
|---|---|
| High-activity mean efficiency | 0.898 |
| Low-activity mean efficiency | 0.930 |
| High-activity median efficiency | 0.929 |
| Low-activity median efficiency | 0.947 |
| Welch's t | -0.98 |
| **p-value (Welch's)** | **0.34** |
| **p-value (Mann-Whitney)** | **0.13** |
| Cohen's d | -0.41 (medium) |
| Shapiro p (high) | 0.001 ← non-normal |
| Shapiro p (low) | <0.001 ← non-normal |
| Levene p | 0.42 (variances equal) |

### CDC split (n = 5 high, n = 9 low) — sensitivity check

| Metric | Value |
|---|---|
| High-activity mean efficiency | 0.949 |
| Low-activity mean efficiency | 0.926 |
| Welch's t | 0.85 |
| p-value | 0.41 |
| Cohen's d | +0.42 (direction flipped) |

### Power analysis

To detect Cohen's d = 0.4 with 80% power at α = 0.05 in a two-sample design, the required sample size is **~96 users per cohort (~192 total)**. The available data (23 users total, 11–12 per cohort) is severely underpowered for the observed effect.

## Interpretation

**We failed to reject the null hypothesis.** There is no statistically significant evidence in this sample that high-activity users have different sleep efficiency than low-activity users.

Three findings are worth surfacing:

1. **Effect-size direction flipped between cohort definitions** (d = -0.41 under median split, +0.42 under CDC split). This indicates either a non-monotonic relationship between activity and sleep efficiency or, more likely, instability driven by small sample sizes.
2. **Means were skewed by outliers** (2-3 users with very low sleep efficiency). The medians of the two cohorts were nearly identical (0.929 vs 0.947), suggesting that on the typical user, the cohorts behave similarly.
3. **The study is severely underpowered.** Detecting a medium-magnitude effect would require roughly 8× more users with paired signals than this dataset contains.

## Product recommendations

If a product analytics team wanted to act on this question:

- **Do not ship a feature** that assumes activity-driven sleep differences based on this dataset.
- **Collect more longitudinal data per user** — within-user comparisons (does *this user* sleep better after high-activity days?) are statistically more powerful than between-user comparisons and avoid confounding by baseline health.
- **Control for confounders** — age, baseline fitness, sleep latency, and device model are likely material.
- **Consider segmenting differently** — sleep consistency, sleep duration, or recovery proxies may show stronger activity relationships than efficiency.

## Files

- `analysis.py` — full analytic pipeline (load → roll up → cohort → test → visualize)
- `cohort_analysis.png` — per-user scatter + cohort boxplot
- `README.md` — this file

## Tech

Python, Pandas, SciPy, Matplotlib, Seaborn.
