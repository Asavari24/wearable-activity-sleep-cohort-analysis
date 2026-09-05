# Stata replication

A second implementation of the cohort analysis in `../analysis.py`, written in
Stata. The question, the cohort definitions, and the tests are unchanged. The
point is to check whether the null result is a property of the data or a
property of pandas/SciPy.

## Files

| file | what it is |
|---|---|
| `01_cohort_analysis.do` | the full pipeline: ingest → curate → test. Self-documenting; read the header block first. |
| `verify_baseline.py` | Python reference run. Confirms the input CSVs reproduce the published README, and pins the targets for every section of the `.do` file. |
| `python_baseline_reference.txt` | saved output of the above — the numbers Stata has to hit. |

## Which data files, exactly

This turned out to matter, and the top-level README does not currently say it.

The published result is produced by the **3.12.16–4.11.16 tranche for BOTH
activity and sleep**:

```
data/raw/dailyActivity_merged_03-12_04-11.csv
data/raw/minuteSleep_merged_03-12_04-11.csv
```

The 4.12.16–5.12.16 `minuteSleep_merged.csv` contains only **459** sessions and
**24** users, so it cannot produce the published 556 sessions / 23 users under
any filter. Pairing tranche-1 sleep with tranche-2 activity gives **p = 0.86**
instead of 0.34 — the same qualitative null, different numbers. The tranche
choice is therefore load-bearing and is documented as decision (2) in the
`.do` header.

The CSVs are not committed here (17 MB, and they are Kaggle's to distribute).
Place them in `data/raw/` using the filenames above.

## Running it

```bash
stata -b do stata/01_cohort_analysis.do
```

Run from the repository root. Output goes to `logs/01_cohort_analysis.log` and
derived `.dta` files to `data/derived/`.

To regenerate the expected numbers first:

```bash
python stata/verify_baseline.py
```

## What the .do file demonstrates

- **append** — the two activity tranches ship from different mirrors with
  different schemas (`participant-id`/unix-epoch vs `Id`/`M/D/YYYY`), so they
  are harmonised, provenance-tagged, stacked, and de-duplicated on the
  overlapping boundary date.
- **merge 1:m** — user-level cohort assignment broadcast back onto the 556
  session rows.
- **reshape long/wide** — minute-level sleep states collapsed to
  session × state counts, pivoted wide to compute efficiency, then
  round-tripped back to long.
- **egen / collapse** — `rowtotal`, `median`, `tag`; collapse to per-user
  aggregates.

and then `summarize`, `ttest`, `regress`, `logit`.

## Verified expectations

Every headline number is pinned in the `.do` header and re-checked inline as
the script runs. The Python reference reproduces the published README on all
16 reported metrics exactly.

| | published | reference run |
|---|---|---|
| users / sessions | 23 / 556 | 23 / 556 |
| median-split threshold | 8135 | 8135 |
| mean efficiency, high / low | 0.898 / 0.930 | 0.898 / 0.930 |
| Welch t / p | −0.98 / 0.34 | −0.98 / 0.34 |
| Cohen's d | −0.41 | −0.41 |
| Mann-Whitney p | 0.13 | 0.13 |
| Levene p | 0.42 | 0.42 |
| CDC split | 5/9, p 0.41, d +0.42 | 5/9, p 0.41, d +0.42 |

## Two things the port surfaced

**Precision.** `participant-id` reaches 8,877,689,391 and `log-id` reaches
11,374,876,178. pandas read both as int64 without being asked. Stata's default
`float` holds ~7 significant digits and would have rounded them, silently
merging distinct users and distinct sleep sessions. Every id is imported as a
string and cast to `double`, and the distinct counts are asserted, because
this failure is invisible in the output but corrupts everything downstream.

**No dose-response.** Regressing efficiency on continuous steps gives
R² = 0.004, *lower* than the cohort-dummy model's 0.042. Dropping the
arbitrary median cut should make a real effect easier to see, not harder.
There is nothing for the split to have discretised.

## Status

The `.do` file has **not yet been executed** — no Stata license was available
on the machine where it was written. It is carefully constructed and its
expected output is fully pinned, but it is unrun code until someone runs it.
Treat any divergence from the table above as a bug in the `.do` file and start
with the double-precision casts in sections 1–2.
