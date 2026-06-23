# FEAST-backed survey budget frontier — handoff

**Goal:** give the EOL / marginal-wells paper a *real FEAST-PtE* basis for the
"how to most efficiently allocate survey money" elbow (Figure 3), replacing the
parametric stub in `ColabMarginalWells/.../figure3_pareto_frontier.py`.

## What was wrong before
Every Bridger script in this folder (`run_bridger_analysis.py`, the "component"
and "bottom-up" variants) used FEAST only as a **data loader** (TMY wind, leak
distribution) and then applied the GML 2.0 POD equation with a Bernoulli
detect/non-detect draw. None of them call FEAST's time-stepping LDAR engine
(`Scenario.run`, `LDARProgram`, `SiteSurvey`/`CompSurvey`, `Repair`). So the
paper could not honestly say "we used FEAST" for the allocation modeling.

## What this adds
`feast_survey_budget_frontier.py` runs the **actual FEAST engine**:

1. Builds a FEAST `GasField` of the PA marginal wells.
2. Encodes four strategies (S1 flat random, S2 tier-stratified, S3 aerial→ground,
   S4 aerial+RS prior) as `LDARProgram`s. Each site-level screen (`SiteSurvey`)
   dispatches flagged sites to a ground follow-up (`CompSurvey`) that localizes
   and repairs components — the plane→OGI→repair pattern from `Example-RunScript.py`.
3. Sweeps a budget grid (budget → wells surveyed via each technique's $/well),
   running a Monte Carlo of `Scenario.run()` per (strategy, budget).
4. Reports per (strategy, budget): **emission reduction %** (vs the auto Null
   baseline, FEAST-native), **realized cost $** (deployment + repair, FEAST-native),
   and **basin emission-rate 95% CI %** (documented coverage/uncertainty estimator).
5. Writes `feast_frontier_results.{csv,json}` and prints the elbow (max curvature
   of CI-vs-budget) per strategy.

The companion `figure3_pareto_frontier_feast.py` (in the Colab repo) consumes that
CSV and renders the frontier + elbow markers and a mitigation-vs-budget plot.

## How to run
```bash
cd <FEAST_PtE root>
pip install -r requirements.txt          # needs numpy, pandas, scipy, matplotlib
# fast wiring check (synthetic ~500-site field, no real data needed):
python MyCodetoRun/feast_survey_budget_frontier.py --smoke
# full run on the real PA marginal wells:
python MyCodetoRun/feast_survey_budget_frontier.py \
    --data MyCodetoRun/feast_emissions.csv --emission-col emission_kgph \
    --n-mc 50 --out ../ColabMarginalWells/WritingPapersClaude/paper_assets/data
# then render Figure 3:
cd ../ColabMarginalWells/WritingPapersClaude/paper_assets/scripts
python figure3_pareto_frontier_feast.py
```
Start with `--n-sites-cap 2000 --n-mc 5` to confirm curve shape, then scale up.

## Validated
The `--smoke` run executes the full FEAST loop in-engine: screens detect sites,
the follow-up CompSurvey localizes components, Repair fires (e.g. 22 repairs,
~93% single-survey mitigation on the synthetic field), cost accrues, and
basin-CI falls monotonically with budget to a detectable elbow. This confirms the
wiring; the **numbers are illustrative until run on the real PA data.**

## Assumptions / knobs to review before quoting numbers
- **Budget sizes the screening campaign**, not total spend; `realized_cost`
  includes follow-up localization + repair on top, so it can exceed the budget
  input. Change in `build_program` if you want budget = total program cost.
- **Aerial strategies (S3/S4) show low mitigation on marginal wells** because
  heli/airborne POD90 ≈ 8–10 kg/h while marginal wells emit <2 kg/h — aerial
  screening is largely blind to them. This is a real finding, not a bug, but
  confirm it matches the instruments you intend to model.
- **`basin_estimate_ci`** is one explicit estimator (measure n wells at the
  technique's `unc_pct`, extrapolate unsurveyed by sample mean, bootstrap). Swap
  in your preferred basin-estimation model if reviewers want a different one.
- **POD curves** are lognormal-CDF anchored at each technique's POD90 from
  `config.py`. `pod_gml2()` (Thorpe 2024 Bridger Combined GML 2.0) is included so
  the aerial curve can be swapped to your exact Bridger model.
- **The real-data GasField loader** scales component count by per-well rate to
  honor the emission distribution while keeping FEAST's native dynamics. For
  exact-truth emissions, construct an `Emission` per site and pass `emissions=`
  to `GasField` (hook noted in `build_pa_gasfield`).

## Files
- `feast_survey_budget_frontier.py` — the FEAST driver (this folder).
- `BridgerResults/frontier/feast_frontier_results.{csv,json}` — outputs.
- `../../ColabMarginalWells/.../scripts/figure3_pareto_frontier_feast.py` — Figure 3.
- `../../ColabMarginalWells/.../data/feast_frontier_results.csv` — currently the
  **smoke-test** CSV; overwrite with the real-data run before using in the paper.
