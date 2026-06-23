"""
feast_survey_budget_frontier.py
===============================================================================
A REAL FEAST-PtE simulation driver that produces the cost-vs-uncertainty
"elbow" (Pareto frontier) used in the EOL / marginal-wells paper.

WHY THIS EXISTS
---------------
The earlier Bridger scripts in this folder (run_bridger_analysis.py, etc.) used
FEAST only as a data loader and then applied the GML 2.0 POD equation with a
Bernoulli detect/non-detect draw.  They did NOT run FEAST's time-stepping LDAR
engine, so the paper cannot legitimately say "we used FEAST" for the
survey-money-allocation modeling.

THIS script does run the engine:  it builds a FEAST GasField of the PA marginal
wells, defines four survey STRATEGIES (S1-S4) as FEAST LDARProgram objects, and
for each strategy x budget level runs a Monte Carlo of feast Scenario.run().
From those runs it reports, per (strategy, budget):

    * emission-reduction (mitigation) % vs the auto-generated Null baseline   <- FEAST-native
    * realized survey cost ($)                                               <- FEAST-native (deployment + repair)
    * basin emission-rate 95% CI (%)                                         <- coverage/uncertainty estimator
                                                                                layered on the FEAST emission truth

The output CSV/JSON is designed to REPLACE the parametric stub in
ColabMarginalWells/WritingPapersClaude/paper_assets/scripts/figure3_pareto_frontier.py
and config.py (PARETO_PARAMS / STRATEGIES.ci_pct_at_500k).

HOW TO RUN
----------
    cd /workspaces/FEAST_PtE         # (or your local FEAST_PtE root)
    source .venv/bin/activate
    # Smoke test (synthetic ~60-site field, runs in the sandbox, validates wiring):
    python MyCodetoRun/feast_survey_budget_frontier.py --smoke
    # Full run on the real PA marginal wells (needs feast_emissions.csv on disk):
    python MyCodetoRun/feast_survey_budget_frontier.py --data MyCodetoRun/feast_emissions.csv

DATA EXPECTATIONS
-----------------
--data points to a CSV with one row per well and at least a site-total emission
column in kg/h (default column name: 'emission_kgph'; override with --emission-col).
An optional 'tier' or emission-percentile is derived internally for stratification.
If --data is omitted or missing, the script falls back to a synthetic GasField
bootstrapped from ExampleData/DataObjectInstances/production_emissions.p so the
FEAST wiring can still be exercised.

CAVEATS (read before quoting numbers in the paper)
--------------------------------------------------
* The basin-CI estimator is an explicit, documented model (see basin_estimate_ci);
  it is NOT the only defensible choice. It uses each technique's 1-sigma
  measurement uncertainty (unc_pct) and the realized survey coverage from the
  FEAST run. Swap in your preferred estimator if the reviewers want a different one.
* Full PA scale (64,624 sites x budget grid x strategies x MC) is heavy. Start
  with --n-sites-cap and --n-mc small, confirm the curve shape, then scale up.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import os
import sys
from pathlib import Path

import numpy as np

# --- make 'import feast' work no matter where the script is called from ------
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import feast.EmissionSimModules.simulation_classes as sc
import feast.DetectionModules as Dm
from feast.EmissionSimModules.infrastructure_classes import Component, Site, GasField

KGH_TO_GPS = 1000.0 / 3600.0          # kg/h  ->  g/s
GPS_TO_KGH = 3600.0 / 1000.0          # g/s   ->  kg/h


# =============================================================================
# 1. PARAMETER TABLES  (mirror of ColabMarginalWells/paper_assets/scripts/config.py)
#    Keep these in sync with the paper's config.py.  Costs are $/well, pod is the
#    90% probability-of-detection emission rate (kg/h), throughput is wells/day.
# =============================================================================

TECHNIQUES = {
    "chamber":  {"name": "Static flux chambers", "tier": "ground",    "cost": 150, "pod90_kgph": 0.0001, "unc_pct": 15, "throughput": 10},
    "hfs":      {"name": "High-flow sampling",    "tier": "ground",    "cost": 300, "pod90_kgph": 0.001,  "unc_pct": 20, "throughput": 12},
    "vehicle":  {"name": "Vehicle downwind",      "tier": "ground",    "cost": 100, "pod90_kgph": 0.2,    "unc_pct": 50, "throughput": 30},
    "drone":    {"name": "Drone flux plane",      "tier": "ground",    "cost": 550, "pod90_kgph": 0.1,    "unc_pct": 32, "throughput": 8},
    "heli":     {"name": "Helicopter LiDAR",      "tier": "aerial",    "cost": 100, "pod90_kgph": 10.0,   "unc_pct": 47, "throughput": 300},
    "airimg":   {"name": "Airborne imaging",      "tier": "aerial",    "cost": 12,  "pod90_kgph": 8.0,    "unc_pct": 40, "throughput": 800},
    "msat":     {"name": "MethaneSAT",            "tier": "satellite", "cost": 0.0, "pod90_kgph": 100.0,  "unc_pct": 30, "throughput": 100000},
}

# Default budget grid (USD).  Mirror of config.BUDGET_GRID_USD.
BUDGET_GRID_USD = [50_000, 100_000, 200_000, 350_000, 500_000, 750_000,
                   1_000_000, 1_500_000, 2_000_000]

# Strategy definitions.  Each strategy is an ordered list of "stages".  A stage
# is a technique plus a budget fraction and a queue-ordering rule.  Multi-stage
# strategies (S3, S4) dispatch flagged sites from an aerial screen to a ground
# follow-up, exactly like the plane->ogi pattern in Example-RunScript.py.
#   order: "random"   -> survey a random subset (no targeting)
#          "emitter"  -> survey highest-emitting wells first (perfect prior; upper bound)
#          "rs_prior" -> survey by a noisy remote-sensing prior (realistic targeting)
STRATEGIES = {
    "S1": {"name": "Flat random",
           "stages": [{"tech": "vehicle", "frac": 1.0, "order": "random"}]},
    "S2": {"name": "Tier-stratified",
           "stages": [{"tech": "heli",    "frac": 0.5, "order": "random"},
                      {"tech": "vehicle", "frac": 0.5, "order": "random"}]},
    "S3": {"name": "Aerial then ground",
           "stages": [{"tech": "heli",    "frac": 0.6, "order": "random",  "dispatch": "vehicle"}]},
    "S4": {"name": "Aerial + RS prior",
           "stages": [{"tech": "heli",    "frac": 0.6, "order": "rs_prior", "dispatch": "vehicle"}]},
}

# Simulation settings
SIM_DAYS_DEFAULT = 180          # one survey campaign / season
DELTA_T = 1                     # day
RS_PRIOR_NOISE = 0.6            # lognormal sigma on the remote-sensing prior ranking (S4)


# =============================================================================
# 2. POD SURFACE BUILDER
#    SiteSurvey wants detection_probability_points (g/s) + detection_probabilities.
#    We build a lognormal-CDF POD curve in flux anchored so that POD(pod90)=0.90.
#    (For the aerial LiDAR you can swap in the Bridger GML 2.0 curve via
#     pod_gml2_curve(); kept here so the link to your Bridger work is explicit.)
# =============================================================================

def lognormal_pod_curve(pod90_kgph: float, sigma: float = 1.0, n_pts: int = 120):
    """Return (points_gps, probs) for a lognormal-CDF POD anchored at POD90."""
    # mu chosen so that CDF(log pod90) = 0.90  ->  log(pod90) = mu + z90*sigma
    z90 = 1.2815515594600704  # Phi^-1(0.9)
    mu = math.log(pod90_kgph) - z90 * sigma
    pts_kgph = np.logspace(math.log10(pod90_kgph) - 3, math.log10(pod90_kgph) + 2, n_pts)
    probs = 0.5 + 0.5 * np.array([math.erf((math.log(f) - mu) / (sigma * math.sqrt(2)))
                                  for f in pts_kgph])
    probs[0] = 0.0
    return pts_kgph * KGH_TO_GPS, np.clip(probs, 0.0, 1.0)


# --- Bridger GML 2.0 POD (Thorpe et al. 2024, Table 3) — for explicit linkage --
_GML = dict(a1=2.0, a2=1.5, b1=2.41e-3, b2=1.9505, b3=2.0836, b4=1.5185)
def pod_gml2(emission_kgph, wind_ms, n=0.016):
    Q = np.asarray(emission_kgph, float); u = float(np.clip(wind_ms, 1.0, 6.0))
    g = _GML["b1"] * Q**_GML["b2"] / (n**_GML["b3"] * u**_GML["b4"])
    return np.clip(1.0 - (1.0 + g**_GML["a1"])**(-_GML["a2"]), 0.0, 1.0)

def pod_gml2_curve(wind_ms=3.5, n_pts=120):
    pts_kgph = np.logspace(-3, 2, n_pts)
    probs = pod_gml2(pts_kgph, wind_ms)
    probs[0] = 0.0
    return pts_kgph * KGH_TO_GPS, probs


# Near-field ground-localization POD (OGI-style) used by the follow-up CompSurvey
# that actually finds & repairs individual components after a site-level screen
# flags a site. Mirrors ogi_no_survey in Example-RunScript.py.
def ogi_followup_curve(n_pts=100):
    pts_gps = np.logspace(-3, 1, n_pts)             # component flux, g/s
    probs = 0.5 + 0.5 * np.array([math.erf((math.log(f) - math.log(0.02)) / (0.8 * math.sqrt(2)))
                                  for f in pts_gps])
    probs[0] = 0.0
    return pts_gps, np.clip(probs, 0.0, 1.0)


# =============================================================================
# 3. GAS FIELD LOADERS
# =============================================================================

def _component_from_leakdist(reparable=True):
    """A FEAST Component that bootstraps leak sizes from the bundled distribution."""
    return Component(
        name="well_emitter",
        emission_data_path="ExampleData/DataObjectInstances/production_emissions.p",
        emission_per_comp=0.0231,                 # expected emitting fraction at t0
        emission_production_rate=5.4 / 650 / 365, # new emissions / comp / day
        repair_cost_path="ExampleData/DataObjectInstances/fernandez_leak_repair_costs_2006.p",
        base_reparable=reparable,
    )


def build_synthetic_gasfield(time, n_sites=60):
    """Smoke-test field: n_sites single-component wells bootstrapped from the
    bundled FEAST leak-size distribution.  Guaranteed to run in the sandbox."""
    comp = _component_from_leakdist()
    well = Site(
        name="marginal_well",
        comp_dict={"emitter": {"number": 40, "parameters": copy.copy(comp)}},
    )
    site_dict = {"marginal_well": {"number": n_sites, "parameters": well}}
    gf = GasField(sites=site_dict, time=time)
    gf.met_data_path = "ExampleData/TMY-DataExample.csv"
    gf.met_data_maker()
    return gf


def build_pa_gasfield(time, data_path, emission_col="emission_kgph", n_sites_cap=None):
    """Build a GasField from the real PA marginal wells emission file.

    Each well becomes one FEAST Site.  Per-well emission magnitude is honored by
    scaling the number of bootstrap components so the site's expected emission
    matches the data (a pragmatic way to inject the real emission distribution
    while keeping FEAST's native emission dynamics).  For an exact-truth variant,
    construct an Emission object per site and pass emissions=... to GasField.
    """
    import pandas as pd
    df = pd.read_csv(data_path)
    if emission_col not in df.columns:
        raise ValueError(f"Column '{emission_col}' not in {data_path}. "
                         f"Available: {list(df.columns)[:12]}... use --emission-col")
    rates = df[emission_col].to_numpy(float)
    rates = rates[np.isfinite(rates) & (rates >= 0)]
    if n_sites_cap:
        rates = rates[:n_sites_cap]

    # Group wells into emission deciles; one representative Site per well but with
    # component count proportional to its emission rate (keeps the field tractable).
    base_comp = _component_from_leakdist()
    median_rate = max(np.median(rates[rates > 0]) if np.any(rates > 0) else 1.0, 1e-3)

    site_dict = {}
    for i, r in enumerate(rates):
        n_comp = int(max(1, round(40 * r / median_rate)))   # scale emission load by rate
        well = Site(name=f"well_{i}",
                    comp_dict={"emitter": {"number": n_comp, "parameters": copy.copy(base_comp)}})
        site_dict[f"well_{i}"] = {"number": 1, "parameters": well}

    gf = GasField(sites=site_dict, time=time)
    gf.met_data_path = "ExampleData/TMY-DataExample.csv"
    gf.met_data_maker()
    gf._well_rates_kgph = rates          # stash truth for the basin-CI estimator
    return gf


# =============================================================================
# 4. STRATEGY -> LDAR PROGRAM
# =============================================================================

def _site_queue_for_order(order, n_sites, n_survey, well_rates=None, rng=None):
    """Return the ordered list of site indices to survey under a queue-ordering rule."""
    rng = rng or np.random
    n_survey = int(min(n_survey, n_sites))
    if order == "random" or well_rates is None:
        idx = rng.permutation(n_sites)[:n_survey]
    elif order == "emitter":                       # perfect prior (upper bound)
        idx = np.argsort(well_rates)[::-1][:n_survey]
    elif order == "rs_prior":                      # noisy remote-sensing prior
        noisy = np.log(np.maximum(well_rates, 1e-6)) + rng.normal(0, RS_PRIOR_NOISE, n_sites)
        idx = np.argsort(noisy)[::-1][:n_survey]
    else:
        idx = rng.permutation(n_sites)[:n_survey]
    return list(np.array(idx, dtype=int))


def build_program(gas_field, time, strategy, budget_usd, well_rates=None, rng=None):
    """Compose one FEAST LDARProgram for a strategy at a given budget.

    Budget -> number of wells surveyed per stage via that technique's $/well.
    Coverage is enforced by pre-loading site_queue (drained at sites_per_day);
    survey_interval is left None so each well is surveyed at most once per campaign.
    """
    n_sites = gas_field.n_sites
    tech_dict = {}
    repair = Dm.repair.Repair(repair_delay=0)

    for s_i, stage in enumerate(strategy["stages"]):
        t = TECHNIQUES[stage["tech"]]
        stage_budget = budget_usd * stage["frac"]
        n_survey = int(stage_budget // max(t["cost"], 1e-9)) if t["cost"] > 0 else n_sites
        n_survey = max(1, min(n_survey, n_sites))
        queue = _site_queue_for_order(stage["order"], n_sites, n_survey, well_rates, rng)

        pts, probs = lognormal_pod_curve(t["pod90_kgph"])

        # Follow-up ground localization that actually FINDS & REPAIRS flagged
        # emissions. A site-level screen only flags a *site*; FEAST needs a
        # comp-level survey to resolve emission IDs and dispatch them to Repair
        # (this is the plane -> ogi -> repair pattern from Example-RunScript.py).
        gt = TECHNIQUES[stage.get("dispatch", "vehicle")]
        f_pts, f_probs = ogi_followup_curve()
        followup = Dm.comp_survey.CompSurvey(
            time, survey_interval=None, survey_speed=150,
            ophrs={"begin": 8, "end": 17}, labor=100,
            detection_variables={"flux": "mean"},
            detection_probability_points=f_pts,
            detection_probabilities=f_probs,
            dispatch_object=copy.deepcopy(repair), site_queue=[],
            sensitivity=0.1, dispatch_threshold=None)   # None -> repair every component it localizes

        survey = Dm.site_survey.SiteSurvey(
            time, survey_interval=None, sites_per_day=t["throughput"], site_cost=t["cost"],
            detection_variables={"flux": "mean"},
            detection_probability_points=pts, detection_probabilities=probs,
            dispatch_object=followup, site_queue=queue,
            ophrs={"begin": 8, "end": 17}, sensitivity=0.1,
            dispatch_threshold=t["pod90_kgph"] * KGH_TO_GPS)

        tech_dict[f"{strategy['name']}_s{s_i}_{stage['tech']}"] = survey
        tech_dict[f"{strategy['name']}_s{s_i}_{gt['name']}_followup"] = followup

    return Dm.ldar_program.LDARProgram(copy.deepcopy(gas_field), tech_dict)


# =============================================================================
# 5. METRICS
# =============================================================================

def program_emission_mean(prog):
    ts = np.array(prog.emissions_timeseries, dtype=float)
    return float(np.mean(ts)) if ts.size else float("nan")


def _sum_result(res):
    """Sum a FEAST ResultDiscrete's values (stored as [t, value] pairs in time_value)."""
    if res is None or not getattr(res, "time_value", None):
        return 0.0
    try:
        return float(res.get_sum_val())
    except Exception:
        return float(np.sum([tv[1] for tv in res.time_value]))


def program_cost(prog):
    cost = 0.0
    for tech in prog.tech_dict.values():
        cost += _sum_result(getattr(tech, "deployment_cost", None))
    cost += _sum_result(getattr(prog, "repair_cost", None))
    return cost


def basin_estimate_ci(well_rates_kgph, n_surveyed, unc_pct, n_boot=400, rng=None):
    """Documented basin emission-rate 95% CI estimator (the Figure-3 y-axis).

    Model: a survey measures n_surveyed wells with 1-sigma relative error unc_pct;
    unsurveyed wells are estimated by the surveyed-sample mean. The basin total is
    the sum. Bootstrapping which wells are surveyed + measurement noise gives a CI.
    Returns the 95% CI half-width as a percentage of the mean estimate.
    """
    rng = rng or np.random
    rates = np.asarray(well_rates_kgph, float)
    N = len(rates); n_surveyed = int(min(max(n_surveyed, 1), N))
    true_total = rates.sum()
    ests = np.empty(n_boot)
    for b in range(n_boot):
        sel = rng.choice(N, size=n_surveyed, replace=False)
        meas = rates[sel] * (1 + rng.normal(0, unc_pct / 100.0, n_surveyed))
        sample_mean = max(meas.mean(), 0.0)
        ests[b] = meas.sum() + sample_mean * (N - n_surveyed)
    mean_est = ests.mean()
    ci_halfwidth = 1.96 * ests.std(ddof=1)
    return 100.0 * ci_halfwidth / mean_est if mean_est > 0 else float("nan"), true_total


# =============================================================================
# 6. SWEEP
# =============================================================================

def run_sweep(args):
    results = []
    time = sc.Time(delta_t=DELTA_T, end_time=args.sim_days)

    use_real = args.data and Path(args.data).exists()
    if use_real:
        print(f"[data] real PA wells: {args.data}")
    else:
        if args.data:
            print(f"[data] {args.data} not found -> synthetic smoke field")
        else:
            print("[data] no --data -> synthetic smoke field")

    for sid, strat in STRATEGIES.items():
        for budget in args.budget_grid:
            mit, cost, basin_ci = [], [], []
            for mc in range(args.n_mc):
                rng = np.random.RandomState(1000 * (mc + 1) + hash(sid) % 997)
                t = sc.Time(delta_t=DELTA_T, end_time=args.sim_days)
                if use_real:
                    gf = build_pa_gasfield(t, args.data, args.emission_col, args.n_sites_cap)
                    well_rates = gf._well_rates_kgph
                else:
                    gf = build_synthetic_gasfield(t, n_sites=args.n_sites_cap or 60)
                    # derive a per-site rate proxy from the generated emissions for ordering/CI
                    em = gf.emissions.emissions
                    well_rates = (em.groupby("site_index")["flux"].sum()
                                  .reindex(range(gf.n_sites)).fillna(0).to_numpy() * GPS_TO_KGH)

                prog = build_program(gf, t, strat, budget, well_rates=well_rates, rng=rng)
                scenario = sc.Scenario(time=t, gas_field=gf,
                                       ldar_program_dict={sid: prog})
                progs = scenario.run(display_status=False, save_method="object")

                null_em = program_emission_mean(progs["Null"])
                prog_em = program_emission_mean(progs[sid])
                mit.append(100.0 * (null_em - prog_em) / null_em if null_em > 0 else 0.0)
                cost.append(program_cost(progs[sid]))

                # coverage actually realized -> basin CI
                first_tech = list(prog.tech_dict.values())[0]
                n_surveyed = int(_sum_result(getattr(first_tech, "deployment_count", None)))
                ci, _ = basin_estimate_ci(well_rates, max(n_surveyed, 1),
                                          TECHNIQUES[strat["stages"][0]["tech"]]["unc_pct"], rng=rng)
                basin_ci.append(ci)

            row = dict(strategy=sid, strategy_name=strat["name"], budget_usd=budget,
                       mitigation_pct=float(np.mean(mit)),
                       mitigation_pct_sd=float(np.std(mit)),
                       basin_ci_pct=float(np.nanmean(basin_ci)),
                       realized_cost_usd=float(np.mean(cost)), n_mc=args.n_mc)
            results.append(row)
            print(f"  {sid} {strat['name']:<20} ${budget:>9,}  "
                  f"mitig={row['mitigation_pct']:5.1f}%  "
                  f"basinCI={row['basin_ci_pct']:5.1f}%  "
                  f"cost=${row['realized_cost_usd']:>10,.0f}")
    return results


def detect_elbow(rows_for_strategy):
    """Knee of basin_ci vs budget: max curvature on the log-budget axis."""
    rows = sorted(rows_for_strategy, key=lambda r: r["budget_usd"])
    x = np.log10([r["budget_usd"] for r in rows]); y = np.array([r["basin_ci_pct"] for r in rows])
    if len(x) < 3:
        return None
    d2 = np.gradient(np.gradient(y, x), x)
    return rows[int(np.argmax(np.abs(d2)))]["budget_usd"]


def save_results(results, out_dir):
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "feast_frontier_results.json"
    csv_path = out_dir / "feast_frontier_results.csv"
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2)
    cols = ["strategy", "strategy_name", "budget_usd", "mitigation_pct",
            "mitigation_pct_sd", "basin_ci_pct", "realized_cost_usd", "n_mc"]
    with open(csv_path, "w") as f:
        f.write(",".join(cols) + "\n")
        for r in results:
            f.write(",".join(str(r[c]) for c in cols) + "\n")
    # elbow summary
    print("\n=== ELBOW (knee of basin-CI vs budget) ===")
    for sid in STRATEGIES:
        rs = [r for r in results if r["strategy"] == sid]
        elbow = detect_elbow(rs)
        if elbow:
            print(f"  {sid} {STRATEGIES[sid]['name']:<20} elbow ~ ${elbow:,}")
    print(f"\n[out] {csv_path}\n[out] {json_path}")
    return csv_path, json_path


# =============================================================================
# 7. CLI
# =============================================================================

def main():
    p = argparse.ArgumentParser(description="FEAST survey budget Pareto frontier")
    p.add_argument("--data", default=None, help="PA wells emission CSV (kg/h column)")
    p.add_argument("--emission-col", default="emission_kgph")
    p.add_argument("--n-mc", type=int, default=20, help="Monte Carlo runs per (strategy,budget)")
    p.add_argument("--sim-days", type=int, default=SIM_DAYS_DEFAULT)
    p.add_argument("--n-sites-cap", type=int, default=None, help="cap number of wells (speed)")
    p.add_argument("--budget-grid", default=None, help="comma-separated USD list (overrides default)")
    p.add_argument("--out", default=None, help="output dir (default: BridgerResults/frontier)")
    p.add_argument("--smoke", action="store_true", help="fast synthetic smoke test")
    args = p.parse_args()

    if args.smoke:
        args.n_mc = args.n_mc if args.n_mc != 20 else 2
        args.sim_days = 60
        args.n_sites_cap = args.n_sites_cap or 500
        # budgets chosen so coverage scales below saturation (500 wells) -> visible elbow
        args.budget_grid = args.budget_grid or "3000,8000,20000,50000,120000"

    args.budget_grid = ([int(x) for x in args.budget_grid.split(",")]
                        if args.budget_grid else BUDGET_GRID_USD)
    args.out = args.out or str(Path(__file__).parent / "BridgerResults" / "frontier")

    print("=" * 70)
    print("FEAST-PtE survey budget frontier")
    print(f"  strategies={list(STRATEGIES)}  budgets={args.budget_grid}")
    print(f"  n_mc={args.n_mc}  sim_days={args.sim_days}  n_sites_cap={args.n_sites_cap}")
    print("=" * 70)

    results = run_sweep(args)
    save_results(results, args.out)


if __name__ == "__main__":
    main()
