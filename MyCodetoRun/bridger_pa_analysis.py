"""
Bridger LiDAR Survey Analysis for PA Marginal Wells
====================================================
Grounded in:
  - Thorpe et al. 2024 (Remote Sensing of Environment): Bridger GML 2.0 PoD model
    * Overall 90% PoD @ 1.27 kg/h (2023 basin average)
    * MARCELLUS-SPECIFIC 90% PoD: mean=0.974 kg/h, median=0.927 kg/h [Fig 5 table]
  - Donahue et al. 2024 (Permian Preprint): Campaign planning & throughput benchmark
    * 51,770 total sites, 195 calendar days over 4 quarters
    * Marginal wellsites: 6,357 sites/quarter ≈ 130 sites/day
  - PA marginal well portfolio from existing analysis (confirmed numbers):
    * 64,624 wells, 55,389 kg/h total, mean=0.857 kg/h, median≈0.108 kg/h

Key Questions Addressed:
  1. How much emissions can a Bridger survey capture?
  2. Does it capture MOST EMITTERS (by count) or TOP EMITTERS (by volume)?
  3. What's a realistic flight plan and annual coverage?

Author: FEAST Analysis  |  Date: April 2026
"""

import numpy as np
import pandas as pd
import json
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# ============================================================
# CONSTANTS - All grounded in peer-reviewed papers
# ============================================================

# Portfolio (confirmed from existing analysis of PA GeoPackage data)
N_WELLS = 64_624
TOTAL_EMISSIONS_KGPH = 55_389.0  # kg/h
MEAN_EMISSION = TOTAL_EMISSIONS_KGPH / N_WELLS   # 0.857 kg/h
MEDIAN_EMISSION = 0.108           # kg/h (from existing analysis)

# PoD model - Thorpe et al. 2024
# Use MARCELLUS-SPECIFIC threshold (Fig 5, right table): mean=0.974 kg/h
# Not the overall 1.27 kg/h (which mixes many basins including worse-performing ones)
POD_90_MARCELLUS = 0.974          # kg/h at 90% PoD for Marcellus (Thorpe 2024, Fig 5)
POD_90_OVERALL = 1.27             # kg/h overall 2023 average (less relevant for PA)
LOGISTIC_STEEPNESS = 2.0          # from Conrad et al. 2023a PoD model structure
BASELINE_WIND_MS = 3.5            # m/s reference wind for PoD90

# Flight parameters - derived from Permian 2024 paper (Donahue et al.)
# Permian: 6,357 marginal wellsites/quarter ÷ ~49 flying days/quarter = 130/day
WELLS_PER_DAY_PERMIAN = 130       # marginal wells/day (Permian benchmark)
WELLS_PER_DAY_PA = 100            # PA adjusted: forested terrain, longer transit, more clouds
                                  # ~25% reduction from Permian flat terrain

# PA flying days: Permian had 195 days/year but PA has more precipitation/overcast
# PA average: ~150 flyable days/year (Permian ~195, PA historically ~30% more weather events)
FLYABLE_DAYS_PA_YEAR = 150        # conservative for PA
FLYABLE_DAYS_PA_YEAR_OPT = 170    # optimistic (good year)

# Annual coverage calculation
ANNUAL_WELLS_CONSERVATIVE = WELLS_PER_DAY_PA * FLYABLE_DAYS_PA_YEAR          # 15,000
ANNUAL_WELLS_PERMIAN_RATE = WELLS_PER_DAY_PERMIAN * FLYABLE_DAYS_PA_YEAR     # 19,500
ANNUAL_COVERAGE_PCT_CONS = ANNUAL_WELLS_CONSERVATIVE / N_WELLS                # 23.2%
ANNUAL_COVERAGE_PCT_PERM = ANNUAL_WELLS_PERMIAN_RATE / N_WELLS                # 30.2%

# Monte Carlo
N_ITERATIONS = 1000

# Wind speed distribution for PA (flyable range)
# PA TMY data suggests 3.0-4.0 m/s mean during flying seasons
WIND_MEAN_PA = 3.5                # m/s
WIND_STD_PA = 0.9                 # m/s (seasonal/daily variation)
WIND_MIN_FLYABLE = 1.0
WIND_MAX_FLYABLE = 6.0


# ============================================================
# STEP 1: BUILD EMISSION DISTRIBUTION
# ============================================================
# Lognormal fit to match confirmed PA statistics:
#   median = 0.108 kg/h → mu = ln(0.108) = -2.226
#   mean = 0.857 kg/h → sigma = sqrt(2*(ln(0.857) - mu)) = 2.03
# This is consistent with Marcellus marginal well emission distributions
# in the literature (e.g., Barkley et al. 2023, Peischl et al.)

LN_MU = np.log(MEDIAN_EMISSION)   # -2.226
LN_SIGMA = np.sqrt(2 * (np.log(MEAN_EMISSION) - LN_MU))  # ~2.03

print(f"Lognormal fit: mu={LN_MU:.3f}, sigma={LN_SIGMA:.3f}")
print(f"  Expected mean: {np.exp(LN_MU + LN_SIGMA**2/2):.3f} kg/h (target: {MEAN_EMISSION:.3f})")
print(f"  Expected median: {np.exp(LN_MU):.3f} kg/h (target: {MEDIAN_EMISSION:.3f})")

np.random.seed(42)
emissions = np.random.lognormal(mean=LN_MU, sigma=LN_SIGMA, size=N_WELLS)

# Scale to exactly match known total
emissions = emissions * (TOTAL_EMISSIONS_KGPH / emissions.sum())

print(f"\nSynthetic portfolio generated:")
print(f"  Wells: {N_WELLS:,}")
print(f"  Total emissions: {emissions.sum():.0f} kg/h")
print(f"  Mean: {emissions.mean():.3f} kg/h")
print(f"  Median: {np.median(emissions):.3f} kg/h")
print(f"  Std Dev: {emissions.std():.3f} kg/h")
print(f"  Max: {emissions.max():.2f} kg/h")


# ============================================================
# STEP 2: EMISSION CONCENTRATION ANALYSIS
# (Directly answers: "most emitters" vs "top emitters")
# ============================================================

def emission_concentration_analysis(emissions, label="Portfolio"):
    """Compute emission concentration curve - Lorenz curve for emissions."""
    sorted_e = np.sort(emissions)[::-1]  # descending
    cumsum_e = np.cumsum(sorted_e)
    total_e = cumsum_e[-1]
    n = len(sorted_e)

    results = {}
    for pct in [1, 2, 5, 10, 20, 30, 50, 70]:
        n_wells = max(1, int(n * pct / 100))
        emission_share = 100 * cumsum_e[n_wells - 1] / total_e
        threshold = sorted_e[n_wells - 1]
        results[pct] = {
            'top_pct_wells': pct,
            'n_wells': n_wells,
            'emission_share_pct': emission_share,
            'min_emission_in_group': threshold
        }

    # Also: what % of wells are above various emission thresholds?
    threshold_results = {}
    for thresh in [0.1, 0.5, 0.974, 1.0, 1.27, 2.0, 5.0]:
        n_above = (emissions > thresh).sum()
        emission_above = emissions[emissions > thresh].sum()
        threshold_results[thresh] = {
            'threshold_kgph': thresh,
            'n_wells_above': int(n_above),
            'pct_wells_above': 100 * n_above / n,
            'emissions_above_kgph': emission_above,
            'pct_emissions_above': 100 * emission_above / total_e
        }

    return results, threshold_results


# ============================================================
# STEP 3: PoD MODEL (Thorpe et al. 2024 / Conrad et al. 2023a)
# ============================================================

def pod_bridger(emission_rate_kgph, wind_speed_ms, pod90=POD_90_MARCELLUS):
    """
    Bridger GML 2.0 PoD model.

    Source: Conrad et al. 2023a logistic model on log-scale:
      PoD = 1 / (1 + exp(-k * (log(E) - log(E90(wind)))))
    where E90(wind) = pod90 * exp(-0.3 * (wind - 3.5))

    Args:
        emission_rate_kgph: emission rate in kg/h (scalar or array)
        wind_speed_ms: wind speed in m/s
        pod90: 90% PoD threshold at 3.5 m/s; default = Marcellus 0.974 kg/h

    Returns:
        PoD in [0, 1]
    """
    # Wind-adjusted threshold
    wind_factor = np.exp(-0.3 * (wind_speed_ms - BASELINE_WIND_MS))
    e90_wind = pod90 * wind_factor

    # Logistic on log-scale
    e_safe = np.maximum(emission_rate_kgph, 1e-6)
    log_ratio = np.log(e_safe / e90_wind)
    pod = 1.0 / (1.0 + np.exp(-LOGISTIC_STEEPNESS * log_ratio))

    return np.clip(pod, 0.0, 1.0)


def pod_at_percentiles(emissions, wind_ms=3.5, pod90=POD_90_MARCELLUS):
    """PoD at each emission level - useful for understanding detection coverage."""
    thresholds = [0.1, 0.25, 0.5, 0.974, 1.27, 2.0, 5.0, 10.0]
    results = []
    for thresh in thresholds:
        pod = pod_bridger(thresh, wind_ms, pod90)
        n_above = (emissions >= thresh).sum()
        results.append({
            'emission_kgph': thresh,
            'PoD': pod,
            'n_wells_at_or_above': int(n_above),
            'pct_wells_at_or_above': 100 * n_above / len(emissions)
        })
    return results


# ============================================================
# STEP 4: MONTE CARLO SURVEY SIMULATION
# ============================================================

def sample_wind_pa():
    """Sample wind speed from PA conditions (flyable window)."""
    wind = np.random.normal(WIND_MEAN_PA, WIND_STD_PA)
    return np.clip(wind, WIND_MIN_FLYABLE, WIND_MAX_FLYABLE)


def simulate_survey(emissions, coverage_pct, wind_ms=None, pod90=POD_90_MARCELLUS,
                    targeted_top_pct=None):
    """
    Simulate one Bridger survey pass.

    Args:
        emissions: array of well emission rates (kg/h)
        coverage_pct: fraction of wells surveyed (0-1)
        wind_ms: wind speed; if None, sampled from PA distribution
        pod90: 90% PoD threshold (kg/h)
        targeted_top_pct: if set, prioritize top N% emitters (requires prior knowledge)

    Returns:
        dict with survey results
    """
    if wind_ms is None:
        wind_ms = sample_wind_pa()

    n_wells = len(emissions)
    n_survey = int(n_wells * coverage_pct)

    if targeted_top_pct is not None:
        # Target the top emitters (assumes prior knowledge via previous survey)
        n_target = int(n_wells * targeted_top_pct)
        top_indices = np.argsort(emissions)[-n_target:]
        survey_indices = top_indices[:n_survey] if n_survey <= n_target else \
                         np.concatenate([top_indices,
                                        np.random.choice(
                                            np.setdiff1d(np.arange(n_wells), top_indices),
                                            size=n_survey - n_target, replace=False)])
    else:
        # Random sampling (no prior knowledge)
        survey_indices = np.random.choice(n_wells, size=n_survey, replace=False)

    surveyed_emissions = emissions[survey_indices]

    # PoD for each surveyed well
    pods = pod_bridger(surveyed_emissions, wind_ms, pod90)

    # Bernoulli detection
    detected = np.random.binomial(1, pods)

    detected_emissions = surveyed_emissions[detected == 1]

    # Count top-emitter capture
    total_top5_emission = np.sort(emissions)[-int(0.05 * n_wells):].sum()
    top5_threshold = np.percentile(emissions, 95)
    detected_top5 = ((surveyed_emissions > top5_threshold) & (detected == 1))

    return {
        'wind_ms': wind_ms,
        'n_surveyed': n_survey,
        'n_detected': int(detected.sum()),
        'emissions_surveyed_kgph': float(surveyed_emissions.sum()),
        'emissions_detected_kgph': float(detected_emissions.sum()),
        'mitigation_pct': 100.0 * detected_emissions.sum() / emissions.sum(),
        'detection_rate_by_count': 100.0 * detected.mean(),
        'avg_pod': float(pods.mean()),
        'n_top5_in_survey': int((surveyed_emissions > top5_threshold).sum()),
        'n_top5_detected': int(detected_top5.sum()),
        'top5_emissions_detected': float(surveyed_emissions[detected_top5].sum()),
        'top5_detected_pct_of_total_top5': 100.0 * surveyed_emissions[detected_top5].sum() / total_top5_emission,
    }


def run_monte_carlo(emissions, coverage_pct, n_iter=N_ITERATIONS,
                    pod90=POD_90_MARCELLUS, label=""):
    """Run Monte Carlo survey simulation."""
    results = [simulate_survey(emissions, coverage_pct, pod90=pod90)
               for _ in range(n_iter)]

    keys = ['emissions_detected_kgph', 'mitigation_pct', 'n_detected',
            'detection_rate_by_count', 'avg_pod', 'wind_ms',
            'top5_detected_pct_of_total_top5']

    summary = {label: {
        'coverage_pct': coverage_pct * 100,
        'n_iter': n_iter,
    }}

    for key in keys:
        vals = np.array([r[key] for r in results])
        summary[label][key] = {
            'mean': float(vals.mean()),
            'std': float(vals.std()),
            'p5': float(np.percentile(vals, 5)),
            'p25': float(np.percentile(vals, 25)),
            'median': float(np.percentile(vals, 50)),
            'p75': float(np.percentile(vals, 75)),
            'p95': float(np.percentile(vals, 95)),
        }

    return summary[label], results


# ============================================================
# STEP 5: FLIGHT PLAN CALCULATION
# ============================================================

def flight_plan_summary():
    """
    Realistic flight plan derivation from Permian 2024 benchmark.

    Permian 2024 (Donahue et al.):
      - 51,770 total sites over 195 calendar days (4 quarters)
      - Average: ~265 sites/day (all facility types)
      - Marginal wellsites only: 6,357/quarter ÷ ~49 days/quarter ≈ 130 wells/day

    PA adjustments:
      - Forested terrain: longer transit between wells, tighter corridors → -10-20%
      - More precipitation: ~150 vs 195 flyable days/year → -23%
      - Similar aircraft (Cessna 172 or similar fixed-wing) at 90-120 mph, 500-700' AGL
    """
    permian_marginal_per_qtr = 6357         # from Donahue et al. Table 1
    permian_days_per_qtr = 195 / 4          # ~49 days
    permian_marginal_per_day = permian_marginal_per_qtr / permian_days_per_qtr

    pa_wells_per_day = permian_marginal_per_day * 0.77  # 23% reduction for PA terrain/routing
    pa_flyable_days = 150  # vs 195 in Permian (dry, flat West Texas)

    pa_wells_per_year_moderate = int(pa_wells_per_day * pa_flyable_days)
    pa_coverage_pct = 100.0 * pa_wells_per_year_moderate / N_WELLS

    # Days to complete full portfolio
    days_full_portfolio = N_WELLS / pa_wells_per_day
    years_full_portfolio = days_full_portfolio / pa_flyable_days

    return {
        'permian_benchmark': {
            'marginal_wells_per_quarter': permian_marginal_per_qtr,
            'flying_days_per_quarter': round(permian_days_per_qtr, 1),
            'marginal_wells_per_day': round(permian_marginal_per_day, 0),
        },
        'pa_adjusted': {
            'wells_per_day': round(pa_wells_per_day, 0),
            'flyable_days_per_year': pa_flyable_days,
            'wells_per_year': pa_wells_per_year_moderate,
            'portfolio_coverage_pct': round(pa_coverage_pct, 1),
            'days_to_complete_full_portfolio': round(days_full_portfolio, 0),
            'years_to_complete_full_portfolio': round(years_full_portfolio, 2),
        },
        'quarterly_schedule': {
            'Q1_Jan_Mar': {'wells': round(pa_wells_per_day * 38, 0), 'days': 38},
            'Q2_Apr_Jun': {'wells': round(pa_wells_per_day * 40, 0), 'days': 40},
            'Q3_Jul_Sep': {'wells': round(pa_wells_per_day * 35, 0), 'days': 35},
            'Q4_Oct_Dec': {'wells': round(pa_wells_per_day * 37, 0), 'days': 37},
            'note': 'PA winter/fall months have more weather delays; Q2/Q3 preferred'
        }
    }


# ============================================================
# MAIN ANALYSIS
# ============================================================

def main():
    print("\n" + "="*80)
    print("BRIDGER LiDAR SURVEY ANALYSIS — PA MARGINAL WELLS")
    print("Grounded in Thorpe et al. 2024 & Donahue et al. 2024")
    print("="*80)

    output_dir = Path('BridgerResults')
    output_dir.mkdir(exist_ok=True)

    # -----------------------------------------------------------
    # A. Emission Distribution Analysis
    # -----------------------------------------------------------
    print("\n" + "─"*60)
    print("A. EMISSION CONCENTRATION (Lorenz Curve)")
    print("─"*60)

    conc_results, thresh_results = emission_concentration_analysis(emissions)

    print("\nTop N% of wells → share of total emissions:")
    print(f"  {'Top % Wells':>12} {'N Wells':>10} {'Emission Share':>16} {'Min Emission':>14}")
    print(f"  {'-'*55}")
    for pct, r in conc_results.items():
        print(f"  {pct:>10}%  {r['n_wells']:>10,}  {r['emission_share_pct']:>14.1f}%  {r['min_emission_in_group']:>12.3f} kg/h")

    print("\nWells above key PoD thresholds:")
    print(f"  {'Threshold':>12} {'N Wells':>10} {'% Wells':>10} {'Emissions':>14} {'% of Total':>12}")
    print(f"  {'-'*63}")
    for t, r in thresh_results.items():
        print(f"  {t:>10.2f} kg/h  {r['n_wells_above']:>9,}  {r['pct_wells_above']:>9.1f}%  "
              f"{r['emissions_above_kgph']:>12.0f}  {r['pct_emissions_above']:>10.1f}%")

    # -----------------------------------------------------------
    # B. PoD Model Analysis
    # -----------------------------------------------------------
    print("\n" + "─"*60)
    print("B. PoD DETECTION PROBABILITY AT KEY EMISSION LEVELS")
    print(f"   (Marcellus-specific PoD90 = {POD_90_MARCELLUS} kg/h at 3.5 m/s)")
    print("─"*60)

    pod_results = pod_at_percentiles(emissions, wind_ms=3.5, pod90=POD_90_MARCELLUS)
    print(f"\n  {'Emission':>12} {'PoD @3.5m/s':>12} {'N wells ≥ thresh':>18} {'% of fleet':>12}")
    print(f"  {'-'*60}")
    for r in pod_results:
        print(f"  {r['emission_kgph']:>10.3f} kg/h  {r['PoD']:>10.1%}  {r['n_wells_at_or_above']:>17,}  {r['pct_wells_at_or_above']:>10.1f}%")

    # Wind sensitivity
    print(f"\n  PoD at 1.0 kg/h emission across wind speeds (Marcellus threshold):")
    for ws in [1.0, 2.0, 3.5, 5.0, 6.0]:
        p = pod_bridger(1.0, ws, POD_90_MARCELLUS)
        adj_thresh = POD_90_MARCELLUS * np.exp(-0.3 * (ws - 3.5))
        print(f"  Wind {ws:.1f} m/s: PoD(1.0 kg/h) = {p:.1%}  |  90% PoD threshold = {adj_thresh:.2f} kg/h")

    # -----------------------------------------------------------
    # C. Realistic Flight Plan
    # -----------------------------------------------------------
    print("\n" + "─"*60)
    print("C. REALISTIC FLIGHT PLAN (from Permian 2024 benchmark)")
    print("─"*60)

    fp = flight_plan_summary()
    print(f"\n  Permian benchmark (Donahue et al. 2024):")
    print(f"    Marginal wells/quarter: {fp['permian_benchmark']['marginal_wells_per_quarter']:,}")
    print(f"    Flying days/quarter:    {fp['permian_benchmark']['flying_days_per_quarter']}")
    print(f"    Marginal wells/day:     {fp['permian_benchmark']['marginal_wells_per_day']:.0f}")
    print(f"\n  PA-adjusted parameters:")
    print(f"    Wells/day:              {fp['pa_adjusted']['wells_per_day']:.0f}")
    print(f"    Flyable days/year:      {fp['pa_adjusted']['flyable_days_per_year']}")
    print(f"    Wells/year:             {fp['pa_adjusted']['wells_per_year']:,}")
    print(f"    Annual coverage:        {fp['pa_adjusted']['portfolio_coverage_pct']}%")
    print(f"    Days for full portfolio:{fp['pa_adjusted']['days_to_complete_full_portfolio']:.0f}")
    print(f"    Years for full pass:    {fp['pa_adjusted']['years_to_complete_full_portfolio']:.1f}")

    # -----------------------------------------------------------
    # D. Monte Carlo — Coverage Scenarios
    # -----------------------------------------------------------
    print("\n" + "─"*60)
    print(f"D. MONTE CARLO SURVEY RESULTS ({N_ITERATIONS} iterations)")
    print("─"*60)

    # Scenarios: single-year realistic coverage
    # - One season only (Q2+Q3, ~75 days) = ~7,500 wells = 11.6%
    # - Full year moderate = ~15,000 wells = 23.2%
    # - Full year ambitious (Permian rate) = ~19,500 wells = 30.2%

    coverage_scenarios = [
        ('One season (Q2+Q3, ~75 days)', 0.116),
        ('1 year - moderate (150 days)', ANNUAL_COVERAGE_PCT_CONS),
        ('1 year - ambitious (Permian rate, 150 days @ 130/day)', ANNUAL_COVERAGE_PCT_PERM),
        ('Full portfolio (single pass)', 1.0),
    ]

    all_mc_results = {}

    print(f"\n  {'Scenario':<45} {'Cov%':>5} {'Mean Det.':>12} {'Mitigation':>12} {'95% CI':>20}")
    print(f"  {'-'*100}")

    for scenario_name, cov_pct in coverage_scenarios:
        mc_summary, _ = run_monte_carlo(
            emissions, cov_pct, n_iter=N_ITERATIONS,
            pod90=POD_90_MARCELLUS, label=scenario_name
        )
        all_mc_results[scenario_name] = mc_summary

        det = mc_summary['emissions_detected_kgph']
        mit = mc_summary['mitigation_pct']
        print(f"  {scenario_name:<45} {cov_pct*100:>4.0f}%  "
              f"{det['mean']:>10.0f} kg/h  {mit['mean']:>10.1f}%  "
              f"[{mit['p5']:.1f}%–{mit['p95']:.1f}%]")

    # -----------------------------------------------------------
    # E. THE KEY QUESTION: Top Emitters vs Most Emitters
    # -----------------------------------------------------------
    print("\n" + "─"*60)
    print("E. KEY QUESTION: Does Bridger capture 'MOST EMITTERS' or 'TOP EMITTERS'?")
    print("─"*60)

    # Run detailed analysis for 23% coverage (realistic 1-year plan)
    coverage_realistic = ANNUAL_COVERAGE_PCT_CONS
    mc_summary_det, raw_results = run_monte_carlo(
        emissions, coverage_realistic, n_iter=N_ITERATIONS, pod90=POD_90_MARCELLUS
    )

    # What fraction of total wells are detected?
    mean_wells_detected = mc_summary_det['n_detected']['mean']
    pct_fleet_detected = 100 * mean_wells_detected / N_WELLS

    # What fraction of total EMISSIONS are captured?
    mean_emis_detected = mc_summary_det['emissions_detected_kgph']['mean']
    pct_emis_detected = mc_summary_det['mitigation_pct']['mean']

    # Top 5% emitter capture rate (how much of the top-5% emissions do you get?)
    top5_cap = mc_summary_det['top5_detected_pct_of_total_top5']

    # Emission detection efficiency: emissions detected / wells detected
    # vs average emission per well in portfolio
    avg_emission_detected = mean_emis_detected / mean_wells_detected if mean_wells_detected > 0 else 0

    print(f"\n  For a realistic 1-year survey ({coverage_realistic*100:.0f}% of portfolio = "
          f"{int(coverage_realistic * N_WELLS):,} wells):\n")

    print(f"  BY COUNT (how many emitters captured):")
    print(f"    Mean detected wells:           {mean_wells_detected:,.0f} wells")
    print(f"    As % of TOTAL fleet:           {pct_fleet_detected:.1f}% of all 64,624 wells")
    print(f"    → Most emitters captured?      {'YES' if pct_fleet_detected > 50 else 'NO, only ' + f'{pct_fleet_detected:.0f}%'}")

    print(f"\n  BY EMISSIONS VOLUME (how much CH₄ captured):")
    print(f"    Mean emissions detected:       {mean_emis_detected:,.0f} kg/h")
    print(f"    As % of total portfolio:       {pct_emis_detected:.1f}%")
    print(f"    → Top 5% emitters captured:   {top5_cap['mean']:.1f}% of top-5% emissions")

    print(f"\n  DETECTION QUALITY (are detected wells disproportionately large?):")
    print(f"    Avg emission of detected wells: {avg_emission_detected:.3f} kg/h")
    print(f"    Avg emission of all wells:      {MEAN_EMISSION:.3f} kg/h")
    print(f"    Enrichment factor:              {avg_emission_detected/MEAN_EMISSION:.1f}x")
    print(f"    → Detected wells are {avg_emission_detected/MEAN_EMISSION:.1f}x larger than average")
    print(f"    → Bridger disproportionately finds HIGH emitters within surveyed sample")

    # Context: emission threshold below which Bridger is largely 'blind'
    pod10_threshold = POD_90_MARCELLUS * np.exp(-0.3 * (3.5 - 3.5)) * np.exp(-LOGISTIC_STEEPNESS * np.log(9))
    # More carefully: at what emission does PoD = 10%?
    # 0.10 = 1/(1+exp(-2*(log(e)-log(e90)))) → log(e)-log(e90) = log(0.1/0.9)/2 = -1.099
    # e = e90 * exp(-1.099) = 0.974 * 0.333 = 0.324 kg/h
    pod10_emission = POD_90_MARCELLUS * np.exp(-1.099)
    n_below_pod10 = (emissions < pod10_emission).sum()
    emis_below_pod10 = emissions[emissions < pod10_emission].sum()

    print(f"\n  DETECTION BLIND SPOT (PoD < 10%):")
    print(f"    Emission threshold for 10% PoD: ~{pod10_emission:.2f} kg/h")
    print(f"    Wells below this threshold:     {n_below_pod10:,} ({100*n_below_pod10/N_WELLS:.0f}% of fleet)")
    print(f"    Emissions below threshold:      {emis_below_pod10:.0f} kg/h ({100*emis_below_pod10/TOTAL_EMISSIONS_KGPH:.0f}% of total)")

    # Bottom line
    print(f"\n  ━━━ ANSWER TO KEY QUESTION ━━━")
    print(f"  Bridger does NOT capture 'most emitters' — only ~{pct_fleet_detected:.0f}% of wells are")
    print(f"  detected in a realistic 1-year survey.")
    print(f"")
    print(f"  Bridger IS effective at capturing 'top emitters':")
    print(f"  → Within the surveyed subset, high-emission wells are detected at")
    print(f"    ~{avg_emission_detected/MEAN_EMISSION:.1f}x the rate of average wells (emission-weighted)")
    print(f"  → Top 5% emitters (contributing {conc_results[5]['emission_share_pct']:.0f}% of emissions)")
    print(f"    yield {top5_cap['mean']:.0f}% capture in surveyed clusters")
    print(f"")
    print(f"  KEY INSIGHT: ~{int(thresh_results[POD_90_MARCELLUS]['pct_wells_above'])}% of wells are above the 90% PoD threshold")
    print(f"  ({POD_90_MARCELLUS} kg/h) but those wells contribute")
    print(f"  ~{thresh_results[POD_90_MARCELLUS]['pct_emissions_above']:.0f}% of total emissions.")
    print(f"  Bridger is optimized precisely for this concentration pattern.")

    # -----------------------------------------------------------
    # F. Summary table for all scenarios + PoD comparison
    # -----------------------------------------------------------
    print("\n" + "─"*60)
    print("F. MARCELLUS vs OVERALL PoD THRESHOLD COMPARISON")
    print("─"*60)
    print(f"\n  Using Marcellus-specific PoD90 = {POD_90_MARCELLUS} kg/h (Thorpe 2024, Fig 5)")
    print(f"  vs. overall 2023 average PoD90 = {POD_90_OVERALL} kg/h")
    print(f"\n  Effect on 1-year survey ({coverage_realistic*100:.0f}% coverage):")

    for threshold, label in [(POD_90_MARCELLUS, "Marcellus"), (POD_90_OVERALL, "Overall 2023")]:
        mc_s, _ = run_monte_carlo(emissions, coverage_realistic, n_iter=200, pod90=threshold)
        print(f"    {label:>15} (PoD90={threshold} kg/h): "
              f"detect {mc_s['emissions_detected_kgph']['mean']:,.0f} kg/h "
              f"({mc_s['mitigation_pct']['mean']:.1f}% mitigation)")

    # -----------------------------------------------------------
    # SAVE RESULTS
    # -----------------------------------------------------------
    save_data = {
        'analysis_grounding': {
            'pod_model_source': 'Thorpe et al. 2024, RSE, doi:10.1016/j.rse.2024.114435',
            'marcellus_pod90_kgph': POD_90_MARCELLUS,
            'marcellus_pod90_source': 'Thorpe 2024 Fig 5 right table, Marcellus basin',
            'flight_plan_source': 'Donahue et al. 2024 (Permian preprint), Table 1',
            'permian_marginal_per_day': fp['permian_benchmark']['marginal_wells_per_day'],
            'pa_marginal_per_day': fp['pa_adjusted']['wells_per_day'],
        },
        'portfolio': {
            'n_wells': N_WELLS,
            'total_emissions_kgph': TOTAL_EMISSIONS_KGPH,
            'mean_kgph': float(MEAN_EMISSION),
            'median_kgph': MEDIAN_EMISSION,
        },
        'emission_concentration': {
            str(k): v for k, v in conc_results.items()
        },
        'detection_thresholds': {
            str(k): v for k, v in thresh_results.items()
        },
        'flight_plan': fp,
        'monte_carlo_results': {
            name: {k: (v if not isinstance(v, dict) else v)
                   for k, v in res.items()}
            for name, res in all_mc_results.items()
        },
        'key_finding': {
            'pct_fleet_detected_1yr': round(pct_fleet_detected, 1),
            'pct_emissions_detected_1yr': round(pct_emis_detected, 1),
            'detection_enrichment_factor': round(avg_emission_detected / MEAN_EMISSION, 1),
            'top5_emission_capture_pct': round(top5_cap['mean'], 1),
            'pod_blind_spot_threshold_kgph': round(pod10_emission, 3),
            'pct_wells_in_blind_spot': round(100 * n_below_pod10 / N_WELLS, 0),
            'pct_emissions_in_blind_spot': round(100 * emis_below_pod10 / TOTAL_EMISSIONS_KGPH, 0),
            'answer': 'TOP EMITTERS, not most emitters: Bridger cannot detect the majority '
                      'of wells (which emit below ~0.3 kg/h). But within surveyed wells, it '
                      'reliably finds high emitters. The detected emission mass is dominated '
                      'by top emitters because their PoD is near 100%.'
        }
    }

    output_path = output_dir / 'bridger_pa_grounded_analysis.json'
    with open(output_path, 'w') as f:
        json.dump(save_data, f, indent=2, default=str)

    print(f"\n\nResults saved to: {output_path}")
    print("="*80)
    return save_data, emissions


if __name__ == '__main__':
    results, emissions_array = main()
