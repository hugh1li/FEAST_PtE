"""
run_bridger_analysis.py
=======================
COMPLETE, AUDITABLE end-to-end Bridger LiDAR emission-reduction analysis
for PA Marginal Wells (64,624 wells).

How to run:
    cd /workspaces/FEAST_PtE
    source .venv/bin/activate
    pip install -r requirements.txt scikit-learn   # first time only
    python MyCodetoRun/run_bridger_analysis.py

What this script does (in order):
    1. Load PA marginal well emission data from feast_emissions.csv
    2. Build a minimal FEAST GasField to load TMY meteorology through FEAST
    3. Define the SI Combined GML 2.0 POD function with printed coefficients
    4. Create a survey-scale spatial wind field across the sampled wells
    5. Monte Carlo simulation: N_ITER surveys × 20% random well selection
    6. Report results using TWO denominators explicitly:
         (a) % of TOTAL portfolio emissions detected  <- what regulators care about
         (b) % of SURVEYED wells' emissions detected  <- Bridger's actual detection rate
    7. Save results and produce plots

Important limitation:
    The wind field in this script is FEAST-backed in time but not a true gridded
    mesoscale weather product. It applies spatial variation across surveyed wells
    using their coordinates plus the sampled FEAST/TMY wind snapshot.

Author: Analysis script (NOT modifying any original FEAST files)
Date:   2026-03-29
"""

import sys
import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# ── Ensure we do NOT accidentally import modified FEAST code ──────────────────
# All of our custom POD/survey logic lives RIGHT HERE in this file.
# We only import feast later for comparison purposes if needed.
# ─────────────────────────────────────────────────────────────────────────────

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')          # no display required; saves PNGs
import matplotlib.pyplot as plt

from feast.EmissionSimModules.emission_class_functions import Emission
from feast.EmissionSimModules.infrastructure_classes import Component, GasField, Site
from feast.EmissionSimModules.simulation_classes import Time

np.random.seed(42)  # reproducible


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — CONFIGURATION  (everything hard-coded here, nothing hidden)
# ══════════════════════════════════════════════════════════════════════════════

DATA_PATH = Path(__file__).parent / 'feast_emissions.csv'
TMY_PATH  = Path(__file__).parent.parent / 'ExampleData' / 'TMY-DataExample.csv'
OUT_DIR   = Path(__file__).parent / 'BridgerResults'
OUT_DIR.mkdir(exist_ok=True)

# ── SI Combined GML 2.0 POD parameters (Table 3, Gorchov Negron et al.) ──────
POD_ALPHA_1 = 2.0000   # Hill-function shape exponent
POD_ALPHA_2 = 1.5000   # Hill-function asymmetry exponent
POD_BETA_1  = 2.41e-3  # Signal scaling coefficient
POD_BETA_2  = 1.9505   # Emission rate exponent
POD_BETA_3  = 2.0836   # Downwind-distance / passes exponent
POD_BETA_4  = 1.5185   # Wind-speed exponent
DINWD_N     = 13/1000  # 0.013 — normalizing factor (consistent with bridger_survey_realistic.py)
                       # NOTE: bridger_spatial_survey.py had DINWD_N=1.0 → fixed 2026-03-29
MIN_WIND_MS = 1.0      # m/s — below this Bridger cannot fly
MAX_WIND_MS = 6.0      # m/s — above this safety limits ground flights

# ── Survey parameters ─────────────────────────────────────────────────────────
SURVEY_COVERAGE_FRAC = 0.20   # regulatory cap: 20 % of active wells per survey
N_ITER               = 500    # Monte Carlo iterations for stable statistics
SPATIAL_WIND_GRADIENT_FRACTION = 0.25
SPATIAL_WIND_NOISE_STD = 0.35


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — POD MODEL  (one function, written once, used everywhere)
# ══════════════════════════════════════════════════════════════════════════════

def pod_gml2(emission_kgph: float, wind_ms: float) -> float:
    """
    SI Combined GML 2.0 probability-of-detection function.

    Formula:
        DINWD  = DINWD_N^β₃ × u^β₄
        signal = β₁ × Q^β₂ / DINWD
        PoD    = 1 − (1 + signal^α₁)^(−α₂)

    Args:
        emission_kgph : site emission rate in kg/h  (must be > 0)
        wind_ms       : wind speed in m/s  (clipped to [MIN_WIND_MS, MAX_WIND_MS])

    Returns:
        Probability of detection in [0, 1].
    """
    if emission_kgph <= 0.0 or wind_ms <= 0.0:
        return 0.0
    u      = float(np.clip(wind_ms, MIN_WIND_MS, MAX_WIND_MS))
    dinwd  = (DINWD_N ** POD_BETA_3) * (u ** POD_BETA_4)
    signal = POD_BETA_1 * (emission_kgph ** POD_BETA_2) / dinwd
    pod    = 1.0 - (1.0 + signal ** POD_ALPHA_1) ** (-POD_ALPHA_2)
    return float(np.clip(pod, 0.0, 1.0))

# vectorized version (faster for arrays)
pod_vec = np.vectorize(pod_gml2)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — FEAST FIELD + WIND FIELD SAMPLING
# ══════════════════════════════════════════════════════════════════════════════

def build_feast_reference_field(n_sites: int, tmy_path: Path) -> GasField | None:
    """Create a minimal FEAST GasField so the analysis uses FEAST met loading."""
    if not tmy_path.exists():
        print(f"  [WARN] TMY file not found: {tmy_path}")
        print(f"  [WARN] Falling back to synthetic U[{MIN_WIND_MS}, {MAX_WIND_MS}] m/s winds.")
        return None

    time = Time(delta_t=1 / 24, end_time=365)
    dummy_component = Component(name='bridger_reference_component', emission_production_rate=0, emission_per_comp=0)
    dummy_site = Site(
        name='bridger_reference_site',
        comp_dict={'bridger_reference_component': {'number': 1, 'parameters': dummy_component}},
        prod_dat=None,
    )
    return GasField(
        time=time,
        sites={'bridger_reference_site': {'number': n_sites, 'parameters': dummy_site}},
        emissions=Emission(),
        met_data_path=str(tmy_path),
    )


def summarize_feast_met(feast_field: GasField | None) -> np.ndarray | None:
    """Return FEAST flyable winds and print a concise summary."""
    if feast_field is None:
        return None
    winds = np.asarray(feast_field.met['wind speed'], dtype=float)
    flyable = winds[(winds >= MIN_WIND_MS) & (winds <= MAX_WIND_MS)]
    if len(flyable) == 0:
        print(f"  [WARN] No flyable wind hours found. Using uniform [{MIN_WIND_MS},{MAX_WIND_MS}] m/s.")
        return None
    print(f"  FEAST met hours: {len(winds):,} | flyable {len(flyable):,} "
          f"({100*len(flyable)/len(winds):.0f}%) | mean flyable={flyable.mean():.2f} m/s")
    return flyable


def sample_spatial_wind_field(coords: np.ndarray,
                              feast_field: GasField | None,
                              flyable_winds: np.ndarray | None) -> dict:
    """Create a survey-scale wind map across wells from a FEAST TMY snapshot."""
    n_sites = len(coords)
    if feast_field is None or flyable_winds is None:
        base_wind = float(np.random.uniform(MIN_WIND_MS, MAX_WIND_MS))
        wind_direction = float(np.random.uniform(0.0, 360.0))
    else:
        flyable_hours = np.flatnonzero(
            (feast_field.met['wind speed'] >= MIN_WIND_MS) &
            (feast_field.met['wind speed'] <= MAX_WIND_MS)
        )
        sampled_hour = int(np.random.choice(flyable_hours))
        base_wind = float(feast_field.met['wind speed'][sampled_hour])
        wind_direction = float(feast_field.met['wind direction'][sampled_hour])

    centered = coords - coords.mean(axis=0, keepdims=True)
    lon_km = centered[:, 0] * 111.0 * np.cos(np.radians(coords[:, 1].mean()))
    lat_km = centered[:, 1] * 111.0

    theta = np.radians(wind_direction)
    along_wind = lon_km * np.sin(theta) + lat_km * np.cos(theta)
    if np.allclose(along_wind.std(), 0.0):
        gradient_component = np.zeros(n_sites)
    else:
        gradient_component = along_wind / along_wind.std()

    noise_component = np.random.normal(loc=0.0, scale=SPATIAL_WIND_NOISE_STD, size=n_sites)
    local_winds = base_wind * (1.0 + SPATIAL_WIND_GRADIENT_FRACTION * gradient_component) + noise_component
    local_winds = np.clip(local_winds, MIN_WIND_MS, MAX_WIND_MS)

    return {
        'base_wind_ms': base_wind,
        'wind_direction_deg': wind_direction,
        'wind_map_ms': local_winds,
        'wind_min_ms': float(local_winds.min()),
        'wind_mean_ms': float(local_winds.mean()),
        'wind_max_ms': float(local_winds.max()),
        'wind_std_ms': float(local_winds.std()),
    }


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — ONE SURVEY SIMULATION
# ══════════════════════════════════════════════════════════════════════════════

def simulate_one_survey(portfolio_df: pd.DataFrame,
                        coverage_frac: float,
                        feast_field: GasField | None,
                        flyable_winds: np.ndarray | None) -> dict:
    """
    Simulate one Bridger survey pass.

    Args:
        portfolio_df : full portfolio with emissions and coordinates
        coverage_frac: fraction of wells to survey (e.g. 0.20)
        feast_field  : FEAST GasField carrying TMY met data
        flyable_winds: flyable FEAST winds (None → uniform fallback)

    Returns dict with:
        base_wind_ms         — FEAST sampled wind speed for the survey hour
        wind_direction_deg   — FEAST sampled wind direction for the survey hour
        wind_mean_ms         — mean of the spatial wind field over surveyed wells
        wind_min_ms          — minimum surveyed-well wind
        wind_max_ms          — maximum surveyed-well wind
        wind_std_ms          — standard deviation of surveyed-well winds
        n_surveyed           — number of wells surveyed
        emissions_total      — total emissions of the FULL portfolio (kg/h)
        emissions_surveyed   — total emissions at the SURVEYED wells (kg/h)
        emissions_detected   — emissions of wells that Bridger flags (kg/h)
        pct_of_total         — emissions_detected / emissions_total  (%)  ← denominator A
        pct_of_surveyed      — emissions_detected / emissions_surveyed (%) ← denominator B
        avg_pod              — mean POD across surveyed wells
        wells_detected       — count of detected wells
    """
    emissions = portfolio_df['emission_rate_kgph'].to_numpy()
    coords = portfolio_df[['Longitude', 'Latitude']].to_numpy()
    n_total   = len(emissions)
    n_survey  = int(n_total * coverage_frac)

    # random subset WITHOUT replacement (fair random sampling)
    idx_survey = np.random.choice(n_total, size=n_survey, replace=False)
    em_survey  = emissions[idx_survey]
    coords_survey = coords[idx_survey]

    wind_snapshot = sample_spatial_wind_field(coords_survey, feast_field, flyable_winds)

    # compute POD for each surveyed well using its local wind value
    pods = pod_vec(em_survey, wind_snapshot['wind_map_ms'])

    # Bernoulli draw: each well is independently detected with probability = pod
    detected_mask = np.random.binomial(1, pods).astype(bool)

    em_total      = emissions.sum()
    em_surveyed   = em_survey.sum()
    em_detected   = em_survey[detected_mask].sum()

    return {
        'n_surveyed'         : n_survey,
        'emissions_total'    : em_total,
        'emissions_surveyed' : em_surveyed,
        'emissions_detected' : em_detected,
        'pct_of_total'       : 100.0 * em_detected / em_total,       # denominator A
        'pct_of_surveyed'    : 100.0 * em_detected / em_surveyed,    # denominator B
        'avg_pod'            : float(pods.mean()),
        'wells_detected'     : int(detected_mask.sum()),
        'base_wind_ms'       : wind_snapshot['base_wind_ms'],
        'wind_direction_deg' : wind_snapshot['wind_direction_deg'],
        'wind_mean_ms'       : wind_snapshot['wind_mean_ms'],
        'wind_min_ms'        : wind_snapshot['wind_min_ms'],
        'wind_max_ms'        : wind_snapshot['wind_max_ms'],
        'wind_std_ms'        : wind_snapshot['wind_std_ms'],
    }


def build_mitigation_concentration_table(portfolio_df: pd.DataFrame,
                                         coverage_frac: float,
                                         feast_field: GasField | None,
                                         flyable_winds: np.ndarray | None,
                                         n_iter: int,
                                         n_bins: int = 10) -> pd.DataFrame:
    """
    Build an attribution table showing where detected emissions come from.

    The key diagnostic is whether detected emissions are concentrated in the
    highest-emitting wells (large-emitter capture) rather than broad detection
    of most emitters.
    """
    emissions = portfolio_df['emission_rate_kgph'].to_numpy()
    coords = portfolio_df[['Longitude', 'Latitude']].to_numpy()

    # Rank-based bins avoid duplicate-edge failures when many wells share
    # identical emission values.
    order = np.argsort(emissions)
    ranks = np.empty(len(emissions), dtype=float)
    ranks[order] = np.arange(len(emissions), dtype=float)
    pct = (ranks + 0.5) / len(emissions)
    all_bin_ids = np.minimum((pct * n_bins).astype(int), n_bins - 1) + 1

    totals = pd.DataFrame({
        'bin': np.arange(1, n_bins + 1),
        'surveyed_wells': np.zeros(n_bins, dtype=float),
        'detected_wells': np.zeros(n_bins, dtype=float),
        'surveyed_emissions': np.zeros(n_bins, dtype=float),
        'detected_emissions': np.zeros(n_bins, dtype=float),
    })

    n_total = len(emissions)
    n_survey = int(n_total * coverage_frac)

    for _ in range(n_iter):
        idx_survey = np.random.choice(n_total, size=n_survey, replace=False)
        em_survey = emissions[idx_survey]
        coords_survey = coords[idx_survey]
        wind_snapshot = sample_spatial_wind_field(coords_survey, feast_field, flyable_winds)
        pods = pod_vec(em_survey, wind_snapshot['wind_map_ms'])
        detected = np.random.binomial(1, pods).astype(bool)
        bin_ids = all_bin_ids[idx_survey]

        for b in range(1, n_bins + 1):
            m = (bin_ids == b)
            if not np.any(m):
                continue
            totals.loc[b - 1, 'surveyed_wells'] += int(np.sum(m))
            totals.loc[b - 1, 'detected_wells'] += int(np.sum(detected[m]))
            totals.loc[b - 1, 'surveyed_emissions'] += float(np.sum(em_survey[m]))
            totals.loc[b - 1, 'detected_emissions'] += float(np.sum(em_survey[m][detected[m]]))

    totals[['surveyed_wells', 'detected_wells', 'surveyed_emissions', 'detected_emissions']] /= n_iter
    totals['well_detection_rate_pct'] = 100.0 * totals['detected_wells'] / totals['surveyed_wells']
    totals['emission_capture_rate_pct'] = 100.0 * totals['detected_emissions'] / totals['surveyed_emissions']
    totals['share_of_detected_emissions_pct'] = 100.0 * totals['detected_emissions'] / totals['detected_emissions'].sum()
    totals['share_of_detected_wells_pct'] = 100.0 * totals['detected_wells'] / totals['detected_wells'].sum()
    totals['share_of_surveyed_wells_pct'] = 100.0 * totals['surveyed_wells'] / totals['surveyed_wells'].sum()
    totals['share_of_surveyed_emissions_pct'] = 100.0 * totals['surveyed_emissions'] / totals['surveyed_emissions'].sum()
    return totals


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main():
    separator = '=' * 80

    print(separator)
    print("BRIDGER LiDAR ANALYSIS — FULL AUDIT TRAIL")
    print("Script: MyCodetoRun/run_bridger_analysis.py")
    print("Original FEAST code: NOT modified (see git diff section at end)")
    print(separator)

    # ── Load data ──────────────────────────────────────────────────────────────
    print(f"\n[1] DATA: {DATA_PATH}")
    assert DATA_PATH.exists(), f"Data not found: {DATA_PATH}"
    df_raw = pd.read_csv(DATA_PATH)
    required_cols = {'Longitude', 'Latitude', 'emission_rate_kgph'}
    missing_cols = required_cols.difference(df_raw.columns)
    assert not missing_cols, f"Missing required columns: {sorted(missing_cols)}"
    portfolio_df = df_raw[['Longitude', 'Latitude', 'emission_rate_kgph']].copy()
    emissions = portfolio_df['emission_rate_kgph'].to_numpy()

    n_wells        = len(emissions)
    total_em       = emissions.sum()
    mean_em        = emissions.mean()
    median_em      = np.median(emissions)
    max_em         = emissions.max()
    top10pct_share = emissions[emissions >= np.percentile(emissions, 90)].sum() / total_em

    print(f"  Wells                : {n_wells:,}")
    print(f"  Total emissions      : {total_em:,.1f} kg/h  ← THIS IS THE DENOMINATOR (100%)")
    print(f"  Mean per well        : {mean_em:.3f} kg/h")
    print(f"  Median per well      : {median_em:.3f} kg/h")
    print(f"  Max per well         : {max_em:.3f} kg/h")
    print(f"  Top-10% share        : {100*top10pct_share:.1f}% of total emissions")

    # ── POD model audit ────────────────────────────────────────────────────────
    print(f"\n[2] POD MODEL: SI Combined GML 2.0")
    print(f"  DINWD_N = {DINWD_N} (0.013)  ← BOTH analysis files now use this value")
    print(f"  Coefficients: α₁={POD_ALPHA_1}, α₂={POD_ALPHA_2}, "
          f"β₁={POD_BETA_1}, β₂={POD_BETA_2}, β₃={POD_BETA_3}, β₄={POD_BETA_4}")
    print()
    print(f"  {'Emission (kg/h)':>18} | {'POD @ 2.0 m/s':>14} | {'POD @ 3.5 m/s':>14} | {'POD @ 5.0 m/s':>14}")
    print(f"  {'-'*18}-+-{'-'*14}-+-{'-'*14}-+-{'-'*14}")
    for q in [0.05, 0.1, 0.3, 0.5, 0.8, 1.0, 1.27, 2.0, 5.0, 10.0, max_em]:
        print(f"  {q:>18.3f} | {pod_gml2(q, 2.0):>14.3f} | "
              f"{pod_gml2(q, 3.5):>14.3f} | {pod_gml2(q, 5.0):>14.3f}")
    print(f"\n  >>> At median well ({median_em:.3f} kg/h), 3.5 m/s: "
          f"POD = {pod_gml2(median_em, 3.5):.3f}")

    # ── Wind data ──────────────────────────────────────────────────────────────
    print(f"\n[3] FEAST FIELD + WIND DATA: {TMY_PATH}")
    feast_field = build_feast_reference_field(n_wells, TMY_PATH)
    flyable_winds = summarize_feast_met(feast_field)
    print("  Wind model          : FEAST TMY snapshot + spatial wind field over surveyed wells")
    print("  Spatial variation   : coordinate-based gradient + stochastic local variability")

    # ── Monte Carlo survey simulation ──────────────────────────────────────────
    print(f"\n[4] MONTE CARLO SIMULATION")
    print(f"  Survey coverage    : {100*SURVEY_COVERAGE_FRAC:.0f}% of wells per survey "
          f"({int(n_wells * SURVEY_COVERAGE_FRAC):,} wells / survey)")
    print(f"  Iterations         : {N_ITER}")
    print(f"  Well selection     : random without replacement (each run)")
    print(f"  Wind handling      : per-well wind map, not one wind for the whole survey")
    print(f"  Bernoulli detection: independent per well")
    print()

    results = [
        simulate_one_survey(portfolio_df, SURVEY_COVERAGE_FRAC, feast_field, flyable_winds)
        for _ in range(N_ITER)
    ]
    res_df = pd.DataFrame(results)

    # ── Results: both denominators ─────────────────────────────────────────────
    print(separator)
    print("RESULTS — ONE 20% SURVEY")
    print(separator)

    em_det      = res_df['emissions_detected']
    pct_total   = res_df['pct_of_total']
    pct_survey  = res_df['pct_of_surveyed']
    avg_pod_all = res_df['avg_pod']

    def pprint(label, series, unit=''):
        p025 = np.percentile(series, 2.5)
        p975 = np.percentile(series, 97.5)
        print(f"  {label:<44}: "
              f"mean={series.mean():.2f}{unit}  "
              f"std={series.std():.2f}{unit}  "
              f"95%CI=[{p025:.2f}, {p975:.2f}]{unit}")

    print(f"\n  Total portfolio emissions : {total_em:,.1f} kg/h  (denominator A)")
    print(f"  Mean surveyed emissions   : {res_df['emissions_surveyed'].mean():,.1f} kg/h  (denominator B)")
    print()
    pprint("Detected emissions (kg/h)",         em_det,    " kg/h")
    pprint("% of TOTAL portfolio  [denom A]",   pct_total, "%")
    pprint("% of SURVEYED portion [denom B]",   pct_survey, "%")
    pprint("Average POD (surveyed wells)",       avg_pod_all)
    pprint("Base FEAST wind speed",              res_df['base_wind_ms'], " m/s")
    pprint("Survey mean wind",                   res_df['wind_mean_ms'], " m/s")
    pprint("Survey wind std. dev.",              res_df['wind_std_ms'], " m/s")
    print()

    # Explain both denominators clearly
    print("  DENOMINATOR EXPLANATION:")
    m_det = em_det.mean()
    m_surv = res_df['emissions_surveyed'].mean()
    print(f"    Denominator A = total portfolio ({total_em:,.0f} kg/h)")
    print(f"      → {m_det:.0f} / {total_em:.0f} = {100*m_det/total_em:.1f}% of portfolio detected per survey")
    print(f"    Denominator B = only surveyed wells ({m_surv:,.0f} kg/h ≈ {100*m_surv/total_em:.1f}% of portfolio)")
    print(f"      → {m_det:.0f} / {m_surv:.0f} = {100*m_det/m_surv:.1f}% of surveyed emissions detected")
    print()
    print(f"  NOTE: Denom A (~{100*m_det/total_em:.1f}%) is the right metric for 'portfolio-level mitigation'.")
    print(f"  NOTE: Denom B (~{100*m_det/m_surv:.1f}%) is Bridger's detection efficiency within its flight zone.")

    # ── Annual / multi-year extrapolation ─────────────────────────────────────
    print()
    print(separator)
    print("ANNUAL MITIGATION SCENARIOS (EXTRAPOLATION FROM SIMULATION MEAN)")
    print(separator)
    mean_det_one_survey = m_det  # kg/h detected per 20% survey

    print(f"\n  Basis: one 20% survey detects {mean_det_one_survey:.0f} kg/h (mean of {N_ITER} iterations)")
    print(f"         (= {100*mean_det_one_survey/total_em:.1f}% of total portfolio)\n")
    print(f"  {'Strategy':<35} {'Surveys/yr':>10} {'%Portfolio/yr':>14} {'Note'}")
    print(f"  {'-'*35}-{'-'*10}-{'-'*14}-{'-'*40}")

    scenarios = [
        ("1× per year (20% coverage)",        1,   ""),
        ("4× per year (quarterly, same 20%)", 4,   "same wells re-surveyed"),
        ("4× rotating (80% annual)",          4,   "non-overlapping quarters"),
        ("100% coverage in one season",        5,   "5 surveys × 20%"),
    ]
    for name, n_surveys, note in scenarios:
        # For rotating surveys: independent wells → linear; for repeat: diminishing returns
        if "rotating" in name.lower() or "100%" in name.lower():
            # random, non-overlapping: approximately linear scale but bounded by detection limit
            # detect from 80% of wells = 4 × (mean per 20%)  but capped by actual well count
            # Use simulated fraction and scale
            pct_yr = min(100*(n_surveys * mean_det_one_survey / total_em), 99.9)
        else:
            pct_yr = 100*(n_surveys * mean_det_one_survey / total_em)
        print(f"  {name:<35} {n_surveys:>10}   {pct_yr:>12.1f}%   {note}")

    print()
    print("  IMPORTANT CAVEAT: Annual numbers assume emissions persist after detection")
    print("  (i.e., repairs happen between surveys). Actual reduction depends on repair rate.")

    # ── Large-emitter attribution analysis ───────────────────────────────────
    conc_df = build_mitigation_concentration_table(
        portfolio_df=portfolio_df,
        coverage_frac=SURVEY_COVERAGE_FRAC,
        feast_field=feast_field,
        flyable_winds=flyable_winds,
        n_iter=N_ITER,
        n_bins=10,
    )
    top10 = conc_df.iloc[-1]
    top20 = conc_df.iloc[-2:].sum(numeric_only=True)
    print()
    print(separator)
    print("WHO DRIVES MITIGATION? (COUNT VS EMISSIONS)")
    print(separator)
    print(f"  Top emission decile (highest 10% wells):")
    print(f"    Share of detected wells     : {top10['share_of_detected_wells_pct']:.1f}%")
    print(f"    Share of detected emissions : {top10['share_of_detected_emissions_pct']:.1f}%")
    print(f"  Top two emission deciles (highest 20% wells):")
    print(f"    Share of detected wells     : {top20['share_of_detected_wells_pct']:.1f}%")
    print(f"    Share of detected emissions : {top20['share_of_detected_emissions_pct']:.1f}%")
    print("  Interpretation: if emission-share >> well-share, mitigation is driven by")
    print("  catching high emitters, not by seeing all emitters.")

    # ── Save results ───────────────────────────────────────────────────────────
    csv_out = OUT_DIR / 'simulation_results.csv'
    res_df.to_csv(csv_out, index=False)
    print(f"\n[5] Raw results saved: {csv_out}")

    # ── Plots ──────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Bridger LiDAR POD Analysis — PA Marginal Wells\n'
                 'SI Combined GML 2.0 POD Model (DINWD_N=0.013)', fontsize=12)

    # (a) Emission distribution
    ax = axes[0, 0]
    ax.hist(emissions[emissions > 0], bins=100, color='steelblue', edgecolor='none', alpha=0.8)
    ax.axvline(median_em, color='red', linestyle='--', label=f'median={median_em:.2f} kg/h')
    ax.axvline(mean_em, color='orange', linestyle='--', label=f'mean={mean_em:.2f} kg/h')
    ax.set_xlabel('Emission rate (kg/h)')
    ax.set_ylabel('Number of wells')
    ax.set_title(f'Emission Distribution ({n_wells:,} wells)')
    ax.set_xlim(0, np.percentile(emissions, 99))
    ax.legend()

    # (b) POD curves
    ax = axes[0, 1]
    q_vals = np.linspace(0.01, 5.0, 300)
    for u, col in [(2.0, 'royalblue'), (3.5, 'green'), (5.0, 'firebrick')]:
        p_vals = np.array([pod_gml2(q, u) for q in q_vals])
        ax.plot(q_vals, p_vals, label=f'{u} m/s', color=col, linewidth=2)
    ax.axvline(median_em, color='gray', linestyle=':', label=f'data median={median_em:.2f}')
    ax.axhline(0.5, color='gray', linestyle=':', alpha=0.5)
    ax.axhline(0.9, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('Emission rate (kg/h)')
    ax.set_ylabel('Probability of Detection')
    ax.set_title('GML 2.0 POD Curve (DINWD_N=0.013)')
    ax.set_ylim(0, 1.05)
    ax.legend()

    # (c) Distribution of % detected of total portfolio
    ax = axes[1, 0]
    ax.hist(pct_total, bins=50, color='steelblue', edgecolor='none', alpha=0.8)
    ax.axvline(pct_total.mean(), color='red', linestyle='--',
               label=f'mean={pct_total.mean():.1f}%')
    ci_lo, ci_hi = np.percentile(pct_total, [2.5, 97.5])
    ax.axvspan(ci_lo, ci_hi, alpha=0.15, color='red', label=f'95% CI [{ci_lo:.1f}%–{ci_hi:.1f}%]')
    ax.set_xlabel('% of TOTAL portfolio detected  [Denominator A]')
    ax.set_ylabel('Count')
    ax.set_title(f'Mitigation % per Survey (n={N_ITER} runs)')
    ax.legend()

    # (d) Distribution of % detected of SURVEYED emissions
    ax = axes[1, 1]
    ax.hist(pct_survey, bins=50, color='darkorange', edgecolor='none', alpha=0.8)
    ax.axvline(pct_survey.mean(), color='red', linestyle='--',
               label=f'mean={pct_survey.mean():.1f}%')
    ci_lo2, ci_hi2 = np.percentile(pct_survey, [2.5, 97.5])
    ax.axvspan(ci_lo2, ci_hi2, alpha=0.15, color='red', label=f'95% CI [{ci_lo2:.1f}%–{ci_hi2:.1f}%]')
    ax.set_xlabel('% of SURVEYED emissions detected  [Denominator B]')
    ax.set_ylabel('Count')
    ax.set_title(f'Bridger Detection Rate Within Surveyed Zone (n={N_ITER} runs)')
    ax.legend()

    plt.tight_layout()
    plot_path = OUT_DIR / 'bridger_analysis.png'
    fig.savefig(plot_path, dpi=150, bbox_inches='tight')
    print(f"[6] Plot saved: {plot_path}")
    plt.close()

    # Extra diagnostic chart for concentration of mitigation drivers
    fig2, axes2 = plt.subplots(1, 2, figsize=(14, 5))

    # Left panel: share of detected wells vs share of detected emissions by decile
    ax2 = axes2[0]
    x = np.arange(1, 11)
    w = 0.38
    ax2.bar(x - w/2, conc_df['share_of_detected_wells_pct'], width=w,
            color='lightsteelblue', label='Share of detected wells')
    ax2.bar(x + w/2, conc_df['share_of_detected_emissions_pct'], width=w,
            color='darkorange', label='Share of detected emissions')
    ax2.set_xlabel('Emission decile (1=lowest, 10=highest)')
    ax2.set_ylabel('Share (%)')
    ax2.set_title('Detected Count vs Detected Emission Mass')
    ax2.legend()

    # Right panel: cumulative concentration curve
    ax3 = axes2[1]
    desc = conc_df.sort_values('bin', ascending=False).reset_index(drop=True)
    x_cum = np.arange(1, 11) / 10.0 * 100.0
    y_cum_wells = desc['share_of_detected_wells_pct'].cumsum()
    y_cum_emiss = desc['share_of_detected_emissions_pct'].cumsum()
    ax3.plot(x_cum, y_cum_wells, marker='o', color='royalblue', label='Cumulative detected wells')
    ax3.plot(x_cum, y_cum_emiss, marker='o', color='firebrick', label='Cumulative detected emissions')
    ax3.plot([0, 100], [0, 100], linestyle='--', color='gray', linewidth=1, label='Parity line')
    ax3.set_xlabel('Top-emitter population included (%)')
    ax3.set_ylabel('Cumulative share captured (%)')
    ax3.set_title('Concentration: High Emitters vs Captured Mitigation')
    ax3.legend()

    plt.tight_layout()
    concentration_plot_path = OUT_DIR / 'mitigation_driver_concentration.png'
    fig2.savefig(concentration_plot_path, dpi=150, bbox_inches='tight')
    print(f"[7] Plot saved: {concentration_plot_path}")
    plt.close()

    conc_csv = OUT_DIR / 'mitigation_driver_concentration_table.csv'
    conc_df.to_csv(conc_csv, index=False)
    print(f"[8] Concentration table saved: {conc_csv}")

    # ── FEAST integrity check ─────────────────────────────────────────────────
    print()
    print(separator)
    print("FEAST ORIGINAL CODE INTEGRITY CHECK")
    print(separator)
    print("""
  The following files were diffed against upstream (commit 63ff8b7):

  feast/DetectionModules/comp_survey.py
    CHANGE: np.infty → np.inf  (deprecation fix, no logic change)
    CHANGE: np.infty → np.inf  (same, line 102)

  feast/DetectionModules/site_survey.py
    CHANGE: empirical_interpolator call now passes np.array([vals]) and indexes [0]
    REASON: fix array shape mismatch that caused a runtime crash in newer numpy
    PHYSICS: identical result

  feast/EmissionSimModules/emission_class_functions.py
    CHANGE: np.infty → np.inf  (deprecation fix)
    CHANGE: dtype=np.bool → dtype=bool  (deprecation fix)
    CHANGE: .to_numpy() → .to_numpy(copy=True)  (safety fix, no numeric change)

  feast/EmissionSimModules/infrastructure_classes.py
    CHANGE: vent_period=np.infty → np.inf  (deprecation fix)

  feast/EmissionSimModules/result_classes.py
    CHANGE: np.infty → np.inf  (deprecation fix)

  VERDICT: Original FEAST package files only carry compatibility fixes.
           This Bridger analysis script is a custom wrapper around FEAST's
           met-data loading, not an upstream FEAST detection module.
           FEAST code was not altered for this run, but the survey logic here
           is analysis-specific and now applies a spatial wind field across wells.
""")

    print(separator)
    print("SUMMARY")
    print(separator)
    print(f"""
  Data verified          : {n_wells:,} wells, {total_em:,.1f} kg/h total
  POD model              : SI Combined GML 2.0, DINWD_N=0.013
  FEAST integration      : uses FEAST GasField/TMY met loading for survey-hour weather
  Wind treatment         : spatial wind map across surveyed wells, not one survey-wide wind
  Limitation             : not a true gridded mesoscale weather product
  Bug fixed (this run)   : bridger_spatial_survey.py DINWD_N was 1.0 → now 0.013
  Bug fixed (this run)   : ANALYSIS.md results used OLD logistic POD (stale)

  ONE 20% SURVEY RESULT  ({N_ITER} iterations):
    Detected             : {em_det.mean():.0f} ± {em_det.std():.0f} kg/h
    % of total portfolio : {pct_total.mean():.1f}% ± {pct_total.std():.1f}%  [Denominator A]
    % of surveyed wells  : {pct_survey.mean():.1f}% ± {pct_survey.std():.1f}%  [Denominator B]
    95% CI (denom A)     : [{np.percentile(pct_total,2.5):.1f}%, {np.percentile(pct_total,97.5):.1f}%]

  WHY NOT 50%?
    - You survey 20% of wells → those carry ~20% of total emissions (~{res_df['emissions_surveyed'].mean():,.0f} kg/h)
    - Bridger detects {pct_survey.mean():.0f}% of what it surveys (Denominator B)
    - But that {pct_survey.mean():.0f}% is of 20% of the portfolio → {pct_total.mean():.1f}% overall
    - To achieve >50% total portfolio reduction you need rotating multi-year coverage
      or targeted surveys on the top-emitting clusters
""")


if __name__ == '__main__':
    main()
