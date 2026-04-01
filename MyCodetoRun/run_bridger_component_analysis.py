"""
run_bridger_component_analysis.py
=================================
FEAST bottom-up Bridger analysis: generates component-level emissions from
FEAST's real emission distribution, then compares per-component POD vs
site-total POD.

Method (consistent with FEAST internals):
- For each site: draw n_leaks ~ Binomial(N_COMPS_PER_SITE, EMISSION_PER_COMP)
- Bootstrap individual leak sizes from production_emissions.p (g/s)
- Site total *emerges* from the sum of individual leaks
- Apply Bridger POD to each individual component leak (bottom-up method)
- Also apply Bridger POD to the site total (site-total method) for comparison

Previous approach (Dirichlet split) was the reverse: it took a fixed site
total and split it proportionally into components. That conserves the total
but is inconsistent with how FEAST generates emissions.
"""

import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


# Paths
DATA_PATH = Path(__file__).parent / "feast_emissions.csv"
TMY_PATH = Path(__file__).parent.parent / "ExampleData" / "TMY-DataExample.csv"
EMISSION_DIST_PATH = (
    Path(__file__).parent.parent
    / "ExampleData"
    / "DataObjectInstances"
    / "production_emissions.p"
)
OUT_DIR = Path(__file__).parent / "BridgerResults"
OUT_DIR.mkdir(exist_ok=True)


# Bridger SI Combined GML 2.0 POD parameters
POD_ALPHA_1 = 2.0000
POD_ALPHA_2 = 1.5000
POD_BETA_1 = 2.41e-3
POD_BETA_2 = 1.9505
POD_BETA_3 = 2.0836
POD_BETA_4 = 1.5185
DINWD_N = 13 / 1000
MIN_WIND_MS = 1.0
MAX_WIND_MS = 6.0


# Simulation settings
RNG_SEED = 42
N_ITER = 500
SURVEY_COVERAGE_FRAC = 0.20
SPATIAL_WIND_GRADIENT_FRACTION = 0.25
SPATIAL_WIND_NOISE_STD = 0.35


# FEAST emission model parameters (from Example-RunScript.py / production_emissions.p calibration)
# N_COMPS_PER_SITE=650, EMISSION_PER_COMP=0.00231 → ~1.5 leaks/site at mean 0.0807 g/s ≈ 0.43 kg/h
# which matches the site totals in feast_emissions.csv.
N_COMPS_PER_SITE = 650       # mechanical components per site
EMISSION_PER_COMP = 0.00231  # fraction of components leaking (binomial P)
GS_TO_KGPH = 3.6             # g/s → kg/h


def pod_gml2_array(emission_kgph: np.ndarray, wind_ms: np.ndarray) -> np.ndarray:
    """Vectorized Bridger SI Combined GML 2.0 POD."""
    emission = np.asarray(emission_kgph, dtype=float)
    wind = np.asarray(wind_ms, dtype=float)
    emission = np.where(emission > 0.0, emission, 0.0)
    wind = np.clip(wind, MIN_WIND_MS, MAX_WIND_MS)
    dinwd = (DINWD_N ** POD_BETA_3) * (wind ** POD_BETA_4)
    signal = POD_BETA_1 * (emission ** POD_BETA_2) / dinwd
    pod = 1.0 - (1.0 + signal ** POD_ALPHA_1) ** (-POD_ALPHA_2)
    return np.clip(pod, 0.0, 1.0)


def load_emission_distribution(path: Path) -> np.ndarray:
    """
    Load bootstrapped leak sizes (g/s) from a FEAST LeakData pickle.

    Pools all detection-method keys so we draw from the combined distribution.
    """
    with open(path, "rb") as f:
        params = pickle.load(f)
    all_sizes: list[float] = []
    for v in params.leak_sizes.values():
        all_sizes.extend(v)
    return np.array(all_sizes, dtype=float)


def load_tmy_flyable_winds(tmy_path: Path) -> np.ndarray | None:
    """
    Return array of wind speeds (m/s) from TMY hours that are within the
    Bridger operational window [MIN_WIND_MS, MAX_WIND_MS].

    FEAST reads TMY with header=1 (second CSV row contains column names).
    """
    if not tmy_path.exists():
        return None
    met = pd.read_csv(tmy_path, header=1)
    speeds = met["wind speed (m/s)"].to_numpy(dtype=float)
    flyable = speeds[(speeds >= MIN_WIND_MS) & (speeds <= MAX_WIND_MS)]
    return flyable if len(flyable) > 0 else None


def generate_site_component_emissions(
    n_sites: int,
    leak_dist: np.ndarray,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Fully vectorized FEAST bottom-up emission generation for one MC iteration.

    For each site:
      1. n_leaks ~ Binomial(N_COMPS_PER_SITE, EMISSION_PER_COMP)  (vectorized)
      2. Bootstrap all leaks at once from the FEAST distribution

    Returns:
      site_totals : (n_sites,) array of site-level emission totals (kg/h)
      leak_kgph   : flat array of all individual leak rates (kg/h)
      leak_site   : flat int array mapping each leak to its site index
    """
    n_leaks_per_site = rng.binomial(N_COMPS_PER_SITE, EMISSION_PER_COMP, size=n_sites)
    total_leaks = int(n_leaks_per_site.sum())

    if total_leaks > 0:
        leak_kgph = rng.choice(leak_dist, size=total_leaks) * GS_TO_KGPH
    else:
        leak_kgph = np.array([], dtype=float)

    leak_site = np.repeat(np.arange(n_sites), n_leaks_per_site)

    site_totals = np.zeros(n_sites, dtype=float)
    if total_leaks > 0:
        np.add.at(site_totals, leak_site, leak_kgph)

    return site_totals, leak_kgph, leak_site


def sample_spatial_wind_field(
    coords: np.ndarray,
    flyable_winds: np.ndarray | None,
    rng: np.random.Generator,
) -> dict:
    """Sample one survey-hour wind speed and build a spatial field over sites."""
    n_sites = len(coords)

    if flyable_winds is not None and len(flyable_winds) > 0:
        base_wind = float(rng.choice(flyable_winds))
    else:
        base_wind = float(rng.uniform(MIN_WIND_MS, MAX_WIND_MS))
    wind_direction = float(rng.uniform(0.0, 360.0))

    centered = coords - coords.mean(axis=0, keepdims=True)
    lon_km = centered[:, 0] * 111.0 * np.cos(np.radians(coords[:, 1].mean()))
    lat_km = centered[:, 1] * 111.0
    theta = np.radians(wind_direction)
    along_wind = lon_km * np.sin(theta) + lat_km * np.cos(theta)
    gradient = (
        np.zeros(n_sites) if np.allclose(along_wind.std(), 0.0)
        else along_wind / along_wind.std()
    )

    noise = rng.normal(loc=0.0, scale=SPATIAL_WIND_NOISE_STD, size=n_sites)
    local_winds = base_wind * (1.0 + SPATIAL_WIND_GRADIENT_FRACTION * gradient) + noise
    local_winds = np.clip(local_winds, MIN_WIND_MS, MAX_WIND_MS)
    return {
        "base_wind_ms": base_wind,
        "wind_direction_deg": wind_direction,
        "local_winds_ms": local_winds,
        "wind_mean_ms": float(local_winds.mean()),
        "wind_std_ms": float(local_winds.std()),
    }


def simulate_one_survey(
    coords: np.ndarray,
    leak_dist: np.ndarray,
    flyable_winds: np.ndarray | None,
    rng: np.random.Generator,
) -> dict:
    """
    One MC survey iteration using the FEAST bottom-up approach (fully vectorized).

    Steps:
      1. Vectorized bottom-up emission generation for ALL sites.
      2. Randomly select SURVEY_COVERAGE_FRAC of sites to survey.
      3a. Component method: apply Bridger POD to each individual leak (vectorized).
      3b. Site-total method: apply Bridger POD to the site total.

    The denominator (emissions_total) varies each iteration because each
    run generates a fresh stochastic emission landscape.
    """
    n_total = len(coords)
    n_survey = int(n_total * SURVEY_COVERAGE_FRAC)

    # Bottom-up: generate all site emissions this iteration (vectorized)
    all_totals, all_leaks, all_leak_site = generate_site_component_emissions(
        n_total, leak_dist, rng
    )
    em_total = float(all_totals.sum())

    # Survey a random 20% subset
    idx = rng.choice(n_total, size=n_survey, replace=False)
    surveyed_totals = all_totals[idx]
    em_surveyed = float(surveyed_totals.sum())

    wind_snap = sample_spatial_wind_field(coords[idx], flyable_winds, rng)
    local_winds = wind_snap["local_winds_ms"]  # shape (n_survey,)

    # --- Method B: site-total POD (simplified, for comparison) ---
    pod_site = pod_gml2_array(surveyed_totals, local_winds)
    detected_site = rng.binomial(1, pod_site).astype(bool)
    em_detected_site_total = float(surveyed_totals[detected_site].sum())

    # --- Method A: per-component POD (FEAST bottom-up, vectorized) ---
    # Build a boolean mask for leaks belonging to surveyed sites.
    # Then broadcast each site's wind to all its leaks via a lookup array.
    if len(all_leaks) > 0:
        is_surveyed = np.zeros(n_total, dtype=bool)
        is_surveyed[idx] = True
        comp_mask = is_surveyed[all_leak_site]      # True for leaks from surveyed sites
        comp_leaks = all_leaks[comp_mask]
        comp_site_global = all_leak_site[comp_mask]

        # Map global site index → local wind (0 for unsurveyed, but those are excluded)
        global_to_wind = np.zeros(n_total, dtype=float)
        global_to_wind[idx] = local_winds
        comp_winds = global_to_wind[comp_site_global]

        comp_pods = pod_gml2_array(comp_leaks, comp_winds)
        comp_detected = rng.binomial(1, comp_pods).astype(bool)
        em_detected_component = float(comp_leaks[comp_detected].sum())
    else:
        em_detected_component = 0.0

    safe_total = max(em_total, 1e-9)
    safe_comp = max(em_detected_component, 1e-9)
    return {
        "n_surveyed": n_survey,
        "emissions_total": em_total,
        "emissions_surveyed": em_surveyed,
        "emissions_detected_site_total": em_detected_site_total,
        "emissions_detected_component": em_detected_component,
        "pct_total_site_total": 100.0 * em_detected_site_total / safe_total,
        "pct_total_component": 100.0 * em_detected_component / safe_total,
        "pct_surveyed_site_total": 100.0 * em_detected_site_total / max(em_surveyed, 1e-9),
        "pct_surveyed_component": 100.0 * em_detected_component / max(em_surveyed, 1e-9),
        "overestimate_pct_points": 100.0 * (em_detected_site_total - em_detected_component) / safe_total,
        "overestimate_relative_pct": 100.0 * (em_detected_site_total - em_detected_component) / safe_comp,
        "base_wind_ms": wind_snap["base_wind_ms"],
        "wind_mean_ms": wind_snap["wind_mean_ms"],
        "wind_std_ms": wind_snap["wind_std_ms"],
    }

def main() -> None:
    print("=" * 80)
    print("BRIDGER BOTTOM-UP COMPONENT ANALYSIS (FEAST-CONSISTENT)")
    print("=" * 80)

    if not DATA_PATH.exists():
        raise FileNotFoundError(f"Missing data file: {DATA_PATH}")
    if not EMISSION_DIST_PATH.exists():
        raise FileNotFoundError(f"Missing emission distribution: {EMISSION_DIST_PATH}")

    df = pd.read_csv(DATA_PATH)
    coords = df[["Longitude", "Latitude"]].to_numpy()
    n_sites = len(coords)

    leak_dist = load_emission_distribution(EMISSION_DIST_PATH)
    flyable_winds = load_tmy_flyable_winds(TMY_PATH)

    expected_leaks = N_COMPS_PER_SITE * EMISSION_PER_COMP
    expected_total_kgph = expected_leaks * leak_dist.mean() * GS_TO_KGPH

    print(f"Sites: {n_sites:,}")
    print(f"Emission distribution: {len(leak_dist):,} bootstrap samples  "
          f"median={np.median(leak_dist):.4f} g/s  mean={leak_dist.mean():.4f} g/s")
    print(f"Components per site: {N_COMPS_PER_SITE}, emission_per_comp: {EMISSION_PER_COMP}")
    print(f"Expected leaks/site: {expected_leaks:.2f}  "
          f"expected site total: {expected_total_kgph:.3f} kg/h")
    if flyable_winds is not None:
        print(f"TMY flyable-window hours: {len(flyable_winds):,}  "
              f"mean {flyable_winds.mean():.2f} m/s")
    print(f"Iterations: {N_ITER}, survey coverage: {100 * SURVEY_COVERAGE_FRAC:.0f}%")
    print()

    rng = np.random.default_rng(RNG_SEED)
    rows = [simulate_one_survey(coords, leak_dist, flyable_winds, rng) for _ in range(N_ITER)]
    res = pd.DataFrame(rows)

    def ci(s: pd.Series) -> tuple[float, float]:
        return float(np.percentile(s, 2.5)), float(np.percentile(s, 97.5))

    lo_s, hi_s = ci(res["pct_total_site_total"])
    lo_c, hi_c = ci(res["pct_total_component"])
    lo_b, hi_b = ci(res["overestimate_pct_points"])

    print("RESULTS (FEAST BOTTOM-UP, ONE 20% SURVEY)")
    print("-" * 80)
    print(f"  Mean total emissions / iter:  {res['emissions_total'].mean():.1f} kg/h")
    print(f"  Mean surveyed emissions:       {res['emissions_surveyed'].mean():.1f} kg/h")
    print()
    print(
        f"  Site-total POD:   mean={res['pct_total_site_total'].mean():.2f}%  "
        f"95%CI=[{lo_s:.2f}, {hi_s:.2f}]"
    )
    print(
        f"  Component POD:    mean={res['pct_total_component'].mean():.2f}%  "
        f"95%CI=[{lo_c:.2f}, {hi_c:.2f}]"
    )
    print(
        f"  Overestimate:     mean={res['overestimate_pct_points'].mean():.2f} pp  "
        f"95%CI=[{lo_b:.2f}, {hi_b:.2f}]"
    )
    print(
        f"  Relative bias:    mean={res['overestimate_relative_pct'].mean():.1f}% "
        f"(site-total vs component-level)"
    )

    out_csv = OUT_DIR / "simulation_results_component_level.csv"
    out_summary = OUT_DIR / "component_vs_sitelevel_summary.json"
    res.to_csv(out_csv, index=False)

    summary = {
        "approach": "feast_bottom_up",
        "n_sites": int(n_sites),
        "n_comps_per_site": int(N_COMPS_PER_SITE),
        "emission_per_comp": float(EMISSION_PER_COMP),
        "emission_dist_n_samples": int(len(leak_dist)),
        "emission_dist_mean_gs": float(leak_dist.mean()),
        "emission_dist_median_gs": float(np.median(leak_dist)),
        "iterations": int(N_ITER),
        "survey_coverage_frac": float(SURVEY_COVERAGE_FRAC),
        "site_total_method": {
            "mean_pct_total": float(res["pct_total_site_total"].mean()),
            "ci95_pct_total": [lo_s, hi_s],
        },
        "component_method": {
            "mean_pct_total": float(res["pct_total_component"].mean()),
            "ci95_pct_total": [lo_c, hi_c],
        },
        "overestimate": {
            "mean_pct_points": float(res["overestimate_pct_points"].mean()),
            "ci95_pct_points": [lo_b, hi_b],
            "mean_relative_pct": float(res["overestimate_relative_pct"].mean()),
        },
    }
    with open(out_summary, "w") as f:
        json.dump(summary, f, indent=2)

    print()
    print(f"Saved: {out_csv}")
    print(f"Saved: {out_summary}")


if __name__ == "__main__":
    main()
