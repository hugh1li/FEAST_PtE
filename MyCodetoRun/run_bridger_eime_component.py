"""
run_bridger_eime_component.py
=============================
Bridger LiDAR survey simulation for PA marginal wells — Path A:
"EIME-anchored component decomposition."

Key idea
--------
We keep the EIME 250-iteration per-well site total as the authoritative
snapshot "truth" (`ch4_kgh_0..249` in the GeoPackage), then use FEAST's
component model to disaggregate that site total into individual leaks
at detection time:

  1. Pick a random EIME iteration k  →  each well i has total T_{i,k}.
  2. For each surveyed well:
       a. n_leaks ~ Binomial(N_COMPS_PER_SITE, EMISSION_PER_COMP),
          forced to >=1 when T > 0.
       b. Bootstrap n_leaks flux samples from FEAST's production_emissions.p.
       c. **Rescale** the vector so sum == T_{i,k}  (preserves EIME total).
  3. Optionally inject plunger / no-plunger unloading events at the survey
     instant via Bernoulli(p_active) with FEAST parameters from
     Example-RunScript.py. Episodic emissions are additive on top of the
     EIME fugitive total when --add-plungers is passed.
  4. Apply Bridger GML 2.0 POD (P4 predictor + Burr link) to every individual
     leak with its local wind; Bernoulli detect.  Compare against the
     site-total POD method used by bridger_pa_full_analysis.py.

Why this exists
---------------
bridger_pa_full_analysis.py applies POD to the site total directly —
implicitly treating the well as one big leak.  run_bridger_component_
analysis.py uses FEAST's bottom-up machinery but throws away the 250 EIME
iterations.  This script is the honest middle path: EIME tells us *how
much*, FEAST tells us the *component shape*, and we apply POD where the
physics actually happens — at the individual leak.

Outputs
-------
  BridgerResults/bridger_eime_component.csv         per-MC-iter rows
  BridgerResults/bridger_eime_component.json        summary stats
  BridgerResults/bridger_eime_component_compare.png diagnostic plots
"""

from __future__ import annotations

import argparse
import json
import pickle
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd

warnings.filterwarnings("ignore")

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# ── Paths ────────────────────────────────────────────────────────────────────
DATA_ITER = Path("/tmp/pa_data/PAMarginalWellEmissions2023_withEachIteration.gpkg")
EMISSION_DIST_PATH = REPO_ROOT / "ExampleData" / "DataObjectInstances" / "production_emissions.p"
OUT_DIR = Path(__file__).parent / "BridgerResults"
OUT_DIR.mkdir(exist_ok=True)

# ── Bridger GML 2.0 POD (Combined, P4 + Burr; Thorpe et al. 2024 Table 3) ─────
POD_ALPHA_1 = 2.0000
POD_ALPHA_2 = 1.5000
POD_BETA_1  = 2.41e-3
POD_BETA_2  = 1.9505
POD_BETA_3  = 2.0836
POD_BETA_4  = 1.5185
# GCN for Marcellus/PA terrain (ppm-m / 1000)
N_MARCELLUS = 16.0 / 1000.0

# ── Wind envelope (Bridger flyable) ──────────────────────────────────────────
MIN_WIND_MS = 1.0
MAX_WIND_MS = 6.0
WIND_MEAN_PA = 3.5
WIND_STD_PA  = 1.0

# ── FEAST component model (matches Example-RunScript.py fugitive values) ─────
N_COMPS_PER_SITE  = 650
EMISSION_PER_COMP = 0.00231
GS_TO_KGPH        = 3.6

# ── FEAST episodic (plunger / no plunger) from Example-RunScript.py ──────────
PLUNGER_SIZE_GS         = 50.36
PLUNGER_DURATION_DAYS   = 2059.0 / 86400.0
PLUNGER_RATE_PER_DAY    = 0.0165

NOPLUNGER_SIZE_GS       = 76.98
NOPLUNGER_DURATION_DAYS = 4973.5 / 86400.0
NOPLUNGER_RATE_PER_DAY  = 0.0002111

# Default fraction of PA marginal wells with a plunger lift (configurable).
# Public literature estimates range from ~40% to ~70% for low-pressure
# marginal verticals; we use 0.5 as a documented midpoint.  Override with
# --plunger-fraction.
DEFAULT_PLUNGER_FRACTION = 0.50

# ── Simulation controls ──────────────────────────────────────────────────────
DEFAULT_N_MC       = 500
DEFAULT_COVERAGE   = 0.23     # ~1-yr moderate campaign (15,000 wells/yr / 64,624)
GRID_DEG           = 0.1      # spatial survey cell size (~11 × 8.5 km)


# ═════════════════════════════════════════════════════════════════════════════
# POD
# ═════════════════════════════════════════════════════════════════════════════
def pod_gml2(e_kgph: np.ndarray, wind_ms: np.ndarray, n: float = N_MARCELLUS) -> np.ndarray:
    """Bridger GML 2.0 PoD (Thorpe et al. 2024, Combined Table 3)."""
    Q = np.maximum(np.asarray(e_kgph, dtype=float), 1e-9)
    u = np.clip(np.asarray(wind_ms, dtype=float), MIN_WIND_MS, MAX_WIND_MS)
    g = POD_BETA_1 * Q**POD_BETA_2 / (n**POD_BETA_3 * u**POD_BETA_4)
    return np.clip(1.0 - (1.0 + g**POD_ALPHA_1) ** (-POD_ALPHA_2), 0.0, 1.0)


# ═════════════════════════════════════════════════════════════════════════════
# Component decomposition (EIME-anchored)
# ═════════════════════════════════════════════════════════════════════════════
def decompose_sites_to_components(
    site_totals_kgph: np.ndarray,
    leak_pool_kgph: np.ndarray,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Vectorized: for each surveyed well with EIME total T_i, draw n_i leaks
    from Binomial(N_COMPS_PER_SITE, EMISSION_PER_COMP), bootstrap their raw
    flux vector from FEAST's pool, then **proportionally rescale** each
    well's leak vector so that sum(leaks_i) == T_i.

    Wells with T_i == 0 contribute no leaks.

    Returns
    -------
    leak_kgph   : (N_total_leaks,) component fluxes after rescaling, kg/h
    leak_site   : (N_total_leaks,) int array mapping each leak to its well
                  index within `site_totals_kgph`
    n_leaks_per : (len(site_totals_kgph),) leaks per well (0 for T==0)
    """
    n_wells = len(site_totals_kgph)
    emitting = site_totals_kgph > 0.0

    n_draws = rng.binomial(N_COMPS_PER_SITE, EMISSION_PER_COMP, size=n_wells)
    # force at least one leak for wells with T > 0
    n_draws = np.where(emitting & (n_draws < 1), 1, n_draws)
    # ensure zero-emission wells contribute zero leaks
    n_draws = np.where(emitting, n_draws, 0)

    total = int(n_draws.sum())
    if total == 0:
        return np.empty(0), np.empty(0, dtype=int), n_draws

    raw = rng.choice(leak_pool_kgph, size=total)
    leak_site = np.repeat(np.arange(n_wells), n_draws)

    # per-site raw sum and target → rescale factor
    raw_sum = np.zeros(n_wells)
    np.add.at(raw_sum, leak_site, raw)
    # avoid divide-by-zero (can't happen for emitting wells with n>=1 and
    # non-degenerate pool, but guard anyway)
    rescale = np.where(raw_sum > 0, site_totals_kgph / np.where(raw_sum > 0, raw_sum, 1.0), 0.0)
    leak_kgph = raw * rescale[leak_site]
    return leak_kgph, leak_site, n_draws


# ═════════════════════════════════════════════════════════════════════════════
# Episodic (plunger / no-plunger) events
# ═════════════════════════════════════════════════════════════════════════════
def sample_episodic_events(
    n_wells: int,
    rng: np.random.Generator,
    plunger_fraction: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    For each surveyed well, decide (i) whether it's a plunger well and
    (ii) whether an unloading event is currently active at the survey instant.

    P(active at instant) = rate_per_day * duration_days  — the steady-state
    occupancy of the continuous-time alternating-renewal process.

    Returns
    -------
    episodic_kgph  : (n_wells,) active episodic emission at the survey instant
                     (0 for wells with no active event)
    is_plunger     : (n_wells,) bool — well-type assignment
    """
    is_plunger = rng.random(n_wells) < plunger_fraction

    p_pl = PLUNGER_RATE_PER_DAY    * PLUNGER_DURATION_DAYS
    p_np = NOPLUNGER_RATE_PER_DAY  * NOPLUNGER_DURATION_DAYS

    active = np.zeros(n_wells, dtype=bool)
    r = rng.random(n_wells)
    active[is_plunger]  = r[is_plunger]  < p_pl
    active[~is_plunger] = r[~is_plunger] < p_np

    episodic_kgph = np.zeros(n_wells)
    episodic_kgph[is_plunger  & active] = PLUNGER_SIZE_GS   * GS_TO_KGPH
    episodic_kgph[~is_plunger & active] = NOPLUNGER_SIZE_GS * GS_TO_KGPH
    return episodic_kgph, is_plunger


# ═════════════════════════════════════════════════════════════════════════════
# Spatial wind field (re-used idea from bridger_pa_full_analysis.py)
# ═════════════════════════════════════════════════════════════════════════════
def sample_winds(n: int, rng: np.random.Generator) -> np.ndarray:
    base = float(np.clip(rng.normal(WIND_MEAN_PA, WIND_STD_PA), MIN_WIND_MS, MAX_WIND_MS))
    local = base + rng.normal(0, 0.35, size=n)
    return np.clip(local, MIN_WIND_MS, MAX_WIND_MS)


# ═════════════════════════════════════════════════════════════════════════════
# Contiguous spatial cell selection (same method as bridger_pa_full_analysis)
# ═════════════════════════════════════════════════════════════════════════════
def build_cell_map(lats: np.ndarray, lons: np.ndarray):
    lat_edges = np.arange(39.6, 42.4 + GRID_DEG, GRID_DEG)
    lon_edges = np.arange(-80.6, -75.8 + GRID_DEG, GRID_DEG)
    lat_idx = np.clip(np.digitize(lats, lat_edges) - 1, 0, len(lat_edges) - 2)
    lon_idx = np.clip(np.digitize(lons, lon_edges) - 1, 0, len(lon_edges) - 2)
    cell_id = lat_idx * 1000 + lon_idx
    unique_cells = np.unique(cell_id)
    cell_to_wells = {c: np.where(cell_id == c)[0] for c in unique_cells}
    cell_id_set = set(unique_cells.tolist())

    def neighbours(cid: int):
        lat_i, lon_i = divmod(cid, 1000)
        out = []
        for dlat, dlon in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nb = (lat_i + dlat) * 1000 + (lon_i + dlon)
            if nb in cell_id_set:
                out.append(nb)
        return out

    return unique_cells, cell_to_wells, neighbours


def contiguous_selection(
    unique_cells: np.ndarray,
    cell_to_wells: dict,
    neighbours_fn,
    target_wells: int,
    rng: np.random.Generator,
) -> np.ndarray:
    n_cells = len(unique_cells)
    seed = int(unique_cells[rng.integers(0, n_cells)])
    selected: list[int] = []
    visited = {seed}
    frontier = [seed]
    count = 0
    while frontier and count < target_wells:
        idx = int(rng.integers(0, len(frontier)))
        cid = frontier.pop(idx)
        w = cell_to_wells[cid]
        selected.extend(w.tolist())
        count += len(w)
        for nb in neighbours_fn(cid):
            if nb not in visited:
                visited.add(nb)
                frontier.append(nb)
        # ferry if exhausted
        if not frontier and count < target_wells:
            for c in unique_cells:
                c = int(c)
                if c not in visited:
                    visited.add(c)
                    frontier.append(c)
                    break
    return np.array(selected[:target_wells], dtype=int)


# ═════════════════════════════════════════════════════════════════════════════
# One MC iteration
# ═════════════════════════════════════════════════════════════════════════════
def run_one_iter(
    E_ITERS: np.ndarray,
    unique_cells, cell_to_wells, neighbours_fn,
    leak_pool_kgph: np.ndarray,
    target_wells: int,
    rng: np.random.Generator,
    add_plungers: bool,
    plunger_fraction: float,
) -> dict:
    n_wells_total, n_iter_eime = E_ITERS.shape

    # Pick EIME iteration
    k = int(rng.integers(0, n_iter_eime))

    # Pick contiguous survey cells
    sel = contiguous_selection(unique_cells, cell_to_wells, neighbours_fn, target_wells, rng)
    T = E_ITERS[sel, k].astype(np.float64)

    # Decompose to components
    leak_kgph, leak_site, n_leaks_per = decompose_sites_to_components(T, leak_pool_kgph, rng)

    # Winds per well
    winds = sample_winds(len(sel), rng)
    leak_winds = winds[leak_site] if len(leak_kgph) > 0 else np.empty(0)

    # Component-level detection
    if len(leak_kgph) > 0:
        p_comp = pod_gml2(leak_kgph, leak_winds)
        det_comp = rng.binomial(1, p_comp).astype(bool)
        detected_component = float(leak_kgph[det_comp].sum())
    else:
        detected_component = 0.0

    # Site-total-method detection (for comparison)
    p_site = pod_gml2(T, winds)
    det_site = rng.binomial(1, p_site).astype(bool)
    detected_site_total = float(T[det_site].sum())

    # Episodic (plunger) add-on
    if add_plungers:
        epi_kgph, is_pl = sample_episodic_events(len(sel), rng, plunger_fraction)
        # apply POD to each active episodic event individually
        active_mask = epi_kgph > 0
        p_epi = np.zeros_like(epi_kgph)
        if active_mask.any():
            p_epi[active_mask] = pod_gml2(epi_kgph[active_mask], winds[active_mask])
        det_epi = rng.binomial(1, p_epi).astype(bool)
        epi_total_surveyed = float(epi_kgph.sum())
        epi_detected       = float(epi_kgph[det_epi].sum())
        n_active_plunger   = int(((epi_kgph > 0) & is_pl).sum())
        n_active_noplunger = int(((epi_kgph > 0) & ~is_pl).sum())
    else:
        epi_total_surveyed = epi_detected = 0.0
        n_active_plunger = n_active_noplunger = 0

    surveyed_total_fugitive = float(T.sum())
    surveyed_total_all      = surveyed_total_fugitive + epi_total_surveyed
    detected_all_component  = detected_component + epi_detected
    detected_all_sitetotal  = detected_site_total + epi_detected  # same episodic branch

    safe_all = max(surveyed_total_all, 1e-9)
    safe_fug = max(surveyed_total_fugitive, 1e-9)

    return {
        "eime_iter": k,
        "n_surveyed": len(sel),
        "n_leaks": int(len(leak_kgph)),
        "leaks_per_well_mean": float(n_leaks_per[T > 0].mean()) if (T > 0).any() else 0.0,
        "wind_mean_ms": float(winds.mean()),
        "wind_min_ms":  float(winds.min()),
        "wind_max_ms":  float(winds.max()),
        "surveyed_fugitive_kgph":   surveyed_total_fugitive,
        "surveyed_episodic_kgph":   epi_total_surveyed,
        "surveyed_total_kgph":      surveyed_total_all,
        "detected_component_kgph":  detected_component,
        "detected_sitetotal_kgph":  detected_site_total,
        "detected_episodic_kgph":   epi_detected,
        "detected_all_component":   detected_all_component,
        "detected_all_sitetotal":   detected_all_sitetotal,
        # within-survey mitigation fractions
        "mit_component_pct":         100.0 * detected_component / safe_fug,
        "mit_sitetotal_pct":         100.0 * detected_site_total / safe_fug,
        "mit_component_allpct":      100.0 * detected_all_component / safe_all,
        "mit_sitetotal_allpct":      100.0 * detected_all_sitetotal / safe_all,
        "overestimate_pp":           100.0 * (detected_site_total - detected_component) / safe_fug,
        "n_active_plunger":          n_active_plunger,
        "n_active_noplunger":        n_active_noplunger,
    }


# ═════════════════════════════════════════════════════════════════════════════
# Main
# ═════════════════════════════════════════════════════════════════════════════
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-mc", type=int, default=DEFAULT_N_MC)
    ap.add_argument("--coverage", type=float, default=DEFAULT_COVERAGE,
                    help="fraction of wells surveyed per MC run (default 0.23 = 1-yr campaign)")
    ap.add_argument("--add-plungers", action="store_true",
                    help="add plunger / no-plunger unloading episodics on top of EIME fugitive")
    ap.add_argument("--plunger-fraction", type=float, default=DEFAULT_PLUNGER_FRACTION,
                    help=f"fraction of wells that are plunger-lifted (default {DEFAULT_PLUNGER_FRACTION})")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--data-iter", type=str, default=str(DATA_ITER),
                    help="path to PAMarginalWellEmissions2023_withEachIteration.gpkg")
    args = ap.parse_args()

    data_iter_path = Path(args.data_iter)
    if not data_iter_path.exists():
        raise FileNotFoundError(
            f"EIME iteration data not found at {data_iter_path}.\n"
            "Unzip PAMarginalWellEmissions2023.zip into /tmp/pa_data first, or "
            "override with --data-iter."
        )
    if not EMISSION_DIST_PATH.exists():
        raise FileNotFoundError(f"Missing FEAST leak distribution: {EMISSION_DIST_PATH}")

    print("=" * 78)
    print("BRIDGER × EIME COMPONENT DECOMPOSITION  (Path A)")
    print("=" * 78)
    print(f"  EIME data:         {data_iter_path}")
    print(f"  Leak distribution: {EMISSION_DIST_PATH.relative_to(REPO_ROOT)}")
    print(f"  MC iterations:     {args.n_mc}")
    print(f"  Coverage:          {args.coverage*100:.1f}%")
    print(f"  Episodic events:   {'ON' if args.add_plungers else 'off'}"
          f"{'  (plunger_fraction=' + str(args.plunger_fraction) + ')' if args.add_plungers else ''}")

    # ── Load EIME iterations ────────────────────────────────────────────────
    print("\nLoading EIME GeoPackage (large, one-time) …")
    gdf = gpd.read_file(data_iter_path)
    lats = gdf.Latitude.values
    lons = gdf.Longitude.values
    iter_cols = sorted([c for c in gdf.columns if c.startswith("ch4_kgh_")])
    E_ITERS = gdf[iter_cols].values.astype(np.float32)
    n_wells, n_iter_eime = E_ITERS.shape
    print(f"  {n_wells:,} wells × {n_iter_eime} EIME iterations "
          f"({E_ITERS.nbytes/1e6:.0f} MB array)")

    # ── Load FEAST leak distribution ────────────────────────────────────────
    with open(EMISSION_DIST_PATH, "rb") as f:
        leak_obj = pickle.load(f)
    leak_pool_gs = np.concatenate([np.asarray(v, dtype=float) for v in leak_obj.leak_sizes.values()])
    leak_pool_kgph = leak_pool_gs * GS_TO_KGPH
    leak_pool_kgph = leak_pool_kgph[leak_pool_kgph > 0]  # drop zeros
    print(f"  FEAST leak pool: n={len(leak_pool_kgph):,}  "
          f"mean={leak_pool_kgph.mean():.4f} kg/h  "
          f"median={np.median(leak_pool_kgph):.4f}  "
          f"max={leak_pool_kgph.max():.2f}")

    # ── Spatial grid ────────────────────────────────────────────────────────
    unique_cells, cell_to_wells, neighbours_fn = build_cell_map(lats, lons)
    print(f"  Survey grid: {len(unique_cells)} occupied {GRID_DEG}° cells")

    target_wells = int(n_wells * args.coverage)

    # ── Monte Carlo ─────────────────────────────────────────────────────────
    print(f"\nRunning {args.n_mc} MC iterations (contiguous cells, random EIME snapshot)…")
    rng = np.random.default_rng(args.seed)
    rows = []
    for i in range(args.n_mc):
        rows.append(run_one_iter(
            E_ITERS, unique_cells, cell_to_wells, neighbours_fn,
            leak_pool_kgph, target_wells, rng,
            add_plungers=args.add_plungers,
            plunger_fraction=args.plunger_fraction,
        ))
        if (i + 1) % max(1, args.n_mc // 10) == 0:
            print(f"  iter {i+1}/{args.n_mc}", flush=True)
    res = pd.DataFrame(rows)

    # ── Summary ─────────────────────────────────────────────────────────────
    def stats(s: pd.Series) -> dict:
        return dict(
            mean=float(s.mean()), std=float(s.std()),
            p5=float(s.quantile(0.05)), p50=float(s.quantile(0.50)),
            p95=float(s.quantile(0.95)),
        )

    summary = {
        "config": {
            "n_mc": args.n_mc, "coverage": args.coverage,
            "target_wells": target_wells,
            "n_comps_per_site": N_COMPS_PER_SITE,
            "emission_per_comp": EMISSION_PER_COMP,
            "add_plungers": args.add_plungers,
            "plunger_fraction": args.plunger_fraction,
            "seed": args.seed,
            "n_wells_total": int(n_wells),
            "n_iter_eime": int(n_iter_eime),
        },
        "within_survey_mitigation_pct_fugitive": {
            "component_method":  stats(res["mit_component_pct"]),
            "sitetotal_method":  stats(res["mit_sitetotal_pct"]),
            "overestimate_pp":   stats(res["overestimate_pp"]),
        },
        "within_survey_mitigation_pct_all": {
            "component_method":  stats(res["mit_component_allpct"]),
            "sitetotal_method":  stats(res["mit_sitetotal_allpct"]),
        } if args.add_plungers else None,
        "detected_kgph": {
            "component_method": stats(res["detected_component_kgph"]),
            "sitetotal_method": stats(res["detected_sitetotal_kgph"]),
        },
        "surveyed_kgph": {
            "fugitive": stats(res["surveyed_fugitive_kgph"]),
            "episodic": stats(res["surveyed_episodic_kgph"]) if args.add_plungers else None,
        },
        "wind_ms": {
            "mean": stats(res["wind_mean_ms"]),
        },
        "leaks_per_well_mean": stats(res["leaks_per_well_mean"]),
    }

    suffix = "_withplungers" if args.add_plungers else ""
    out_csv  = OUT_DIR / f"bridger_eime_component{suffix}.csv"
    out_json = OUT_DIR / f"bridger_eime_component{suffix}.json"
    res.to_csv(out_csv, index=False)
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print("\n" + "─" * 70)
    print("RESULTS (within-survey mitigation %, EIME fugitive denominator)")
    print("─" * 70)
    m_c = summary["within_survey_mitigation_pct_fugitive"]["component_method"]
    m_s = summary["within_survey_mitigation_pct_fugitive"]["sitetotal_method"]
    m_o = summary["within_survey_mitigation_pct_fugitive"]["overestimate_pp"]
    print(f"  Component-level POD:  {m_c['mean']:5.1f}% ± {m_c['std']:.1f}   90% CI [{m_c['p5']:.0f}–{m_c['p95']:.0f}]")
    print(f"  Site-total POD:       {m_s['mean']:5.1f}% ± {m_s['std']:.1f}   90% CI [{m_s['p5']:.0f}–{m_s['p95']:.0f}]")
    print(f"  Site-total − comp:   +{m_o['mean']:5.1f} pp   (bias from treating a well as one big leak)")
    if args.add_plungers:
        m_all_c = summary["within_survey_mitigation_pct_all"]["component_method"]
        m_all_s = summary["within_survey_mitigation_pct_all"]["sitetotal_method"]
        print()
        print(f"  With episodics (EIME fugitive + unloadings) — component: {m_all_c['mean']:.1f}%, site-total: {m_all_s['mean']:.1f}%")
        epi = summary["surveyed_kgph"]["episodic"]
        print(f"  Surveyed episodic contribution: {epi['mean']:.0f} ± {epi['std']:.0f} kg/h")

    print(f"\n  Mean leaks/emitting-well: {summary['leaks_per_well_mean']['mean']:.2f}  "
          f"(FEAST expected ≈ {N_COMPS_PER_SITE*EMISSION_PER_COMP:.2f})")
    print(f"\n  Saved  {out_csv.relative_to(REPO_ROOT)}")
    print(f"  Saved  {out_json.relative_to(REPO_ROOT)}")

    # ── Diagnostic plot ─────────────────────────────────────────────────────
    if not args.no_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
        fig.suptitle(
            f"EIME × Component Decomposition (n_mc={args.n_mc}, coverage={args.coverage*100:.0f}%, "
            f"plungers={'on' if args.add_plungers else 'off'})",
            fontweight="bold",
        )

        ax = axes[0]
        ax.hist(res["mit_component_pct"], bins=40, alpha=0.7, label="Component POD", color="#1976D2")
        ax.hist(res["mit_sitetotal_pct"], bins=40, alpha=0.55, label="Site-total POD", color="#D32F2F")
        ax.axvline(res["mit_component_pct"].mean(), color="#1976D2", lw=2)
        ax.axvline(res["mit_sitetotal_pct"].mean(), color="#D32F2F", lw=2, ls="--")
        ax.set_xlabel("Within-survey mitigation %")
        ax.set_ylabel("Count")
        ax.set_title("(A) Mitigation % — component vs site-total")
        ax.legend()

        ax = axes[1]
        ax.scatter(res["mit_sitetotal_pct"], res["mit_component_pct"], s=6, alpha=0.5)
        lims = [0, max(res["mit_sitetotal_pct"].max(), res["mit_component_pct"].max()) * 1.05]
        ax.plot(lims, lims, "k--", alpha=0.5, label="y = x")
        ax.set_xlabel("Site-total POD mitigation %")
        ax.set_ylabel("Component POD mitigation %")
        ax.set_title("(B) Per-iteration parity")
        ax.legend()

        ax = axes[2]
        ax.hist(res["leaks_per_well_mean"], bins=30, color="#388E3C", alpha=0.8)
        ax.axvline(N_COMPS_PER_SITE * EMISSION_PER_COMP, color="red", ls="--",
                   label=f"FEAST expected = {N_COMPS_PER_SITE*EMISSION_PER_COMP:.2f}")
        ax.set_xlabel("Mean leaks per emitting well (per iter)")
        ax.set_ylabel("Count")
        ax.set_title("(C) Leak frequency realized")
        ax.legend()

        plt.tight_layout()
        out_png = OUT_DIR / f"bridger_eime_component_compare{suffix}.png"
        plt.savefig(out_png, dpi=140, bbox_inches="tight")
        print(f"  Saved  {out_png.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
