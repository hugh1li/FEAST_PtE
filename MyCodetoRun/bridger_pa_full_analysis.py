"""
Bridger LiDAR Survey Analysis — PA Marginal Wells (Full Real-Data Version)
==========================================================================
Uses actual GeoPackage data with 250 FEAST stochastic emission iterations.

Key methodological choices:
  1. Emissions: each survey simulation draws one random FEAST iteration
     (ch4_kgh_0...249) as the snapshot emission state at time of survey.
     `mean_ch4_kgh` is used only for portfolio characterization.

  2. Spatial survey: wells are selected as geographic grid cells (0.1° x 0.1°),
     mimicking how Bridger plans routes over contiguous areas, not random
     individual wells scattered across PA.

  3. Denominator: mitigation% = detected_emissions / SURVEYED_wells_total_emissions
     (NOT / total portfolio). This is the within-survey detection efficiency —
     how well the instrument found what it flew over.

  4. PoD: Marcellus-specific threshold from Thorpe et al. 2024 (Fig.5):
     PoD90 = 0.974 kg/h at 3.5 m/s wind. Wind sampled from PA climatology.

  5. Flight plan: grounded in Donahue et al. 2024 (Permian paper):
     ~130 marginal wells/day in Permian → ~100/day for PA (terrain adjustment),
     ~150 flyable days/year (vs Permian 195), giving ~23% annual coverage.

Sources:
  Thorpe et al. 2024, RSE 315:114435
  Donahue et al. 2024, preprint ES&T
"""

import numpy as np
import pandas as pd
import geopandas as gpd
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# ── Paths ────────────────────────────────────────────────────────────────────
DATA_MEAN   = '/tmp/pa_data/PAMarginalWellEmissions2023.gpkg'
DATA_ITER   = '/tmp/pa_data/PAMarginalWellEmissions2023_withEachIteration.gpkg'
OUT_DIR     = Path('BridgerResults')
OUT_DIR.mkdir(exist_ok=True)

# ── Constants (all grounded in literature) ────────────────────────────────────
POD90_MARCELLUS   = 0.974   # kg/h  Thorpe 2024 Fig.5 Marcellus mean
POD90_OVERALL     = 1.27    # kg/h  Thorpe 2024 Fig.5 overall 2023 avg
K_LOGISTIC        = 2.0     # logistic steepness (Conrad et al. 2023a)
WIND_REF          = 3.5     # m/s   reference wind for PoD90 calibration
WIND_EXP          = -0.3    # wind sensitivity coefficient (Thorpe 2024)
WIND_MEAN_PA      = 3.5     # m/s   PA seasonal flying-window mean
WIND_STD_PA       = 0.9     # m/s
WIND_MIN          = 1.0
WIND_MAX          = 6.0

WELLS_PER_DAY     = 100     # PA-adjusted (Permian 130 × 0.77 terrain factor)
FLYABLE_DAYS_YR   = 150     # PA vs Permian 195
ANNUAL_WELLS      = WELLS_PER_DAY * FLYABLE_DAYS_YR   # 15,000
GRID_DEG          = 0.1     # geographic survey cell size in degrees

N_MC              = 500     # Monte Carlo iterations (each draws random FEAST iter + wind)


# ═══════════════════════════════════════════════════════════════════════════
# 1. LOAD DATA
# ═══════════════════════════════════════════════════════════════════════════
print("Loading data …")
gdf_mean = gpd.read_file(DATA_MEAN)
gdf_iter = gpd.read_file(DATA_ITER)

N_WELLS = len(gdf_mean)
lats    = gdf_mean.Latitude.values
lons    = gdf_mean.Longitude.values
e_mean  = gdf_mean.mean_ch4_kgh.values          # 250-iteration mean per well

# Stack all 250 FEAST iterations into a (N_WELLS × 250) array
iter_cols = sorted([c for c in gdf_iter.columns if c.startswith('ch4_kgh_')])
N_ITER    = len(iter_cols)
print(f"  {N_WELLS:,} wells, {N_ITER} FEAST iterations. Building array …")
E_ITERS   = gdf_iter[iter_cols].values.astype(np.float32)   # shape (64624, 250)
print(f"  Array shape: {E_ITERS.shape}, size {E_ITERS.nbytes/1e6:.0f} MB")


# ═══════════════════════════════════════════════════════════════════════════
# 2. SPATIAL GRID CELLS  (survey routing units)
# ═══════════════════════════════════════════════════════════════════════════
lat_edges = np.arange(39.6, 42.4 + GRID_DEG, GRID_DEG)
lon_edges = np.arange(-80.6, -75.8 + GRID_DEG, GRID_DEG)

lat_idx = np.clip(np.digitize(lats, lat_edges) - 1, 0, len(lat_edges)-2)
lon_idx = np.clip(np.digitize(lons, lon_edges) - 1, 0, len(lon_edges)-2)
cell_id = lat_idx * 1000 + lon_idx

unique_cells = np.unique(cell_id)
N_CELLS      = len(unique_cells)

# Build cell-to-well index map
cell_to_wells = {c: np.where(cell_id == c)[0] for c in unique_cells}
cell_sizes    = np.array([len(cell_to_wells[c]) for c in unique_cells])

# Cell centroid lat/lon (for mapping)
cell_centroid_lat = np.array([lats[cell_to_wells[c]].mean() for c in unique_cells])
cell_centroid_lon = np.array([lons[cell_to_wells[c]].mean() for c in unique_cells])

# Days per cell at WELLS_PER_DAY throughput
days_per_cell   = cell_sizes / WELLS_PER_DAY
total_days_full = days_per_cell.sum()
years_full_pass = total_days_full / FLYABLE_DAYS_YR

# Fraction of wells in cells that fit in one year
annual_cell_budget = FLYABLE_DAYS_YR   # days available
# Select cells greedily (largest first) until budget used
# For survey planning we pick by total-days weighting
cum_days          = np.cumsum(np.sort(days_per_cell)[::-1])
n_cells_per_year  = np.searchsorted(cum_days, annual_cell_budget) + 1
n_cells_per_year  = min(n_cells_per_year, N_CELLS)
annual_wells_spatial = int(cell_sizes[np.argsort(cell_sizes)[::-1][:n_cells_per_year]].sum())
annual_coverage_pct  = 100 * annual_wells_spatial / N_WELLS

print(f"\nSpatial grid ({GRID_DEG}°): {N_CELLS} occupied cells")
print(f"  Wells/cell: mean={cell_sizes.mean():.0f}, median={np.median(cell_sizes):.0f}, "
      f"max={cell_sizes.max()}, min={cell_sizes.min()}")
print(f"  Days to survey ALL cells: {total_days_full:.0f} ({years_full_pass:.1f} years @ {FLYABLE_DAYS_YR} days/yr)")
print(f"  Cells coverable per year: {n_cells_per_year} → {annual_wells_spatial:,} wells ({annual_coverage_pct:.1f}%)")


# ═══════════════════════════════════════════════════════════════════════════
# 3. PoD MODEL
# ═══════════════════════════════════════════════════════════════════════════
def pod(e, wind, pod90=POD90_MARCELLUS):
    """Vectorised Bridger GML2.0 PoD (Thorpe/Conrad 2024)."""
    e90_wind = pod90 * np.exp(WIND_EXP * (wind - WIND_REF))
    e_safe   = np.maximum(e, 1e-9)
    return np.clip(1.0 / (1.0 + np.exp(-K_LOGISTIC * np.log(e_safe / e90_wind))), 0, 1)

def sample_wind():
    return float(np.clip(np.random.normal(WIND_MEAN_PA, WIND_STD_PA), WIND_MIN, WIND_MAX))


# ═══════════════════════════════════════════════════════════════════════════
# 4. EMISSION DISTRIBUTION ANALYSIS  (uses mean emissions)
# ═══════════════════════════════════════════════════════════════════════════
print("\n── Emission Distribution (mean_ch4_kgh) ──")
e_sorted_desc = np.sort(e_mean)[::-1]
cum_e         = np.cumsum(e_sorted_desc)
total_mean    = e_mean.sum()

conc = {}
for p in [1, 2, 5, 10, 20, 30, 50, 70]:
    n = max(1, int(N_WELLS * p / 100))
    conc[p] = dict(n_wells=n, emission_share=100*cum_e[n-1]/total_mean,
                   threshold_kgph=float(e_sorted_desc[n-1]))

print(f"  {'Top%':>6} {'N wells':>8} {'Emission%':>10} {'Min kgph':>10}")
for p, r in conc.items():
    print(f"  {p:>5}%  {r['n_wells']:>8,}  {r['emission_share']:>9.1f}%  {r['threshold_kgph']:>9.3f}")

# Wells above detection thresholds (mean)
thresh_stats = {}
for t in [0.3, 0.5, 0.974, 1.27, 2.0, 5.0]:
    n_a = (e_mean > t).sum()
    e_a = e_mean[e_mean > t].sum()
    thresh_stats[t] = dict(n_wells=int(n_a), pct_wells=100*n_a/N_WELLS,
                           emissions_kgph=float(e_a), pct_emissions=100*e_a/total_mean)

print(f"\n  Threshold  N_wells  %fleet  emissions  %total")
for t, r in thresh_stats.items():
    print(f"  {t:.3f} kg/h  {r['n_wells']:>7,}  {r['pct_wells']:>5.1f}%  "
          f"{r['emissions_kgph']:>9.0f}  {r['pct_emissions']:>5.1f}%")

# Per-iteration stats
iter_totals = E_ITERS.sum(axis=0).astype(float)
iter_above  = (E_ITERS > POD90_MARCELLUS).mean(axis=0) * 100
print(f"\n  Across {N_ITER} FEAST iterations:")
print(f"  Total emissions: {iter_totals.mean():.0f} ± {iter_totals.std():.0f} kg/h "
      f"[{iter_totals.min():.0f}–{iter_totals.max():.0f}]")
print(f"  % wells > PoD90 threshold: {iter_above.mean():.1f}% ± {iter_above.std():.1f}%")


# ═══════════════════════════════════════════════════════════════════════════
# 5. MONTE CARLO SURVEY SIMULATION
# ═══════════════════════════════════════════════════════════════════════════
# For each MC run:
#   a. Select survey cells (geographic cluster selection)
#   b. Draw random FEAST iteration  → snapshot emission state
#   c. Draw random wind
#   d. Apply PoD → Bernoulli detection
#   e. DENOMINATOR = sum(emissions of SURVEYED wells in that iteration)
#      NUMERATOR   = sum(emissions of DETECTED wells)
#   f. mitigation% = numerator / denominator

def run_survey_mc(target_coverage_pct, n_mc=N_MC, pod90=POD90_MARCELLUS):
    """
    Monte Carlo survey simulation using spatial cell selection.
    target_coverage_pct: fraction of WELLS to cover (0–1).
    """
    target_wells = int(N_WELLS * target_coverage_pct)

    # Pre-compute cell selection sets for speed:
    # Randomly select cells until we have ≈ target_wells wells
    results = []
    for _ in range(n_mc):
        # --- select cells in random order until target reached ---
        cell_order = np.random.permutation(N_CELLS)
        selected_well_idx = []
        for ci in cell_order:
            c = unique_cells[ci]
            selected_well_idx.extend(cell_to_wells[c])
            if len(selected_well_idx) >= target_wells:
                break
        selected_well_idx = np.array(selected_well_idx[:target_wells])

        # --- draw random FEAST iteration ---
        it    = np.random.randint(0, N_ITER)
        e_now = E_ITERS[selected_well_idx, it].astype(float)

        # --- wind ---
        w = sample_wind()

        # --- PoD + Bernoulli ---
        p_det  = pod(e_now, w, pod90)
        det    = np.random.binomial(1, p_det).astype(bool)

        surveyed_total   = e_now.sum()
        detected_total   = e_now[det].sum()
        mit_pct          = (100.0 * detected_total / surveyed_total
                            if surveyed_total > 0 else 0.0)

        # Top-5% capture (within surveyed wells)
        n5 = max(1, int(len(e_now) * 0.05))
        top5_thresh = np.partition(e_now, -n5)[-n5] if len(e_now) > n5 else 0
        top5_mask   = e_now >= top5_thresh
        top5_emit   = e_now[top5_mask].sum()
        top5_det    = e_now[det & top5_mask].sum()
        top5_cap    = 100.0 * top5_det / top5_emit if top5_emit > 0 else 0.0

        # Average emission of detected wells vs all surveyed
        avg_det  = float(e_now[det].mean()) if det.any() else 0.0
        avg_surv = float(e_now.mean()) if len(e_now) > 0 else 0.0

        results.append(dict(
            n_surveyed=len(selected_well_idx),
            n_detected=int(det.sum()),
            surveyed_total_kgph=float(surveyed_total),
            detected_total_kgph=float(detected_total),
            mitigation_pct=float(mit_pct),       # detected / SURVEYED (correct)
            avg_pod=float(p_det.mean()),
            wind_ms=float(w),
            feast_iter=int(it),
            top5_capture_pct=float(top5_cap),
            avg_det_emission=avg_det,
            avg_surv_emission=avg_surv,
            enrichment=avg_det/avg_surv if avg_surv > 0 else 0,
        ))

    df = pd.DataFrame(results)
    summary = {}
    for col in ['mitigation_pct','detected_total_kgph','surveyed_total_kgph',
                'n_detected','avg_pod','top5_capture_pct','enrichment','wind_ms']:
        summary[col] = dict(
            mean=float(df[col].mean()),
            std=float(df[col].std()),
            p5=float(df[col].quantile(0.05)),
            p25=float(df[col].quantile(0.25)),
            median=float(df[col].median()),
            p75=float(df[col].quantile(0.75)),
            p95=float(df[col].quantile(0.95)),
        )
    return summary, df


# ── Coverage scenarios ────────────────────────────────────────────────────
print(f"\n── Monte Carlo Survey ({N_MC} runs each) ──")
print(f"  NOTE: mitigation% = detected / surveyed-wells'-total-emissions\n")

coverage_scenarios = [
    ("One season (Q2+Q3, ~75 days)",         0.116),
    ("1 year moderate (150 days, 100/day)",   ANNUAL_WELLS / N_WELLS),
    ("1 year ambitious (150 days, 130/day)",  130*150/N_WELLS),
    ("Full portfolio pass",                   0.999),
]

all_mc = {}
print(f"  {'Scenario':<45} {'Cov%':>5} | {'Within-survey':>16} | {'Mean det.':>12} | {'Top5-cap':>10} | {'Enrichment':>11}")
print(f"  {'':45} {'':5} | {'mitigation %':>16} | {'kg/h':>12} | {'%':>10} | {'factor':>11}")
print(f"  {'-'*105}")

for name, cov in coverage_scenarios:
    summ, df = run_survey_mc(cov, n_mc=N_MC)
    all_mc[name] = {'coverage_pct': cov*100, 'summary': summ, 'df': df}
    m = summ['mitigation_pct']
    d = summ['detected_total_kgph']
    t = summ['top5_capture_pct']
    en = summ['enrichment']
    print(f"  {name:<45} {cov*100:>4.0f}%"
          f" | {m['mean']:>6.1f}% [{m['p5']:.0f}–{m['p95']:.0f}]"
          f" | {d['mean']:>10.0f} ± {d['std']:.0f}"
          f" | {t['mean']:>8.1f}%"
          f" | {en['mean']:>9.1f}x")


# ── Realistic 1-year scenario deep dive ──────────────────────────────────
FOCAL_NAME = "1 year moderate (150 days, 100/day)"
focal      = all_mc[FOCAL_NAME]
summ_f     = focal['summary']
df_f       = focal['df']
cov_f      = focal['coverage_pct']

print(f"\n── Deep dive: '{FOCAL_NAME}' ──")
print(f"  Coverage: {cov_f:.1f}% of fleet = {int(N_WELLS*cov_f/100):,} wells")
print()
print(f"  Within-survey mitigation (detected / surveyed-total):")
m = summ_f['mitigation_pct']
print(f"    Mean:    {m['mean']:.1f}%")
print(f"    Std:     {m['std']:.1f}%")
print(f"    90% CI:  [{m['p5']:.1f}% – {m['p95']:.1f}%]")
print()
print(f"  Emissions detected (kg/h):")
d = summ_f['detected_total_kgph']
print(f"    Mean ± Std: {d['mean']:.0f} ± {d['std']:.0f} kg/h")
print(f"    90% CI:     [{d['p5']:.0f} – {d['p95']:.0f}] kg/h")
print()
print(f"  Average PoD: {summ_f['avg_pod']['mean']:.3f}")
print(f"  Detection enrichment: {summ_f['enrichment']['mean']:.1f}x  "
      f"(detected wells emit {summ_f['enrichment']['mean']:.1f}× average of surveyed)")
print(f"  Top-5% emitter capture: {summ_f['top5_capture_pct']['mean']:.1f}% of top-5% emissions detected")

# Most-emitters vs top-emitters
n_surveyed_mean = df_f['n_surveyed'].mean()
n_detected_mean = df_f['n_detected'].mean()
pct_fleet_det   = 100 * n_detected_mean / N_WELLS
print()
print(f"  ── 'Most emitters' vs 'Top emitters' ──")
print(f"  Wells surveyed per run:  {n_surveyed_mean:,.0f} ({cov_f:.1f}% of fleet)")
print(f"  Wells detected per run:  {n_detected_mean:,.0f} ({pct_fleet_det:.1f}% of total fleet)")
print(f"  → 'Most emitters captured'?  {'YES' if pct_fleet_det > 50 else f'NO – only {pct_fleet_det:.1f}% of all wells detected'}")
print(f"  → 'Top emitters captured'?   YES – detected wells emit {summ_f['enrichment']['mean']:.1f}× average;")
print(f"                                top-5% emitters yield {summ_f['top5_capture_pct']['mean']:.0f}% capture in surveyed zones")


# ═══════════════════════════════════════════════════════════════════════════
# 6. FLIGHT PLAN SUMMARY
# ═══════════════════════════════════════════════════════════════════════════
print(f"\n── Flight Plan (Permian-grounded) ──")
print(f"  Permian benchmark (Donahue 2024):  130 marginal wells/day")
print(f"  PA adjustment (terrain, weather):  100 wells/day, 150 flyable days/yr")
print(f"  PA wells/year:  {ANNUAL_WELLS:,}  ({100*ANNUAL_WELLS/N_WELLS:.1f}% coverage)")
print(f"  Full portfolio: {total_days_full:.0f} days ({years_full_pass:.1f} years)")
print(f"  Quarterly schedule (PA):")
for q, days, note in [("Q1 Jan–Mar","38","(winter, fewest flyable)"),
                       ("Q2 Apr–Jun","40","(recommended – mild winds)"),
                       ("Q3 Jul–Sep","35","(good – thunderstorm gaps)"),
                       ("Q4 Oct–Dec","37","(declining flyability)")]:
    wells = int(int(days) * WELLS_PER_DAY)
    print(f"    {q}: {days} days × {WELLS_PER_DAY} wells/day = {wells:,} wells  {note}")


# ═══════════════════════════════════════════════════════════════════════════
# 7. PLOTS
# ═══════════════════════════════════════════════════════════════════════════
print("\n── Generating plots …")
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
fig.suptitle("Bridger LiDAR Survey — PA Marginal Wells\n"
             "Grounded in Thorpe et al. 2024 & Donahue et al. 2024",
             fontsize=13, fontweight='bold')

# ── (A) Emission distribution (mean) ──
ax = axes[0, 0]
e_pos = e_mean[e_mean > 0]
ax.hist(e_pos[e_pos <= 3.0], bins=80, color='steelblue', alpha=0.8, density=True)
ax.axvline(POD90_MARCELLUS, color='red',   lw=2, ls='--', label=f'Marcellus PoD90={POD90_MARCELLUS} kg/h')
ax.axvline(POD90_OVERALL,   color='orange',lw=2, ls=':',  label=f'Overall PoD90={POD90_OVERALL} kg/h')
ax.set_xlabel('Mean Emission (kg/h)')
ax.set_ylabel('Density')
ax.set_title('(A) Well Emission Distribution\n(mean across 250 FEAST iterations)')
ax.set_xlim(0, 3.0)
ax.legend(fontsize=8)

# ── (B) Lorenz curve ──
ax = axes[0, 1]
sorted_e = np.sort(e_mean)
cum_wells = np.linspace(0, 100, len(sorted_e))
cum_emis  = np.cumsum(sorted_e) / sorted_e.sum() * 100
ax.plot(cum_wells, cum_emis, 'steelblue', lw=2)
ax.plot([0,100],[0,100],'k--', alpha=0.4, label='Perfect equality')
ax.fill_between(cum_wells, cum_emis, cum_wells, alpha=0.15, color='steelblue')
# Mark 50% emission line
idx50 = np.searchsorted(cum_emis, 50)
ax.axhline(50, color='red', ls=':', lw=1.5)
ax.axvline(cum_wells[idx50], color='red', ls=':', lw=1.5,
           label=f'Bottom {cum_wells[idx50]:.0f}% wells → 50% emissions')
ax.set_xlabel('Cumulative % of Wells (smallest first)')
ax.set_ylabel('Cumulative % of Emissions')
ax.set_title('(B) Emission Lorenz Curve\n(concentration of emissions)')
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

# ── (C) PoD curve ──
ax = axes[0, 2]
e_range = np.logspace(-2, 2, 500)
for ws, col, ls in [(1.0,'blue','-.'), (2.0,'cyan',':'),
                    (3.5,'green','-'), (5.0,'orange','--'), (6.0,'red','-.')]:
    ax.semilogx(e_range, pod(e_range, ws), color=col, ls=ls, label=f'{ws} m/s')
ax.axhline(0.9, color='gray', ls=':', lw=1)
ax.axvline(POD90_MARCELLUS, color='red', ls='--', lw=1.5, label=f'Marcellus PoD90')
ax.set_xlabel('Emission Rate (kg/h)')
ax.set_ylabel('Probability of Detection')
ax.set_title(f'(C) Bridger PoD vs Wind Speed\n(Marcellus, PoD90={POD90_MARCELLUS} kg/h)')
ax.legend(fontsize=7, ncol=2)
ax.grid(alpha=0.3)
ax.set_xlim(0.01, 100)

# ── (D) Within-survey mitigation distribution ──
ax = axes[1, 0]
colors_ = ['#2196F3','#4CAF50','#FF9800','#9C27B0']
for (name, cov), c in zip(coverage_scenarios, colors_):
    df_sc = all_mc[name]['df']
    ax.hist(df_sc['mitigation_pct'], bins=40, alpha=0.55,
            label=f"{cov*100:.0f}% cov (μ={df_sc['mitigation_pct'].mean():.0f}%)",
            color=c, density=True)
ax.set_xlabel('Within-survey mitigation %\n(detected / surveyed-total emissions)')
ax.set_ylabel('Density')
ax.set_title('(D) Within-Survey Detection Rate\n(correct denominator: surveyed emissions)')
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

# ── (E) Emissions detected distribution (1-year scenario) ──
ax = axes[1, 1]
ax.hist(df_f['detected_total_kgph'], bins=50, color='steelblue', alpha=0.75)
ax.axvline(df_f['detected_total_kgph'].mean(), color='red', lw=2,
           label=f"Mean: {df_f['detected_total_kgph'].mean():.0f} kg/h")
ax.axvline(df_f['detected_total_kgph'].quantile(0.05), color='orange', lw=2, ls='--',
           label=f"p5: {df_f['detected_total_kgph'].quantile(0.05):.0f} kg/h")
ax.axvline(df_f['detected_total_kgph'].quantile(0.95), color='orange', lw=2, ls='--',
           label=f"p95: {df_f['detected_total_kgph'].quantile(0.95):.0f} kg/h")
ax.set_xlabel('Emissions Detected (kg/h)')
ax.set_ylabel('Count')
ax.set_title(f'(E) Detected Emissions Distribution\n1-year survey ({cov_f:.0f}% coverage, {N_MC} MC runs)')
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

# ── (F) Spatial well map + cell grid ──
ax = axes[1, 2]
# Thin scatter of all wells
ax.scatter(lons, lats, s=0.1, c='lightgray', alpha=0.3, rasterized=True)
# Colour cells by mean emission
cell_mean_e = np.array([e_mean[cell_to_wells[c]].mean() for c in unique_cells])
sc = ax.scatter(cell_centroid_lon, cell_centroid_lat, c=cell_mean_e,
                cmap='YlOrRd', s=cell_sizes/5, alpha=0.8, vmin=0, vmax=2)
plt.colorbar(sc, ax=ax, label='Cell mean emission (kg/h)', shrink=0.7)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title(f'(F) PA Marginal Wells — {N_CELLS} Survey Grid Cells\n(dot size ∝ well count, colour = mean emission)')

plt.tight_layout()
out_fig = OUT_DIR / 'bridger_pa_full_analysis.png'
plt.savefig(out_fig, dpi=150, bbox_inches='tight')
print(f"  Saved: {out_fig}")


# ── Extra: single-iteration distribution plot ──
fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5))
fig2.suptitle("FEAST Emission Variability: Mean vs Single Iteration",
              fontsize=12, fontweight='bold')

ax2 = axes2[0]
e_i0 = E_ITERS[:, 0].astype(float)
e_i0_pos = e_i0[e_i0 > 0]
ax2.hist(e_i0_pos[e_i0_pos <= 10], bins=80, color='tomato', alpha=0.7,
         density=True, label=f'One iteration (max={e_i0.max():.0f} kg/h)')
e_mean_pos = e_mean[e_mean > 0]
ax2.hist(e_mean_pos[e_mean_pos <= 10], bins=80, color='steelblue', alpha=0.5,
         density=True, label='Mean of 250 iterations')
ax2.axvline(POD90_MARCELLUS, color='red', ls='--', lw=1.5, label='Marcellus PoD90')
ax2.set_xlabel('Emission (kg/h)')
ax2.set_ylabel('Density')
ax2.set_title('Emission distribution:\nTime-averaged vs instantaneous snapshot')
ax2.legend(fontsize=9)
ax2.set_xlim(0, 10)

ax2 = axes2[1]
# Show how total emissions vary across iterations
ax2.hist(iter_totals, bins=40, color='mediumseagreen', alpha=0.75)
ax2.axvline(iter_totals.mean(), color='red', lw=2,
            label=f'Mean: {iter_totals.mean():.0f} kg/h')
ax2.axvline(np.percentile(iter_totals,5), color='orange', ls='--', lw=1.5,
            label=f'p5: {np.percentile(iter_totals,5):.0f} kg/h')
ax2.axvline(np.percentile(iter_totals,95), color='orange', ls='--', lw=1.5,
            label=f'p95: {np.percentile(iter_totals,95):.0f} kg/h')
ax2.set_xlabel('Total Portfolio Emissions (kg/h)')
ax2.set_ylabel('Count (iterations)')
ax2.set_title('Total PA marginal well emissions\nacross 250 FEAST stochastic iterations')
ax2.legend(fontsize=9)

plt.tight_layout()
out_fig2 = OUT_DIR / 'feast_iteration_variability.png'
plt.savefig(out_fig2, dpi=150, bbox_inches='tight')
print(f"  Saved: {out_fig2}")


# ═══════════════════════════════════════════════════════════════════════════
# 8. SAVE JSON RESULTS
# ═══════════════════════════════════════════════════════════════════════════
output = {
    'metadata': {
        'n_wells': N_WELLS, 'n_feast_iterations': N_ITER,
        'pod90_marcellus_kgph': POD90_MARCELLUS,
        'wells_per_day_pa': WELLS_PER_DAY,
        'flyable_days_pa_year': FLYABLE_DAYS_YR,
        'annual_coverage_pct': round(100*ANNUAL_WELLS/N_WELLS, 1),
        'years_for_full_portfolio': round(years_full_pass, 1),
        'n_survey_cells': N_CELLS,
        'survey_denominator': 'surveyed_wells_total_emissions (NOT full portfolio)',
    },
    'emission_concentration': {str(k): v for k,v in conc.items()},
    'threshold_stats': {str(k): v for k,v in thresh_stats.items()},
    'feast_iteration_variability': {
        'total_kgph_mean': float(iter_totals.mean()),
        'total_kgph_std': float(iter_totals.std()),
        'total_kgph_p5': float(np.percentile(iter_totals, 5)),
        'total_kgph_p95': float(np.percentile(iter_totals, 95)),
    },
    'survey_scenarios': {
        name: {'coverage_pct': v['coverage_pct'],
               'summary': {k: vv for k, vv in v['summary'].items()}}
        for name, v in all_mc.items()
    },
    'key_findings': {
        'pct_wells_detected_1yr': round(pct_fleet_det, 1),
        'within_survey_mitigation_pct_mean': round(summ_f['mitigation_pct']['mean'], 1),
        'within_survey_mitigation_pct_ci90': [
            round(summ_f['mitigation_pct']['p5'], 1),
            round(summ_f['mitigation_pct']['p95'], 1)],
        'detection_enrichment_factor': round(summ_f['enrichment']['mean'], 1),
        'top5_capture_pct': round(summ_f['top5_capture_pct']['mean'], 1),
        'answer': (
            'Bridger captures TOP EMITTERS, not most emitters. '
            f'Only {pct_fleet_det:.1f}% of all wells detected in a 1-year survey. '
            f'But detected wells emit {summ_f["enrichment"]["mean"]:.1f}x the average '
            'of surveyed wells — the instrument disproportionately finds large sources. '
            f'Within each surveyed zone, {summ_f["mitigation_pct"]["mean"]:.0f}% of '
            'that zone\'s estimated emissions are detected on average.'
        ),
    },
}

out_json = OUT_DIR / 'bridger_pa_full_analysis.json'
with open(out_json, 'w') as f:
    json.dump(output, f, indent=2, default=str)
print(f"\n  Results saved: {out_json}")
print("\n" + "="*70)
print("DONE")
print("="*70)
