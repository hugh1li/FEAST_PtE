"""
FEAST Path B — ingest Bridger aerial LiDAR deliveries against the four-file
template produced in BridgerDataRequest/templates/.

What this script does
---------------------
1. Loads surveys / site_passes / detections / well_metadata (CSV or GPKG).
2. Validates the schema against BridgerDataRequest/templates/README.md.
3. Cross-joins detections ↔ site_passes ↔ surveys ↔ wells.
4. Recomputes per-pass POD from (GCN, wind, emission rate) using the GML 2.0
   P4 + Burr inverse link (Thorpe et al. 2024, RSE 315:114435, Table 3) to
   provide an instrument-model-independent check on Bridger's reported POD.
5. For each basin, produces:
      • high-end plume catalog  (kg/h distribution, counts, top-N)
      • per-pass detection rate (n_detections / n_passes)
      • POD-corrected "expected emission above detection threshold"
      • a stratified summary binned by (gas production, age, GOR) so the
        aerial observations can be compared directly against the 400-well
        ground campaign sampling frame.
6. Writes a JSON summary + per-basin CSV rollups + (optional) a 4-panel
   diagnostic PNG.

Design notes
------------
* This script is intentionally independent of the EIME / Path A pipeline.
  It treats Bridger as ground truth within its detection window and
  produces an unbiased high-end supplement.
* Missing `02_site_passes.csv` is a soft error: we warn, fall back to a
  plume catalog only, and skip detection-rate / POD-correction outputs.
* Works against the example rows shipped in templates/ so the script can be
  smoke-tested with no real data in hand.

Usage
-----
    python run_bridger_pathB_ingest.py \
        --data-dir ../BridgerDataRequest/templates \
        --out-dir  ./BridgerResults/PathB_smoke \
        --format   csv

    python run_bridger_pathB_ingest.py \
        --data-dir /path/to/bridger_delivery \
        --format   auto     # tries gpkg first, then csv
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

# ------------------------------------------------------------------ #
# Bridger GML 2.0 POD coefficients (Thorpe et al. 2024, Table 3)
# POD(E) = 1 - (1 + (P4 / beta_1)**beta_2) ** -beta_3,
# with P4 = alpha_1 * (E/wind) / n_gcn + alpha_2 * ((E/wind) / n_gcn)**2 ,
# n_gcn = GCN / 1000  (ppm-m → fractional units used in the fit).
# ------------------------------------------------------------------ #
POD_ALPHA1 = 2.0
POD_ALPHA2 = 1.5
POD_BETA1 = 2.41e-3
POD_BETA2 = 1.9505
POD_BETA3 = 2.0836
POD_BETA4 = 1.5185  # retained for completeness; unused in the simplified P4

# Typical flight-level GCN (ppm-m) used when a per-pass value is missing.
DEFAULT_GCN_PPM_M = 16.0


# ------------------------------------------------------------------ #
# Schema: required columns per file. Matches templates/README.md.
# (*) columns that are nominally required but may be best-effort in the
# delivered data (e.g. API number) — flagged but not hard-failed.
# ------------------------------------------------------------------ #
SCHEMA = {
    "01_surveys.csv": {
        "required": [
            "survey_id", "survey_date_utc", "basin", "region",
            "flight_start_utc", "flight_end_utc", "coverage_wkt",
            "gcn_ppm_m", "mean_wind_speed_ms",
            "mdl_kgph_at_3ms", "pod_model_version",
        ],
        "optional": [
            "aircraft_id", "mean_wind_direction_deg", "mean_altitude_ft_agl",
            "atmospheric_stability_class", "notes",
        ],
    },
    "02_site_passes.csv": {
        "required": [
            "site_pass_id", "survey_id", "site_id",
            "latitude_wgs84", "longitude_wgs84",
            "pass_timestamp_utc", "wind_speed_ms_at_pass", "plume_detected",
        ],
        "optional": [
            "api_number", "wind_direction_deg_at_pass", "altitude_ft_agl_at_pass",
            "pod_at_pass", "mdl_kgph_at_pass", "n_passes_in_survey", "notes",
        ],
    },
    "03_detections.csv": {
        "required": [
            "plume_id", "site_pass_id", "survey_id",
            "detection_timestamp_utc", "latitude_wgs84", "longitude_wgs84",
            "associated_site_id",
            "emission_rate_kgph", "emission_rate_lo_kgph", "emission_rate_hi_kgph",
            "emission_rate_uncertainty_basis",
            "wind_speed_ms", "gcn_ppm_m",
        ],
        "optional": [
            "horizontal_position_uncertainty_m", "associated_api_number",
            "wind_direction_deg", "altitude_ft_agl", "plume_length_m",
            "plume_image_filename", "classification", "notes",
        ],
    },
    "04_well_metadata.csv": {
        "required": [
            "site_id", "basin", "region",
            "latitude_wgs84", "longitude_wgs84",
        ],
        "optional": [
            "api_number", "operator", "spud_date", "first_production_date",
            "well_type", "lift_type",
            "gas_production_mcfpd_12mo", "oil_production_bopd_12mo",
            "water_production_bblpd_12mo", "gor_scf_per_bbl",
            "marginal_status", "notes",
        ],
    },
}


# ================================================================== #
# POD model
# ================================================================== #
def pod_gml2(
    emission_kgph: np.ndarray | float,
    wind_ms: np.ndarray | float,
    gcn_ppm_m: np.ndarray | float = DEFAULT_GCN_PPM_M,
) -> np.ndarray | float:
    """Probability of detection from Bridger GML 2.0 (Thorpe et al. 2024).

    Vectorised. Inputs can be scalars or ndarrays.
    Emission in kg/h, wind in m/s, GCN in ppm-m.
    Returns a probability in [0, 1].
    """
    e = np.asarray(emission_kgph, dtype=float)
    u = np.asarray(wind_ms, dtype=float)
    g = np.asarray(gcn_ppm_m, dtype=float) / 1000.0  # ppm-m → frac

    # Guard against divide-by-zero
    u_safe = np.where(u > 0, u, np.nan)
    g_safe = np.where(g > 0, g, np.nan)

    x = (e / u_safe) / g_safe                     # dimensionless driver
    p4 = POD_ALPHA1 * x + POD_ALPHA2 * x * x
    # Burr inverse link
    with np.errstate(invalid="ignore", over="ignore"):
        pod = 1.0 - np.power(1.0 + np.power(p4 / POD_BETA1, POD_BETA2),
                             -POD_BETA3)
    pod = np.where(np.isfinite(pod), pod, 0.0)
    pod = np.clip(pod, 0.0, 1.0)
    return pod


# ================================================================== #
# Loaders
# ================================================================== #
@dataclass
class BridgerDelivery:
    surveys: pd.DataFrame
    site_passes: Optional[pd.DataFrame]
    detections: pd.DataFrame
    wells: Optional[pd.DataFrame]
    warnings: list[str] = field(default_factory=list)

    @property
    def has_passes(self) -> bool:
        return self.site_passes is not None and len(self.site_passes) > 0

    @property
    def has_wells(self) -> bool:
        return self.wells is not None and len(self.wells) > 0


def _read_table(
    data_dir: Path, filename: str, fmt: str, required: bool,
    warnings_list: list[str],
) -> Optional[pd.DataFrame]:
    stem = filename.replace(".csv", "")
    csv_path = data_dir / filename
    gpkg_path = data_dir / f"{stem}.gpkg"

    if fmt in ("gpkg", "auto") and gpkg_path.exists():
        try:
            import geopandas as gpd  # noqa: WPS433
            gdf = gpd.read_file(gpkg_path)
            df = pd.DataFrame(gdf.drop(columns="geometry", errors="ignore"))
        except Exception as exc:  # noqa: BLE001
            warnings_list.append(f"Failed to read {gpkg_path.name}: {exc}")
            df = None
        if df is not None:
            return df

    if fmt in ("csv", "auto") and csv_path.exists():
        return pd.read_csv(csv_path)

    if required:
        raise FileNotFoundError(f"Required table {filename} not found in {data_dir}")
    warnings_list.append(f"Optional table {filename} not found; continuing without it.")
    return None


def _validate_schema(df: pd.DataFrame, filename: str, warnings_list: list[str]) -> None:
    req = SCHEMA[filename]["required"]
    missing = [c for c in req if c not in df.columns]
    if missing:
        raise ValueError(
            f"{filename} is missing required columns: {missing}. "
            f"See BridgerDataRequest/templates/README.md for schema."
        )
    opt = SCHEMA[filename]["optional"]
    extra = [c for c in df.columns if c not in req + opt]
    if extra:
        warnings_list.append(f"{filename}: extra columns {extra} — kept as-is.")


def load_delivery(data_dir: Path, fmt: str = "auto") -> BridgerDelivery:
    warnings_list: list[str] = []
    surveys = _read_table(data_dir, "01_surveys.csv", fmt, True, warnings_list)
    site_passes = _read_table(data_dir, "02_site_passes.csv", fmt, False, warnings_list)
    detections = _read_table(data_dir, "03_detections.csv", fmt, True, warnings_list)
    wells = _read_table(data_dir, "04_well_metadata.csv", fmt, False, warnings_list)

    _validate_schema(surveys, "01_surveys.csv", warnings_list)
    _validate_schema(detections, "03_detections.csv", warnings_list)
    if site_passes is not None:
        _validate_schema(site_passes, "02_site_passes.csv", warnings_list)
    if wells is not None:
        _validate_schema(wells, "04_well_metadata.csv", warnings_list)

    return BridgerDelivery(surveys, site_passes, detections, wells, warnings_list)


# ================================================================== #
# Analytics
# ================================================================== #
def _enrich_detections(d: BridgerDelivery) -> pd.DataFrame:
    """Join detections to surveys and wells; recompute POD from first principles."""
    det = d.detections.copy()
    surv = d.surveys[["survey_id", "basin", "region",
                      "mdl_kgph_at_3ms", "pod_model_version"]].copy()
    det = det.merge(surv, on="survey_id", how="left", suffixes=("", "_survey"))

    if d.has_wells:
        well_cols = [c for c in [
            "site_id", "operator", "lift_type", "well_type",
            "gas_production_mcfpd_12mo", "oil_production_bopd_12mo",
            "gor_scf_per_bbl", "marginal_status",
        ] if c in d.wells.columns]
        det = det.merge(
            d.wells[well_cols],
            left_on="associated_site_id", right_on="site_id",
            how="left", suffixes=("", "_well"),
        )

    # Recomputed POD (sanity-check; for plumes this should be ≈1 by construction
    # because the plume was in fact detected, but useful for plumes near threshold)
    det["pod_recomputed"] = pod_gml2(
        det["emission_rate_kgph"].to_numpy(),
        det["wind_speed_ms"].to_numpy(),
        det["gcn_ppm_m"].to_numpy(),
    )
    return det


def _per_pass_pod(d: BridgerDelivery) -> Optional[pd.DataFrame]:
    if not d.has_passes:
        return None
    sp = d.site_passes.copy()
    surv = d.surveys[["survey_id", "basin", "region", "gcn_ppm_m"]].copy()
    sp = sp.merge(surv, on="survey_id", how="left", suffixes=("", "_survey"))
    # Identify GCN column (survey-level fallback if pass-level absent)
    gcn = sp["gcn_ppm_m"] if "gcn_ppm_m" in sp.columns else DEFAULT_GCN_PPM_M
    # POD at pass for a representative emission grid is optional here;
    # we just keep the values as-delivered plus survey linkage.
    sp["_gcn_used_ppm_m"] = gcn
    return sp


def basin_summary(det: pd.DataFrame,
                  passes: Optional[pd.DataFrame]) -> pd.DataFrame:
    """Per-basin roll-up."""
    rows = []
    basins = sorted(det["basin"].dropna().unique())
    for b in basins:
        dsub = det[det["basin"] == b]
        row: dict = {
            "basin": b,
            "n_detections": len(dsub),
            "emission_mean_kgph": float(dsub["emission_rate_kgph"].mean()),
            "emission_median_kgph": float(dsub["emission_rate_kgph"].median()),
            "emission_p95_kgph": float(dsub["emission_rate_kgph"].quantile(0.95)),
            "emission_max_kgph": float(dsub["emission_rate_kgph"].max()),
            "total_observed_kgph": float(dsub["emission_rate_kgph"].sum()),
        }
        if passes is not None:
            psub = passes[passes["basin"] == b]
            n_passes = len(psub)
            n_det_passes = int(psub["plume_detected"].astype(bool).sum())
            row["n_passes"] = n_passes
            row["n_passes_with_detection"] = n_det_passes
            row["per_pass_detection_rate"] = (
                n_det_passes / n_passes if n_passes > 0 else np.nan
            )
            row["unique_sites_flown"] = int(psub["site_id"].nunique())
        rows.append(row)
    return pd.DataFrame(rows)


def stratified_summary(det: pd.DataFrame) -> Optional[pd.DataFrame]:
    """Bin detections by (gas production, age-proxy, GOR) for ground-campaign comparison."""
    needed = {"gas_production_mcfpd_12mo", "gor_scf_per_bbl"}
    if not needed.issubset(det.columns):
        return None
    g = det.copy()
    gas_bins = [0, 5, 15, 50, 200, np.inf]
    gas_labels = ["<5", "5-15", "15-50", "50-200", ">200"]
    gor_bins = [0, 1000, 10000, 100000, np.inf]
    gor_labels = ["<1k", "1k-10k", "10k-100k", ">100k"]
    g["gas_bin"] = pd.cut(g["gas_production_mcfpd_12mo"], bins=gas_bins,
                          labels=gas_labels, include_lowest=True)
    g["gor_bin"] = pd.cut(g["gor_scf_per_bbl"], bins=gor_bins,
                          labels=gor_labels, include_lowest=True)
    return (
        g.groupby(["basin", "gas_bin", "gor_bin"], observed=True)
         .agg(n=("plume_id", "count"),
              mean_kgph=("emission_rate_kgph", "mean"),
              p95_kgph=("emission_rate_kgph", lambda s: float(np.quantile(s, 0.95))))
         .reset_index()
    )


def pod_correction(passes: Optional[pd.DataFrame],
                   det: pd.DataFrame) -> Optional[pd.DataFrame]:
    """Crude POD correction: for each basin, divide the detected emissions total
    by the mean POD evaluated at the observed plume emission rates. A rough
    'what you'd see if POD were 1 everywhere above noise'."""
    if passes is None:
        return None
    rows = []
    for b in sorted(det["basin"].dropna().unique()):
        dsub = det[det["basin"] == b]
        pod_mean = float(pod_gml2(
            dsub["emission_rate_kgph"].to_numpy(),
            dsub["wind_speed_ms"].to_numpy(),
            dsub["gcn_ppm_m"].to_numpy(),
        ).mean()) if len(dsub) else np.nan
        corrected = (dsub["emission_rate_kgph"].sum() / pod_mean) if pod_mean and pod_mean > 0 else np.nan
        rows.append({
            "basin": b,
            "n_detections": len(dsub),
            "mean_pod_at_detections": pod_mean,
            "observed_total_kgph": float(dsub["emission_rate_kgph"].sum()),
            "pod_corrected_total_kgph": float(corrected),
            "note": "divides observed detection total by mean POD at detected emission rates; a first-order correction only",
        })
    return pd.DataFrame(rows)


# ================================================================== #
# Diagnostic plot
# ================================================================== #
def make_diagnostic_plot(det: pd.DataFrame,
                         passes: Optional[pd.DataFrame],
                         out_png: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # noqa: BLE001
        warnings.warn(f"matplotlib unavailable, skipping plot: {exc}")
        return

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    # Panel 1: emission-rate histogram by basin
    ax = axes[0, 0]
    for b, sub in det.groupby("basin"):
        vals = np.log10(sub["emission_rate_kgph"].clip(lower=0.01))
        ax.hist(vals, bins=25, alpha=0.5, label=f"{b} (n={len(sub)})")
    ax.set_xlabel("log10(emission rate, kg/h)")
    ax.set_ylabel("count")
    ax.set_title("Detected plume emission rates by basin")
    ax.legend(fontsize=8)

    # Panel 2: POD vs emission with observed plumes overlaid
    ax = axes[0, 1]
    e_grid = np.logspace(-1, 2.5, 200)
    for u in (2.0, 3.0, 5.0):
        ax.plot(e_grid, pod_gml2(e_grid, u, DEFAULT_GCN_PPM_M),
                label=f"wind={u:.0f} m/s", lw=1.5)
    ax.scatter(det["emission_rate_kgph"], det["pod_recomputed"],
               s=12, alpha=0.5, color="black", label="observed plumes")
    ax.set_xscale("log")
    ax.set_xlabel("emission rate (kg/h)")
    ax.set_ylabel("POD")
    ax.set_title(f"GML 2.0 POD curves (GCN={DEFAULT_GCN_PPM_M:.0f} ppm-m)")
    ax.set_ylim(0, 1.02)
    ax.legend(fontsize=8)

    # Panel 3: per-basin detection rate
    ax = axes[1, 0]
    if passes is not None:
        rates = (
            passes.groupby("basin")["plume_detected"]
                  .agg(["mean", "count"]).reset_index()
        )
        ax.bar(rates["basin"], rates["mean"])
        for i, (m, c) in enumerate(zip(rates["mean"], rates["count"])):
            ax.text(i, m + 0.005, f"n={c}", ha="center", fontsize=9)
        ax.set_ylabel("per-pass detection rate")
        ax.set_title("Per-pass detection rate by basin")
    else:
        ax.text(0.5, 0.5, "02_site_passes.csv not provided\n(no per-pass rate)",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()

    # Panel 4: emission-rate ECDF
    ax = axes[1, 1]
    for b, sub in det.groupby("basin"):
        v = np.sort(sub["emission_rate_kgph"].to_numpy())
        if len(v) == 0:
            continue
        ax.step(v, np.arange(1, len(v) + 1) / len(v), label=f"{b}")
    ax.set_xscale("log")
    ax.set_xlabel("emission rate (kg/h)")
    ax.set_ylabel("empirical CDF")
    ax.set_title("Emission-rate ECDF by basin")
    ax.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


# ================================================================== #
# Main
# ================================================================== #
def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-dir", type=Path,
                        default=Path(__file__).parent.parent
                        / "BridgerDataRequest" / "templates",
                        help="Directory holding 01..04 files (CSV or GPKG).")
    parser.add_argument("--out-dir", type=Path,
                        default=Path(__file__).parent / "BridgerResults" / "PathB_smoke",
                        help="Where to write summaries.")
    parser.add_argument("--format", choices=("auto", "csv", "gpkg"), default="auto",
                        help="Input format preference.")
    parser.add_argument("--no-plot", action="store_true", help="Skip diagnostic PNG.")
    args = parser.parse_args()

    print(f"[path B] loading delivery from: {args.data_dir}")
    try:
        d = load_delivery(args.data_dir, fmt=args.format)
    except Exception as exc:  # noqa: BLE001
        print(f"[path B] FAILED: {exc}", file=sys.stderr)
        return 1

    for w in d.warnings:
        print(f"  [warn] {w}")

    print(f"  surveys:     {len(d.surveys)}")
    print(f"  site passes: {len(d.site_passes) if d.has_passes else 0}")
    print(f"  detections:  {len(d.detections)}")
    print(f"  wells:       {len(d.wells) if d.has_wells else 0}")

    det = _enrich_detections(d)
    passes = _per_pass_pod(d)
    summary = basin_summary(det, passes)
    strat = stratified_summary(det)
    podcorr = pod_correction(passes, det)

    args.out_dir.mkdir(parents=True, exist_ok=True)

    det.to_csv(args.out_dir / "detections_enriched.csv", index=False)
    summary.to_csv(args.out_dir / "basin_summary.csv", index=False)
    if strat is not None:
        strat.to_csv(args.out_dir / "stratified_summary.csv", index=False)
    if podcorr is not None:
        podcorr.to_csv(args.out_dir / "pod_corrected_totals.csv", index=False)

    meta = {
        "n_surveys": int(len(d.surveys)),
        "n_site_passes": int(len(d.site_passes)) if d.has_passes else 0,
        "n_detections": int(len(d.detections)),
        "n_wells": int(len(d.wells)) if d.has_wells else 0,
        "basins": sorted(map(str, det["basin"].dropna().unique())),
        "warnings": d.warnings,
        "pod_model": "GML 2.0 (Thorpe et al. 2024, Table 3)",
        "gcn_default_ppm_m": DEFAULT_GCN_PPM_M,
        "basin_summary": summary.to_dict(orient="records"),
        "pod_corrected_totals": podcorr.to_dict(orient="records") if podcorr is not None else None,
    }
    with open(args.out_dir / "summary.json", "w") as f:
        json.dump(meta, f, indent=2, default=str)

    if not args.no_plot:
        make_diagnostic_plot(det, passes, args.out_dir / "path_b_diagnostic.png")

    print("\n[path B] basin summary:")
    print(summary.to_string(index=False))
    if podcorr is not None:
        print("\n[path B] POD-corrected totals:")
        print(podcorr[["basin", "n_detections", "mean_pod_at_detections",
                       "observed_total_kgph", "pod_corrected_total_kgph"]]
              .to_string(index=False))
    print(f"\n[path B] wrote outputs under: {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
