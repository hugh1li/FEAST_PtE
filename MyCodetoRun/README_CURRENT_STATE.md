# MyCodetoRun Folder — Current State (Mar 31, 2026)

## 🎯 Active Analysis Script
- **`run_bridger_analysis.py`** — Main Bridger LiDAR analysis pipeline
  - **Status:** ✅ Current (updated Mar 31, 2026)
  - **Approach:** FEAST-backed TMY meteorology + spatial wind field across surveyed wells
  - **Output:** 500 MC iterations of 20% surveys on PA marginal wells (64,624 total)
  - **Latest Result:** ~15.4% portfolio mitigation per survey (95% CI: 11.6%–19.5%)

## 📊 Current Output Files
All in `BridgerResults/` with timestamps Mar 31, 2026:
- **`bridger_analysis.png`** — 4-panel diagnostic plot (distributions, POD curves, detection rates)
- **`mitigation_driver_concentration.png`** — 2-panel concentration analysis (well count vs emission mass)
- **`simulation_results.csv`** — Raw MC iteration data (wind speeds, detection counts, POD)
- **`mitigation_driver_concentration_table.csv`** — Emission decile binning results

## 📁 Archive Folders
**`OldScripts/`** — Superseded implementations:
- `bridger_spatial_survey.py` (has DINWD_N=1.0 bug)
- `bridger_survey_realistic.py` (simple uniform-wind approach)
- `survey_strategy_summary.py` (utility script)

**`BridgerResults/OldRuns/`** — Obsolete outputs:
- `bridger_spatial_survey.json` (Feb 10)
- `bridger_survey_realistic.json` (Mar 18, large 2.4 MB)
- `strategy_comparison.json` (Feb 10)

**`Archive/`** — Duplicate/auxiliary data:
- `feast_emissions.pkl` (redundant pickle format of CSV)
- `spike0spiegel-methaneair-ea316b306af7.json` (unused metadata)

## 📋 Data & Documentation
- **`feast_emissions.csv`** (2.5 MB) — Current emission data (64,624 PA wells)
- **`ingest_bridger.py`** — Script for ingesting well emission data from raw sources
- **`BRIDGER_SURVEY_ANALYSIS.md`** — Spatial clustering analysis & strategy notes
- **`SPATIAL_VS_NONSPATIAL_ANALYSIS.md`** — Comparison of approaches
- **`pan_permian_manuscript_eartharxiv_v1.pdf`** (2.9 MB) — Reference paper on Bridger LiDAR
- **`PAMarginalWellEmissions2023.zip`** (119 MB) — Source emission archive (keep for reproducibility)

## 🚀 How to Run Current Analysis
```bash
cd /workspaces/FEAST_PtE
source .venv/bin/activate
python MyCodetoRun/run_bridger_analysis.py
```

## ✍️ Key Updates (Mar 31, 2026)
1. **FEAST Integration:** Script now builds a minimal FEAST GasField to load TMY meteorology
2. **Spatial Wind Field:** Each surveyed well's wind is sampled from a coordinate-based field, not one survey-wide value
3. **Output Clarity:** Wind statistics now reported (base, mean, min, max, std across surveyed wells)
4. **Documentation:** Script docstring and output summary clarify the FEAST-backed but analysis-specific approach

---
**Last updated:** Mar 31, 2026 | **Status:** Ready to use
