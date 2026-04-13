# Bridger LiDAR Survey — PA Marginal Wells: Full Real-Data Analysis
### Using actual GeoPackage data (64,624 wells, 250 FEAST stochastic iterations, real lat/lon)

**Date:** April 2026  
**Code:** `bridger_pa_full_analysis.py`  
**Data:** `PAMarginalWellEmissions2023.gpkg` + `_withEachIteration.gpkg`  
**Sources:** Thorpe et al. 2024 (RSE 315:114435) · Donahue et al. 2024 (Permian preprint)

---

## Critical Methodological Corrections from Prior Analysis

Three important fixes were applied in this version:

**1. Correct denominator.** Mitigation % now = *detected emissions / total emissions of the surveyed wells*, not / total portfolio. This measures within-survey detection efficiency — how much of what was actually flown over did Bridger find. The old approach (/ total portfolio) conflated survey coverage with detection performance.

**2. Use individual FEAST iterations, not time-averaged means.** Each Monte Carlo run draws one random FEAST iteration (0–249) as the emission snapshot at time of survey. This captures the stochastic, intermittent nature of emissions — some wells will be in high-emission states, others near-zero. Using the iteration mean flattens out super-emitters and underestimates detection efficiency.

**3. Spatial grid-cell survey routing.** Wells are selected as 0.1°×0.1° geographic grid cells (~11×8.5 km), mimicking how an aircraft surveys contiguous areas rather than random scattered individual wells.

---

## 1. The Emission Data: Mean vs. Stochastic Snapshot

The FEAST model provides 250 independent stochastic realizations per well:

| Quantity | Mean emissions (time-avg) | Single FEAST iteration (snapshot) |
|----------|--------------------------|-----------------------------------|
| Total portfolio | 55,389 kg/h | 33,844–89,207 kg/h (±19%) |
| Max single well | 2.96 kg/h | up to ~80 kg/h (super-emitters) |
| % wells > 0.974 kg/h | 40% | 7.3% ± 6.7% per snapshot |
| Median emission | 0.50 kg/h | ~0.45 kg/h |

**This is the most important structural point.** The time-averaged data makes PA marginal wells look like they have a smooth, moderate distribution peaking at ~0.44–0.50 kg/h with a long tail to ~2.96 kg/h (the FEAST maximum). But any *given* survey day resembles the stochastic snapshot: most wells are at near-baseline (~0.4–0.6 kg/h) and a small fraction are in high-emission states (>2 kg/h, up to ~80 kg/h for upset events). These snapshot super-emitters are exactly what Bridger excels at finding.

---

## 2. Emission Concentration (time-averaged mean)

| Top % of Wells | # Wells | Share of Emissions | Min emission |
|----------------|---------|-------------------|-------------|
| Top 1% | 646 | 2.5% | 1.97 kg/h |
| Top 5% | 3,231 | 10.9% | 1.72 kg/h |
| Top 10% | 6,462 | 20.6% | 1.59 kg/h |
| **Top 30%** | **19,387** | **54.1%** | **1.30 kg/h** |
| Top 50% | 32,312 | 74.3% | 0.50 kg/h |

Using the time-averaged mean, the concentration is *much less extreme* than in the Permian or typical O&G super-emitter studies. The top 30% of wells (emitting ≥1.30 kg/h on average) account for >50% of emissions. This is because the PA data is driven by FEAST's component-level emission factors for conventional vertical wells — a relatively homogeneous fleet.

**Using individual iterations (more realistic for a survey):** the top 5% of wells per snapshot contribute ~40% ± 8% of that iteration's total emissions, and those top-5% wells are changing composition between snapshots (intermittent super-emitters dominate in any given pass).

---

## 3. Bridger Detection Against This Distribution

At the **Marcellus-specific PoD₉₀ = 0.974 kg/h** (Thorpe 2024, Fig. 5):

| Threshold | % of fleet (time-avg) | % of emissions | Per-snapshot average |
|-----------|----------------------|----------------|---------------------|
| > 0.5 kg/h | 47.9% | 73.0% | varies by iteration |
| **> 0.974 kg/h** | **40.0%** | **68.1%** | **~7.3% per snapshot** |
| > 1.27 kg/h | 32.1% | 57.3% | ~5.5% per snapshot |

The time-averaged view shows 40% of wells above the PoD₉₀ threshold. But on any given survey day, only ~7.3% of wells are in a high-emission state. The key insight is that **Bridger will reliably catch whichever wells happen to be emitting heavily on survey day** — and those wells carry ~40% of that day's total emissions (top-5% share ~40% per iteration).

---

## 4. Spatial Survey Plan (real lat/lon)

Using actual well coordinates, a 0.1°×0.1° (~11×8.5 km) grid produces **474 occupied survey cells**:

| Metric | Value |
|--------|-------|
| Survey cells | 474 |
| Mean wells/cell | 136 |
| Median wells/cell | 53 |
| Largest cell (Indiana Co.) | 1,799 wells |
| Days to survey all cells @ 100/day | **646 days** |
| Years for full portfolio pass | **4.3 years** |

**Annual coverage at 100 wells/day, 150 flyable days:**

| Quarter | Days | Wells | Notes |
|---------|------|-------|-------|
| Q1 Jan–Mar | 38 | 3,800 | Winter — fewest flyable |
| **Q2 Apr–Jun** | **40** | **4,000** | **Recommended — mild, consistent winds** |
| **Q3 Jul–Sep** | **35** | **3,500** | **Good — moderate thunderstorm gaps** |
| Q4 Oct–Dec | 37 | 3,700 | Declining flyability |
| **Annual total** | **150** | **15,000** | **23.2% of portfolio** |

The throughput benchmark is directly from Donahue et al. 2024: Bridger scanned 6,357 marginal wellsites per quarter in the Permian across ~49 flying days = 130 sites/day. A 23% terrain/weather adjustment gives ~100/day for PA.

---

## 5. Monte Carlo Results (500 runs, real FEAST iterations)

**Key change from prior analysis:** denominator = surveyed-wells' total emissions per iteration.

| Scenario | Coverage | Within-Survey Mitigation % | Detected (kg/h) | Top-5% Capture | Enrichment |
|----------|----------|---------------------------|----------------|---------------|-----------|
| One season (Q2+Q3, ~75 days) | 12% | **56% [39–71%]** | 3,653 ± 1,222 | 95% | 2.3× |
| **1-year moderate (150 days)** | **23%** | **56.5% [39–72%]** | **7,353 ± 2,446** | **95%** | **2.3×** |
| 1-year ambitious (130/day) | 30% | 57% [38–71%] | 9,637 ± 3,402 | 95% | 2.3× |
| Full portfolio pass | 100% | 57% [39–72%] | 32,330 ± 10,321 | 95% | 2.3× |

**What the 56% within-survey mitigation means:** In each surveyed 0.1°×0.1° cell, Bridger detects on average **56% of that cell's estimated emissions** on survey day. The remaining 44% is either below the detection threshold (low-emission wells, PoD ~1–20%) or missed due to wind conditions.

**The wide 90% CI ([39–72%])** comes almost entirely from FEAST emission variability (some iterations have most wells at baseline, others have many super-emitters) combined with wind variability. When many wells are in high-emission states on survey day, detection is high; when most are at baseline (~0.45 kg/h, PoD ~20%), detection is moderate.

**Detection enrichment 2.3×:** Detected wells emit on average 2.3× the mean of their surveyed-zone peers — Bridger disproportionately flags the highest emitters within each surveyed area.

---

## 6. Answer to the Core Question

> *Does a Bridger survey capture 'most emitters' or 'top emitters since top emitters contribute more than half of emissions'?*

**The answer is: TOP EMITTERS — and with very high efficiency within surveyed zones.**

| Claim | Verdict | Evidence |
|-------|---------|----------|
| "We captured most emitters (by well count)" | ❌ No | Only 6.3% of all 64,624 wells detected in a 1-year survey |
| "We captured most emissions of the portfolio" | ❌ No | 23% coverage × 56% within-survey detection ≈ 13% of total portfolio emissions |
| "We captured the top emitters in areas we surveyed" | ✅ Yes | Top-5% emitters in surveyed zones: 95% capture; enrichment 2.3× |
| "Detected emissions are dominated by top emitters" | ✅ Yes | Detected wells emit 2.3× average; top-5% of each surveyed zone nearly all found |
| "Top emitters account for >50% of emissions" | ✅ Yes (but at p30, not p5) | Top 30% of wells (≥1.30 kg/h mean) = 54% of emissions |

**Nuanced but correct framing:**
> *"A Bridger survey reliably detects the highest-emission wells within each surveyed area, capturing approximately 56% of that zone's estimated emissions. Since emission concentration is moderate (top 30% of wells = 54% of portfolio emissions), Bridger is well-suited to finding the predominant sources when it flies over them. The primary limitation is coverage: a 1-year campaign covers ~23% of the portfolio, so 77% of wells — including their top emitters — remain unvisited."*

**On the "top emitters > half of total" argument:** With this specific PA dataset (FEAST conventional vertical wells), the emission distribution is less extreme than in Permian shale or production-weighted inventories. Top 5% of wells = only 11% of emissions (time-averaged). Top 30% = 54%. So the "top emitters = majority of emissions" threshold is at the **30th percentile**, not the 5th. Using individual FEAST iterations, top 5% by snapshot emission = ~40% of that day's emissions (more concentrated because of intermittent super-emitters). The claim that "top emitters contribute more than half" is most defensible when applied to the stochastic snapshot view.

---

## 7. Key Numbers for Reporting

For a **1-year, 23% coverage, spatial-routing** Bridger survey of PA marginal wells:

| Metric | Value |
|--------|-------|
| Wells surveyed/year | ~15,000 |
| Wells detected/year | ~4,050 (6.3% of total fleet) |
| Within-survey mitigation | **56% ± 10%** [39–72% 90% CI] |
| Emissions detected/year | ~7,350 ± 2,450 kg/h |
| As % of total portfolio | ~13% |
| Top-5% emitter capture (surveyed zones) | **~95%** |
| Detection enrichment vs average | **2.3×** |
| Full portfolio pass duration | **4.3 years** |

The large variability (±2,450 kg/h std) is dominated by FEAST emission stochasticity — whether super-emitter events occur in surveyed wells on survey day. This is not reducible by better survey design; it reflects genuine temporal variability in the emission field.

---

## 8. Outputs

| File | Description |
|------|-------------|
| `bridger_pa_full_analysis.py` | Full analysis code (real data, spatial survey, FEAST iterations) |
| `BridgerResults/bridger_pa_full_analysis.png` | 6-panel figure: distribution, Lorenz, PoD, mitigation, detection, spatial map |
| `BridgerResults/feast_iteration_variability.png` | Mean vs. snapshot emission distribution + iteration variability |
| `BridgerResults/bridger_pa_full_analysis.json` | Machine-readable results |
| `BRIDGER_PA_FULL_ANALYSIS.md` | This document |

---

*Thorpe et al. 2024, RSE 315:114435 — Bridger GML 2.0 PoD; Marcellus basin PoD₉₀ = 0.974 kg/h*  
*Donahue et al. 2024, ES&T preprint — Permian campaign: 6,357 marginal wellsites/quarter over ~49 days = 130/day*
