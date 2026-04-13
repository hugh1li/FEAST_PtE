# Bridger LiDAR Survey of PA Marginal Wells: Grounded Analysis
### Addressing: Realistic Flight Plans, Emission Capture, and "Top Emitter vs. Most Emitter"

**Date:** April 2026  
**Sources:** Thorpe et al. 2024 (RSE 315:114435) | Donahue et al. 2024 (Permian preprint, ES&T submitted)  
**Portfolio:** Pennsylvania marginal wells — 64,624 wells, 55,389 kg/h total baseline

---

## The Core Question

> *If I do a Bridger survey of PA marginal wells, can I say that the survey captures **most emitters**, or does it simply capture the **top emitters** since top emitters contribute more than half of total emissions?*

**Short Answer:** Bridger captures the **top emitters**, not *most* emitters — and that's actually the right framing to use. Here's why, in detail.

---

## 1. The PA Marginal Well Emission Distribution

The portfolio is **highly right-skewed** (lognormal, typical for oil & gas):

| Top % of Wells | # Wells | Share of Total Emissions | Min Emission in Group |
|----------------|---------|--------------------------|----------------------|
| Top 1% | 646 | **37%** | 12.4 kg/h |
| Top 2% | 1,292 | **48%** | 7.2 kg/h |
| Top 5% | 3,231 | **64%** | 3.1 kg/h |
| Top 10% | 6,462 | **77%** | 1.5 kg/h |
| Top 20% | 12,924 | **88%** | 0.6 kg/h |
| Bottom 50% | 32,312 | **2%** | — |

The bottom **70% of wells** emit less than 0.32 kg/h each — Bridger's 10% PoD threshold. Those wells contribute just **7%** of total emissions. You cannot effectively find them with aerial surveys. This isn't a flaw of the survey — it's a fundamental property of the emission distribution.

---

## 2. The Bridger PoD Model (Marcellus-Specific)

**Critical correction from existing analysis:** The existing code used the **overall 2023 average** PoD₉₀ = 1.27 kg/h. Thorpe et al. (2024) Fig. 5 right panel gives **Marcellus-specific** performance: **mean PoD₉₀ = 0.974 kg/h**. Marcellus terrain (forested, vegetated) is characterized in the Bridger data with a corresponding GCN, and the sensitivity is modestly better than the overall average (0.974 vs 1.27 kg/h). This matters for PA analysis.

**PoD curve at 3.5 m/s wind (Marcellus-specific):**

| Emission Rate | PoD | # PA Wells at This Level | % of Fleet |
|--------------|-----|--------------------------|------------|
| 0.1 kg/h | 1% | 33,500 | 52% |
| 0.25 kg/h | 6% | 22,200 | 34% |
| 0.5 kg/h | 21% | 14,900 | 23% |
| **0.974 kg/h** | **50%** | **9,200** | **14%** |
| 1.27 kg/h | 63% | 7,400 | 11% |
| 2.0 kg/h | 81% | 5,000 | 8% |
| 5.0 kg/h | 96% | 2,000 | 3% |
| 10.0 kg/h | 99% | 830 | 1% |

**Key finding:** Only ~14% of wells (9,200 wells) emit at or above the 50% PoD level. But those 9,200 wells contribute ~83% of total emissions.

**Wind sensitivity** — how the 90% PoD threshold shifts:

| Wind Speed | 90% PoD Threshold | PoD at 1 kg/h |
|-----------|-------------------|----------------|
| 1.0 m/s (calm) | 2.06 kg/h | 19% |
| 2.0 m/s | 1.53 kg/h | 30% |
| **3.5 m/s (optimal)** | **0.97 kg/h** | **51%** |
| 5.0 m/s | 0.62 kg/h | 72% |
| 6.0 m/s | 0.46 kg/h | 83% |

Wind speed between 3.5–6.0 m/s is the sweet spot. PA summer afternoons typically deliver this. Calm mornings (<2 m/s) significantly degrade detection.

---

## 3. Realistic Flight Plan (Derived from Permian 2024)

The Permian preprint (Donahue et al. 2024, Table 1) provides the best benchmark for Bridger throughput in a real campaign:

**Permian 2024 Campaign Benchmark:**
- Total sites scanned: 51,770 across 195 calendar days (4 quarters)
- Marginal wellsites specifically: **6,357/quarter** over ~49 flying days/quarter
- → **~130 marginal wells/day** in the Permian

**PA Adjustments:**
- Permian: flat terrain, low vegetation, simple routing → PA has forested ridgelines, rural roads, and more airspace complexity → ~**-20% throughput adjustment**
- Permian: 195 flyable days/year (West Texas desert) → PA: more precipitation, cloud cover, winter icing → ~**150 flyable days/year** (conservative)

**Resulting PA flight plan:**

| Parameter | Permian Benchmark | PA Estimate |
|-----------|-------------------|-------------|
| Wells/day (marginal) | 130 | **100** |
| Flyable days/year | 195 | **150** |
| Wells/year | ~19,500 | **~15,000** |
| Annual coverage | — | **~23% of portfolio** |
| Days to survey full portfolio | — | **~644 days (~4.3 years)** |

**Quarterly allocation (PA):**
- Q1 Jan–Mar: ~38 days (cold/ice, fewest flyable)
- **Q2 Apr–Jun: ~40 days (best – mild winds, spring conditions)**
- **Q3 Jul–Sep: ~35 days (good – but thunderstorm interruptions)**
- Q4 Oct–Dec: ~37 days (declining flyability)

**Recommendation:** Prioritize Q2 and Q3 for maximum survey productivity. These 75 combined days cover ~7,500 wells in a single season.

---

## 4. Monte Carlo Survey Results (1,000 iterations)

| Scenario | Coverage | Mean Detected | % Mitigation | 90% CI |
|----------|----------|---------------|-------------|---------|
| One season (Q2+Q3, ~75 days) | 12% | 5,168 kg/h | 9.3% | [8.0%–11.2%] |
| **1-year moderate (150 days)** | **23%** | **10,366 kg/h** | **18.7%** | **[16.5%–21.1%]** |
| 1-year ambitious (Permian rate) | 30% | 13,539 kg/h | 24.4% | [21.8%–27.2%] |
| Full portfolio (one pass) | 100% | 44,674 kg/h | 80.7% | [74.8%–86.1%] |

**Why does a full single-pass not reach 100%?** Even with 100% coverage, wells below ~0.3 kg/h have near-zero PoD. Those ~45,000 wells contribute ~7% of emissions — they're largely undetectable with aerial flyover. You'd need ground-based OGI or component-level inspection to find them.

---

## 5. The "Most Emitters vs. Top Emitters" Question — Direct Answer

For a realistic **1-year, 23% coverage** survey:

**By count of emitters:**
- ~15,000 wells surveyed → ~2,465 wells detected (3.8% of all 64,624)
- **You cannot claim "most emitters captured."** You find the leaks in the 23% you fly over, not in the other 77%.

**By emission mass:**
- ~10,350 kg/h detected → 18.7% of portfolio emissions
- The ~2,465 detected wells emit an **average of 4.2 kg/h** — vs. the portfolio average of 0.857 kg/h
- **Enrichment factor: 4.9×** — detected wells are nearly 5× larger than average

**Top-emitter capture rate:**
- Top 5% emitters = 3,231 wells emitting ≥3.1 kg/h, representing **64% of total emissions**
- In a random 23% survey, you'll visit ~743 of them (23%)
- Their PoD is ~91%+ at 3.5 m/s wind → you detect ~675 of them
- That's ~21% of the top-5% emissions captured from just 1% of surveyed wells

**Correct claim you CAN make:**
> *"Our Bridger survey detected emissions that are disproportionately concentrated in high-emitting wells — the detected wells averaged 4.9× the portfolio mean emission. The top emitters in our surveyed sample were reliably detected (>90% PoD), making this an effective tool for identifying the largest contributors to total emissions."*

**Claim you CANNOT make:**
> ~~"We captured most emitters"~~ — This would imply you found most of the individual leaking wells, which requires surveying much more of the portfolio.

---

## 6. Detection Blind Spot

| Below threshold | # Wells | % of Fleet | Emissions | % of Total |
|----------------|---------|-----------|-----------|-----------|
| < 0.32 kg/h (PoD < 10%) | 45,380 | **70%** | 3,734 kg/h | **7%** |
| < 0.1 kg/h | ~31,000 | ~48% | ~1,100 kg/h | ~2% |

**The good news:** The majority of undetectable wells (those emitting <0.32 kg/h) collectively contribute only ~7% of total emissions. Bridger misses the individual count of small emitters but doesn't miss much of the actual methane.

This is why the correct framing is: **Bridger is emission-mass efficient, not emitter-count efficient.**

---

## 7. Comparison of Survey Claim Accuracy

| Claim | Accurate? | Notes |
|-------|-----------|-------|
| "We surveyed X% of wells" | ✅ Yes | Straightforward |
| "We detected the top emitters" | ✅ Yes, *within surveyed wells* | PoD >90% for large emitters; detected wells 4.9× average size |
| "Our detected emissions dominate the survey's emission mass" | ✅ Yes | Top emitters account for most of what Bridger flags |
| "We captured most emitters in the portfolio" | ❌ No | Only 3-4% of wells detected in 1-year survey |
| "We captured the wells responsible for most emissions" | ⚠️ Partial | You captured a proportional random sample; top emitters in that sample are reliably found |
| "A Bridger survey finds the worst leakers" | ✅ Yes | This is the strongest legitimate claim |

---

## 8. Practical Implications

**If your goal is to identify and fix the largest sources:** Bridger is extremely well-suited. Within any surveyed subset, it will find the high emitters reliably. The limiting factor is coverage — you need to fly over the right wells to find their leaks.

**Recommended strategy for maximum emissions impact:**
1. **Year 1 (23% random survey):** Establish a baseline, identify ~675 top-5% emitters, repair them → ~7,400 kg/h removed (~13% of portfolio)
2. **Year 2 (revisit repaired + survey new 23%):** Confirm repairs, scan new wells → cumulative ~25-30% reduction
3. **Years 3–5 (rotating coverage):** By year 4–5, approach full portfolio coverage with ~80% emissions identified

**If the PA data is available**, the analysis can be sharpened considerably by:
- Using actual lat/lon for spatial clustering (target well-pad density areas first)
- Stratifying by county or basin area to optimize routing
- Prioritizing wells with known production histories (higher production → higher emission probability)

---

## 9. Code Files

| File | Purpose |
|------|---------|
| `bridger_pa_analysis.py` | **NEW** — Grounded analysis with Marcellus PoD, Permian flight plan, 1000-iteration MC |
| `bridger_survey_realistic.py` | Original — Wind-dependent single survey (uses 1.27 kg/h overall threshold) |
| `bridger_spatial_survey.py` | Original — DBSCAN spatial clustering (has a Q3 outlier bug due to DBSCAN eps) |
| `survey_strategy_summary.py` | Original — Strategy comparison |

**Key fix in new analysis:** Uses Marcellus-specific PoD₉₀ = 0.974 kg/h (not 1.27 overall avg), and flight throughput grounded in Permian paper Table 1 (~130/day Permian → ~100/day PA-adjusted).

---

## Summary

A Bridger survey of PA marginal wells captures **top emitters**, not *most* emitters:

- A realistic 1-year survey (~23% portfolio coverage) detects ~10,350 kg/h — **18.7% mitigation**
- Detected wells emit on average **4.9× the portfolio mean** — they are the large sources
- **70% of wells** are essentially invisible to Bridger (emit <0.32 kg/h) — but they only contribute **7% of emissions**
- The correct claim is: *"Bridger identifies the highest-emitting wells in surveyed areas, which are responsible for the majority of emissions within those areas"*
- To capture *most* of the **top-emitting wells** (those in the top 5% by emission), you would need to survey ~80–90% of the portfolio or implement stratified targeted sampling

The Permian 2024 campaign (Donahue et al.) provides the clearest precedent: they achieved 29% sampling fraction for marginal wellsites per quarter, explicitly targeting them across strata. A PA campaign following this design would achieve similar coverage rates.

---

*Analysis grounded in:*  
- *Thorpe et al. 2024, Remote Sensing of Environment 315:114435 — GML 2.0 PoD model & Marcellus basin detection sensitivity*  
- *Donahue et al. 2024, preprint submitted to ES&T — Permian Basin aerial campaign: 51,770 sites, 130 marginal wells/flying day*
