"""
Bridger LiDAR Survey Simulation - Realistic Approach
Models wind-dependent POD with random wind sampling from TMY data.

This script simulates a single Bridger survey of 20% of PA marginal wells,
accounting for seasonal wind variation effects on detection probability.

Author: FEAST Analysis
Date: 2026
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path
import pickle


# === CONSTANTS ===
POD_90_GPS = 1.27  # kg/h at which Bridger has 90% POD (Thorpe et al. 2024)
LOGISTIC_STEEPNESS = 2.0
BASELINE_WIND_MS = 3.5
WIND_EXPONENT = -0.3
MIN_WIND_MS = 1.0
MAX_WIND_MS = 6.0
SURVEY_COVERAGE_PCT = 0.20
SURVEY_DAYS_ASSUMED = 108  # ~120 wells/day → 12,925 wells / 120 = 108 days
N_ITERATIONS = 5


def pod_bridger(emission_rate_kgph, wind_speed_ms):
    """
    Calculate Bridger POD for given emission and wind speed.
    
    Based on Thorpe et al. (2024) Remote Sensing of Environment:
    - 90% POD at 1.27 kg/h (3.5 m/s wind)
    - Wind adjustment: exp(-0.3 * wind_delta)
    
    Args:
        emission_rate_kgph: Well emission in kg/h
        wind_speed_ms: Wind speed in m/s
    
    Returns:
        Detection probability [0, 1]
    """
    # Adjusted threshold based on wind
    wind_factor = np.exp(WIND_EXPONENT * (wind_speed_ms - BASELINE_WIND_MS))
    adjusted_threshold = POD_90_GPS * wind_factor
    
    # Logistic curve on log-scale
    if emission_rate_kgph <= 0:
        return 0.0
    
    log_e = np.log(emission_rate_kgph)
    log_threshold = np.log(adjusted_threshold)
    exponent = LOGISTIC_STEEPNESS * (log_e - log_threshold)
    pod = 1.0 / (1.0 + np.exp(-exponent))
    
    return np.clip(pod, 0.0, 1.0)


def sample_wind_from_tmy(tmy_path='../ExampleData/TMY-DataExample.csv'):
    """
    Sample wind speed from TMY (Typical Meteorological Year) data.
    Filtered to flyable conditions [1.0, 6.0] m/s.
    """
    if not Path(tmy_path).exists():
        print(f"Warning: TMY file not found at {tmy_path}. Using uniform [1-6] m/s distribution.")
        return np.random.uniform(MIN_WIND_MS, MAX_WIND_MS)
    
    try:
        tmy_df = pd.read_csv(tmy_path)
        if 'WindSpeed' in tmy_df.columns:
            wind_data = tmy_df['WindSpeed'].dropna().values
        elif 'wind_speed' in tmy_df.columns:
            wind_data = tmy_df['wind_speed'].dropna().values
        else:
            print("Warning: Could not find wind column in TMY. Columns:", tmy_df.columns.tolist())
            return np.random.uniform(MIN_WIND_MS, MAX_WIND_MS)
        
        # Filter to flyable range
        flyable = wind_data[(wind_data >= MIN_WIND_MS) & (wind_data <= MAX_WIND_MS)]
        if len(flyable) == 0:
            print("Warning: No flyable wind speeds in TMY. Using [1-6] m/s.")
            return np.random.uniform(MIN_WIND_MS, MAX_WIND_MS)
        
        return np.random.choice(flyable)
    except Exception as e:
        print(f"Error loading TMY: {e}. Using uniform [1-6] m/s.")
        return np.random.uniform(MIN_WIND_MS, MAX_WIND_MS)


def simulate_survey(emissions_df, coverage_pct=SURVEY_COVERAGE_PCT):
    """
    Simulate single Bridger survey with wind-dependent POD.
    
    Args:
        emissions_df: DataFrame with 'emission_rate_kgph' column
        coverage_pct: Fraction of wells to survey (default 20%)
    
    Returns:
        dict with detection results
    """
    n_wells = len(emissions_df)
    n_survey = int(n_wells * coverage_pct)
    
    # Random well selection
    survey_indices = np.random.choice(n_wells, size=n_survey, replace=False)
    survey_wells = emissions_df.iloc[survey_indices].copy()
    
    # Sample wind for this iteration
    wind_speed = sample_wind_from_tmy()
    
    # Calculate POD for each well
    survey_wells['pod'] = survey_wells['emission_rate_kgph'].apply(
        lambda e: pod_bridger(e, wind_speed)
    )
    
    # Bernoulli detection: P(detect) = POD, P(not detect) = 1-POD
    survey_wells['detected'] = np.random.binomial(1, survey_wells['pod'])
    
    # Results
    results = {
        'wind_speed_ms': wind_speed,
        'pods': survey_wells['pod'].values,
        'detected': survey_wells['detected'].values,
        'emissions_sampled': survey_wells['emission_rate_kgph'].sum(),
        'emissions_detected': (survey_wells['emission_rate_kgph'] * survey_wells['detected']).sum(),
        'wells_sampled': len(survey_wells),
        'wells_detected': int(survey_wells['detected'].sum()),
        'avg_pod': survey_wells['pod'].mean(),
    }
    
    return results


def main():
    """Main simulation."""
    # Load emissions data
    data_path = 'PAMarginalWellEmissions2023.pkl'
    if not Path(data_path).exists():
        print(f"Error: {data_path} not found")
        return
    
    print("Loading emissions data...")
    with open(data_path, 'rb') as f:
        gdf = pickle.load(f)
    
    if 'mean_ch4_kgh' in gdf.columns:
        emissions_col = 'mean_ch4_kgh'
    elif 'emission_rate_kgph' in gdf.columns:
        emissions_col = 'emission_rate_kgph'
    else:
        print(f"Error: Could not find emission column. Available: {gdf.columns.tolist()}")
        return
    
    emissions_df = pd.DataFrame({
        'emission_rate_kgph': gdf[emissions_col].values
    })
    
    total_emissions = emissions_df['emission_rate_kgph'].sum()
    
    print(f"\n{'='*80}")
    print("BRIDGER LiDAR SURVEY - PA MARGINAL WELLS (REALISTIC)")
    print(f"{'='*80}\n")
    
    print(f"Portfolio:")
    print(f"  Total wells: {len(emissions_df):,}")
    print(f"  Total emissions: {total_emissions:.1f} kg/h\n")
    
    # Run multiple iterations
    print(f"Running {N_ITERATIONS} survey iterations (20% coverage each)...\n")
    
    all_results = []
    
    for i in range(N_ITERATIONS):
        results = simulate_survey(emissions_df, SURVEY_COVERAGE_PCT)
        all_results.append(results)
        
        mitigation_pct = 100.0 * results['emissions_detected'] / total_emissions
        
        print(f"Iteration {i+1}:")
        print(f"  Wind speed: {results['wind_speed_ms']:.2f} m/s")
        print(f"  Wells surveyed: {results['wells_sampled']:,}")
        print(f"  Emissions surveyed: {results['emissions_sampled']:.1f} kg/h")
        print(f"  Average PoD: {results['avg_pod']:.3f}")
        print(f"  Wells detected: {results['wells_detected']:,} ({100*results['wells_detected']/results['wells_sampled']:.1f}%)")
        print(f"  Emissions detected: {results['emissions_detected']:.1f} kg/h")
        print(f"  Mitigation potential: {mitigation_pct:.1f}% of portfolio")
        print()
    
    # Summary statistics
    detected_emissions = np.array([r['emissions_detected'] for r in all_results])
    mean_detected = detected_emissions.mean()
    std_detected = detected_emissions.std()
    ci_lower = np.percentile(detected_emissions, 2.5)
    ci_upper = np.percentile(detected_emissions, 97.5)
    
    print(f"{'='*80}")
    print("SUMMARY (5 iterations)")
    print(f"{'='*80}")
    print(f"Mean emissions detected: {mean_detected:.0f} kg/h (±{std_detected:.0f} std)")
    print(f"Range: {detected_emissions.min():.0f} - {detected_emissions.max():.0f} kg/h")
    print(f"95% CI: {ci_lower:.0f} - {ci_upper:.0f} kg/h")
    print(f"Mean mitigation: {100*mean_detected/total_emissions:.1f}%")
    print(f"Range: {100*detected_emissions.min()/total_emissions:.1f}% - {100*detected_emissions.max()/total_emissions:.1f}%")
    
    # Save results
    output_dir = Path('BridgerResults')
    output_dir.mkdir(exist_ok=True)
    
    output_data = {
        'methodology': 'Single 20% survey with wind-dependent POD',
        'portfolio': {
            'total_wells': len(emissions_df),
            'total_emissions_kgph': float(total_emissions),
        },
        'survey_params': {
            'coverage_pct': SURVEY_COVERAGE_PCT,
            'iterations': N_ITERATIONS,
        },
        'iterations': all_results,
        'summary': {
            'mean_detected_kgph': float(mean_detected),
            'std_detected_kgph': float(std_detected),
            'ci_lower_kgph': float(ci_lower),
            'ci_upper_kgph': float(ci_upper),
            'mean_mitigation_pct': float(100*mean_detected/total_emissions),
        }
    }
    
    output_path = output_dir / 'bridger_survey_realistic.json'
    with open(output_path, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\nResults saved: {output_path}")


if __name__ == '__main__':
    main()
