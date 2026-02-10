"""
Bridger LiDAR Survey - Spatial Clustering Approach
Uses DBSCAN to cluster wells by geographic proximity, then plans surveys
based on Permian 2024 quarterly schedule.

Author: FEAST Analysis
Date: 2026
"""

import pandas as pd
import numpy as np
import pickle
import json
from pathlib import Path
from sklearn.cluster import DBSCAN
from scipy.spatial.distance import cdist


# === CONSTANTS ===
DBSCAN_EPS_KM = 0.8  # Neighborhood radius (facility pad scale)
DBSCAN_MIN_SAMPLES = 2
POD_90_KGPH = 1.27
LOGISTIC_STEEPNESS = 2.0
BASELINE_WIND_MS = 3.5
WIND_EXPONENT = -0.3
MIN_WIND_MS = 1.0
MAX_WIND_MS = 6.0

# Permian 2024 schedule parameters
PERMIAN_QUARTERLY_COVERAGE = 0.20  # 20% of clusters per quarter
PERMIAN_DAYS_PER_YEAR = 195
PERMIAN_SENSORS = 2.7
FLIGHT_SPEED_KMH = 70
SURVEY_TIME_BASELINE_MIN = 2  # min per cluster
SURVEY_TIME_PER_WELL_MIN = 0.5  # min per well


def haversine_distance(lat1, lon1, lat2, lon2):
    """Calculate distance between two lat/lon points in km."""
    R = 6371  # Earth radius in km
    phi1, phi2 = np.radians(lat1), np.radians(lat2)
    delta_phi = np.radians(lat2 - lat1)
    delta_lambda = np.radians(lon2 - lon1)
    
    a = np.sin(delta_phi/2)**2 + np.cos(phi1)*np.cos(phi2)*np.sin(delta_lambda/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    return R * c


def cluster_wells_dbscan(gdf, eps_km=DBSCAN_EPS_KM, min_samples=DBSCAN_MIN_SAMPLES):
    """
    Cluster wells using DBSCAN on geographic coordinates.
    
    Args:
        gdf: GeoDataFrame with Latitude, Longitude
        eps_km: Search radius in km
        min_samples: Minimum cluster size
    
    Returns:
        cluster_labels: Array of cluster IDs (-1 for noise)
        centroids: Dict mapping cluster_id -> (lat, lon)
    """
    coords = gdf[['Latitude', 'Longitude']].values
    
    # Convert km to degrees (rough: 1 degree ≈ 111 km at equator)
    eps_degrees = eps_km / 111.0
    
    # Run DBSCAN
    clustering = DBSCAN(eps=eps_degrees, min_samples=min_samples).fit(coords)
    labels = clustering.labels_
    
    # Compute centroids
    centroids = {}
    for cluster_id in set(labels):
        if cluster_id == -1:  # Skip noise
            continue
        mask = labels == cluster_id
        centroid_lat = coords[mask, 0].mean()
        centroid_lon = coords[mask, 1].mean()
        centroids[cluster_id] = (centroid_lat, centroid_lon)
    
    return labels, centroids


def estimate_survey_time_per_cluster(gdf, cluster_labels, cluster_id):
    """Estimate survey time for a cluster in minutes."""
    mask = cluster_labels == cluster_id
    n_wells = mask.sum()
    survey_time = SURVEY_TIME_BASELINE_MIN + n_wells * SURVEY_TIME_PER_WELL_MIN
    return survey_time


def estimate_travel_time(lat1, lon1, lat2, lon2):
    """Estimate travel time between clusters in minutes."""
    distance_km = haversine_distance(lat1, lon1, lat2, lon2)
    travel_time_hours = distance_km / FLIGHT_SPEED_KMH
    return travel_time_hours * 60


def pod_bridger(emission_rate_kgph, wind_speed_ms):
    """Calculate Bridger POD."""
    wind_factor = np.exp(WIND_EXPONENT * (wind_speed_ms - BASELINE_WIND_MS))
    adjusted_threshold = POD_90_KGPH * wind_factor
    
    if emission_rate_kgph <= 0:
        return 0.0
    
    log_e = np.log(emission_rate_kgph)
    log_threshold = np.log(adjusted_threshold)
    exponent = LOGISTIC_STEEPNESS * (log_e - log_threshold)
    pod = 1.0 / (1.0 + np.exp(-exponent))
    
    return np.clip(pod, 0.0, 1.0)


def sample_wind_from_tmy(tmy_path='../ExampleData/TMY-DataExample.csv'):
    """Sample wind speed from TMY data."""
    if not Path(tmy_path).exists():
        return np.random.uniform(MIN_WIND_MS, MAX_WIND_MS)
    
    try:
        tmy_df = pd.read_csv(tmy_path)
        if 'WindSpeed' in tmy_df.columns:
            wind_data = tmy_df['WindSpeed'].dropna().values
        elif 'wind_speed' in tmy_df.columns:
            wind_data = tmy_df['wind_speed'].dropna().values
        else:
            return np.random.uniform(MIN_WIND_MS, MAX_WIND_MS)
        
        flyable = wind_data[(wind_data >= MIN_WIND_MS) & (wind_data <= MAX_WIND_MS)]
        if len(flyable) == 0:
            return np.random.uniform(MIN_WIND_MS, MAX_WIND_MS)
        
        return np.random.choice(flyable)
    except:
        return np.random.uniform(MIN_WIND_MS, MAX_WIND_MS)


def simulate_quarterly_survey(gdf, cluster_labels, cluster_id_list, 
                              emission_col='mean_ch4_kgh'):
    """
    Simulate survey of selected clusters for one quarter.
    
    Returns:
        dict with detection results
    """
    total_wells = 0
    total_emissions_sampled = 0.0
    detected_wells = 0
    detected_emissions = 0.0
    pod_list = []
    
    wind_speed = sample_wind_from_tmy()
    
    for cluster_id in cluster_id_list:
        mask = cluster_labels == cluster_id
        cluster_wells = gdf[mask].copy()
        
        # Calculate POD for each well
        cluster_wells['pod'] = cluster_wells[emission_col].apply(
            lambda e: pod_bridger(e, wind_speed)
        )
        
        # Bernoulli detection
        cluster_wells['detected'] = np.random.binomial(1, cluster_wells['pod'])
        
        total_wells += len(cluster_wells)
        total_emissions_sampled += cluster_wells[emission_col].sum()
        detected_wells += int(cluster_wells['detected'].sum())
        detected_emissions += (cluster_wells[emission_col] * cluster_wells['detected']).sum()
        pod_list.extend(cluster_wells['pod'].values)
    
    result = {
        'clusters_surveyed': len(cluster_id_list),
        'wells_surveyed': total_wells,
        'emissions_surveyed_kgph': total_emissions_sampled,
        'wind_speed_ms': wind_speed,
        'avg_pod': np.mean(pod_list) if pod_list else 0.0,
        'wells_detected': detected_wells,
        'emissions_detected_kgph': detected_emissions,
    }
    
    return result


def main():
    """Main spatial survey simulation."""
    # Load data
    data_path = 'PAMarginalWellEmissions2023.pkl'
    if not Path(data_path).exists():
        print(f"Error: {data_path} not found")
        return
    
    print("Loading emissions data...")
    with open(data_path, 'rb') as f:
        gdf = pickle.load(f)
    
    # Detect emission column
    if 'mean_ch4_kgh' in gdf.columns:
        emission_col = 'mean_ch4_kgh'
    elif 'emission_rate_kgph' in gdf.columns:
        emission_col = 'emission_rate_kgph'
    else:
        print(f"Error: Could not find emission column")
        return
    
    total_emissions = gdf[emission_col].sum()
    
    print(f"\n{'='*80}")
    print("BRIDGER LIDA LiDAR SURVEY - PA MARGINAL WELLS (SPATIAL CLUSTERING)")
    print(f"{'='*80}\n")
    
    print("1. PORTFOLIO BASELINE:")
    print(f"   Total wells: {len(gdf):,}")
    print(f"   Total emissions: {total_emissions:.1f} kg/h\n")
    
    # Cluster wells
    print("2. SPATIAL CLUSTERING (DBSCAN eps=0.8 km):")
    cluster_labels, centroids = cluster_wells_dbscan(gdf, DBSCAN_EPS_KM, DBSCAN_MIN_SAMPLES)
    n_clusters = len(centroids)
    n_noise = (cluster_labels == -1).sum()
    avg_wells_per_cluster = (len(gdf) - n_noise) / n_clusters if n_clusters > 0 else 0
    
    print(f"   Clusters identified: {n_clusters}")
    print(f"   Noise wells (isolated): {n_noise}")
    print(f"   Avg wells per cluster: {avg_wells_per_cluster:.1f}\n")
    
    # Survey time estimation
    print("3. SURVEY METRICS:")
    total_survey_time = 0
    for cluster_id in range(len(centroids)):
        survey_time = estimate_survey_time_per_cluster(gdf, cluster_labels, cluster_id)
        total_survey_time += survey_time
    
    avg_survey_time = total_survey_time / n_clusters if n_clusters > 0 else 0
    
    print(f"   Total survey time (clusters only): {total_survey_time/60:.1f} hours")
    print(f"   Avg survey time per cluster: {avg_survey_time:.1f} min\n")
    
    # Permian schedule
    print("4. PERMIAN-BASED SCHEDULE:")
    quarterly_clusters = int(n_clusters * PERMIAN_QUARTERLY_COVERAGE)
    time_per_quarter = quarterly_clusters * avg_survey_time / 60  # hours
    days_per_quarter = time_per_quarter / 8  # 8-hour flight days
    timeline_years = (n_clusters * PERMIAN_QUARTERLY_COVERAGE) / (quarterly_clusters / 4) if quarterly_clusters > 0 else np.inf
    
    print(f"   Survey coverage: {PERMIAN_QUARTERLY_COVERAGE*100:.0f}% of clusters (~{quarterly_clusters} clusters)")
    print(f"   Survey time per quarter: {days_per_quarter:.1f} days")
    print(f"   Days available per quarter: {PERMIAN_DAYS_PER_YEAR/4:.0f} days")
    print(f"   Timeline to 100% coverage: {timeline_years:.1f} years\n")
    
    # Quarterly surveys
    print("5. QUARTERLY SURVEY PLAN:")
    print(f"   Q1 (Jan-Mar): Survey 20% of clusters = {quarterly_clusters} clusters")
    print(f"   Q2 (Apr-Jun): Survey 20% of clusters = {quarterly_clusters} clusters")
    print(f"   Q3 (Jul-Sep): Survey 20% of clusters = {quarterly_clusters} clusters")
    print(f"   Q4 (Oct-Dec): Survey 20% of clusters = {quarterly_clusters} clusters\n")
    
    # Wind conditions
    print("6. WIND CONDITIONS (Flyable: 1.0-6.0 m/s):")
    try:
        tmy_df = pd.read_csv('../ExampleData/TMY-DataExample.csv')
        if 'WindSpeed' in tmy_df.columns:
            wind_data = tmy_df['WindSpeed'].dropna().values
        else:
            wind_data = tmy_df['wind_speed'].dropna().values
        
        flyable = wind_data[(wind_data >= MIN_WIND_MS) & (wind_data <= MAX_WIND_MS)]
        flyable_days = len(flyable)
        flyable_pct = 100 * len(flyable) / len(wind_data)
        
        print(f"   Mean: {wind_data.mean():.2f} m/s")
        print(f"   Days flyable: ~{flyable_days} days/year ({flyable_pct:.0f}%)\n")
    except:
        print(f"   Mean: 3.32 m/s (estimated)")
        print(f"   Days flyable: ~219 days/year\n")
    
    # Run quarterly simulations
    print("7. QUARTERLY DETECTION RESULTS:")
    print(f"{'='*80}\n")
    
    all_results = []
    
    for quarter in range(4):
        quarter_names = ['Q1 (Jan-Mar)', 'Q2 (Apr-Jun)', 'Q3 (Jul-Sep)', 'Q4 (Oct-Dec)']
        
        # Random cluster selection
        cluster_ids = list(centroids.keys())
        selected_clusters = np.random.choice(cluster_ids, size=quarterly_clusters, replace=False)
        
        # Simulate survey
        result = simulate_quarterly_survey(gdf, cluster_labels, selected_clusters, emission_col)
        all_results.append(result)
        
        mitigation_pct = 100.0 * result['emissions_detected_kgph'] / total_emissions
        
        print(f"{quarter_names[quarter]}:")
        print(f"  Clusters surveyed: {result['clusters_surveyed']}")
        print(f"  Wells surveyed: {result['wells_surveyed']:,}")
        print(f"  Emissions surveyed: {result['emissions_surveyed_kgph']:.1f} kg/h")
        print(f"  Wind speed: {result['wind_speed_ms']:.2f} m/s")
        print(f"  Avg PoD: {result['avg_pod']:.3f}")
        print(f"  Wells detected: {result['wells_detected']:,} ({100*result['wells_detected']/max(1,result['wells_surveyed']):.1f}%)")
        print(f"  Emissions detected: {result['emissions_detected_kgph']:.1f} kg/h")
        print(f"  Mitigation potential: {mitigation_pct:.2f}% of portfolio")
        print()
    
    print(f"{'='*80}")
    print("ANNUAL SUMMARY:")
    total_detected = sum(r['emissions_detected_kgph'] for r in all_results)
    annual_mitigation = 100.0 * total_detected / total_emissions
    avg_quarterly = total_detected / 4
    
    print(f"  Total detected across 4 quarters: {total_detected:.1f} kg/h")
    print(f"  Portfolio reduction potential: {annual_mitigation:.1f}%")
    print(f"  Average per quarter: {avg_quarterly:.1f} kg/h")
    
    # Save results
    output_dir = Path('BridgerResults')
    output_dir.mkdir(exist_ok=True)
    
    output_data = {
        'methodology': 'Spatial clustering (DBSCAN) + Permian 2024 schedule',
        'clustering': {
            'eps_km': DBSCAN_EPS_KM,
            'n_clusters': n_clusters,
            'n_noise_wells': int(n_noise),
        },
        'portfolio': {
            'total_wells': len(gdf),
            'total_emissions_kgph': float(total_emissions),
        },
        'survey_schedule': {
            'model': 'Permian 2024: 4 quarterly surveys, 195 days/year, 2.7 sensors',
            'quarterly_coverage_pct': PERMIAN_QUARTERLY_COVERAGE,
            'days_per_quarter': float(days_per_quarter),
        },
        'annual_results': all_results,
        'summary': {
            'total_detected_kgph': float(total_detected),
            'total_mitigation_pct': float(annual_mitigation),
            'avg_detected_per_quarter_kgph': float(avg_quarterly),
        }
    }
    
    output_path = output_dir / 'bridger_spatial_survey.json'
    with open(output_path, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\nResults saved: {output_path}")


if __name__ == '__main__':
    main()
