"""
Ingest Bridger LiDAR survey data from ZIP file (GeoPackage or pickle format)
into FEAST-compatible CSV and pickle formats.

Author: FEAST Analysis
Date: 2026
"""

import pandas as pd
import geopandas as gpd
import pickle
import os
from pathlib import Path
import zipfile
import json


def extract_known_from_zip(zip_path):
    """
    Extract data from ZIP containing GeoPackage (.gpkg) or pickle (.pkl) files.
    Returns GeoDataFrame with emissions and coordinates.
    """
    extracted_files = {}
    
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        for file_info in zip_ref.filelist:
            if file_info.filename.endswith('.gpkg'):
                print(f"Found GeoPackage: {file_info.filename}")
                with zip_ref.open(file_info.filename) as f:
                    temp_path = f'/tmp/{file_info.filename}'
                    with open(temp_path, 'wb') as tmp:
                        tmp.write(f.read())
                    extracted_files['gpkg'] = temp_path
            elif file_info.filename.endswith('.pkl'):
                print(f"Found pickle: {file_info.filename}")
                with zip_ref.open(file_info.filename) as f:
                    temp_path = f'/tmp/{file_info.filename}'
                    with open(temp_path, 'wb') as tmp:
                        tmp.write(f.read())
                    extracted_files['pkl'] = temp_path
    
    return extracted_files


def load_data(file_path):
    """Load data from .gpkg or .pkl file."""
    if file_path.endswith('.gpkg'):
        return gpd.read_file(file_path)
    elif file_path.endswith('.pkl'):
        with open(file_path, 'rb') as f:
            return pickle.load(f)
    else:
        raise ValueError(f"Unsupported file format: {file_path}")


def detect_emission_column(gdf):
    """Auto-detect emission rate column from GeoDataFrame."""
    emission_keywords = ['emission', 'ch4', 'methane', 'ch4_kg', 'mean_ch4', 
                         'mean_emission', 'flux', 'rate', 'kgph', 'kg/h']
    
    for col in gdf.columns:
        col_lower = col.lower()
        if any(keyword in col_lower for keyword in emission_keywords):
            print(f"Detected emission column: {col}")
            return col
    
    raise ValueError(f"Could not auto-detect emission column. Available columns: {list(gdf.columns)}")


def ensure_lonlat(gdf):
    """Ensure GeoDataFrame has Latitude and Longitude columns."""
    if 'Latitude' not in gdf.columns or 'Longitude' not in gdf.columns:
        if gdf.geometry.name is not None:
            gdf['Longitude'] = gdf.geometry.x
            gdf['Latitude'] = gdf.geometry.y
    
    return gdf


def convert_to_feast(gdf, emission_col, output_dir='./'):
    """Convert to FEAST-ready CSV and pickle formats."""
    # Select relevant columns
    feast_cols = ['Longitude', 'Latitude', emission_col]
    feast_df = gdf[feast_cols].copy()
    feast_df.columns = ['Longitude', 'Latitude', 'emission_rate_kgph']
    
    # Save CSV
    csv_path = os.path.join(output_dir, 'feast_emissions.csv')
    feast_df.to_csv(csv_path, index=False)
    print(f"Saved CSV: {csv_path}")
    
    # Save pickle
    pkl_path = os.path.join(output_dir, 'feast_emissions.pkl')
    feast_df.to_pickle(pkl_path)
    print(f"Saved pickle: {pkl_path}")
    
    return feast_df


def main():
    """Main ingestion pipeline."""
    zip_path = 'PAMarginalWellEmissions2023.zip'
    
    if not os.path.exists(zip_path):
        print(f"Error: {zip_path} not found. Looking for .pkl directly...")
        if os.path.exists('PAMarginalWellEmissions2023.pkl'):
            zip_path = 'PAMarginalWellEmissions2023.pkl'
        else:
            print("Error: Could not find input file")
            return
    
    # Extract from ZIP if needed
    if zip_path.endswith('.zip'):
        extracted = extract_known_from_zip(zip_path)
        if 'gpkg' in extracted:
            data_path = extracted['gpkg']
        elif 'pkl' in extracted:
            data_path = extracted['pkl']
        else:
            print("Error: Could not find .gpkg or .pkl in ZIP")
            return
    else:
        data_path = zip_path
    
    # Load data
    print(f"Loading data from {data_path}...")
    gdf = load_data(data_path)
    print(f"Loaded {len(gdf)} records")
    
    # Detect emission column
    emission_col = detect_emission_column(gdf)
    
    # Ensure coordinates
    gdf = ensure_lonlat(gdf)
    
    # Convert to FEAST format
    feast_df = convert_to_feast(gdf, emission_col)
    
    # Summary statistics
    print("\n=== DATA SUMMARY ===")
    print(f"Total wells: {len(feast_df)}")
    print(f"Total emissions: {feast_df['emission_rate_kgph'].sum():.1f} kg/h")
    print(f"Mean emission: {feast_df['emission_rate_kgph'].mean():.3f} kg/h")
    print(f"Median emission: {feast_df['emission_rate_kgph'].median():.3f} kg/h")
    print(f"Max emission: {feast_df['emission_rate_kgph'].max():.3f} kg/h")
    print(f"Min emission: {feast_df['emission_rate_kgph'].min():.3f} kg/h")


if __name__ == '__main__':
    main()
