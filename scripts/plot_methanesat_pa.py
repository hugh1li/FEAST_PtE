#!/usr/bin/env python3
"""
Plot MethaneSat points and overlay the 4x4 km grid boxes for the PA Marcellus region.

Reads:
 - MyCodetoRun/feast_emissions.csv (satellite point emissions)
 - ExampleRunScriptResults/PA_Marginal_GridBoxes_methaneSat.csv (4x4 km aggregated boxes)

Output:
 - ExampleRunScriptResults/PA_MethaneSat_4x4km_overlay.png
"""

from pathlib import Path
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point, box
import numpy as np


def main():
    out_dir = Path('ExampleRunScriptResults')
    out_dir.mkdir(exist_ok=True)

    sat_csv = Path('MyCodetoRun/feast_emissions.csv')
    grid_csv = Path('ExampleRunScriptResults/PA_Marginal_GridBoxes_methaneSat.csv')

    if not sat_csv.exists():
        print('Satellite CSV not found at', sat_csv)
        return
    if not grid_csv.exists():
        print('Grid CSV not found at', grid_csv)
        return

    print('Loading satellite points...')
    sat_df = pd.read_csv(sat_csv)
    sat_gdf = gpd.GeoDataFrame(sat_df, geometry=gpd.points_from_xy(sat_df.Longitude, sat_df.Latitude), crs='EPSG:4326')

    print('Loading aggregated grid boxes...')
    grid_df = pd.read_csv(grid_csv)

    # compute constant lat_deg (approx) and per-row lon_deg using lat_center
    box_km = 4.0
    lat_deg = box_km / 111.0
    # Create polygons for each grid row
    polys = []
    for _, row in grid_df.iterrows():
        latc = float(row['lat_center'])
        lonc = float(row['lon_center'])
        lon_deg = box_km / (111.0 * np.cos(np.radians(latc)))
        minx = lonc - lon_deg/2.0
        maxx = lonc + lon_deg/2.0
        miny = latc - lat_deg/2.0
        maxy = latc + lat_deg/2.0
        polys.append(box(minx, miny, maxx, maxy))

    grid_gdf = gpd.GeoDataFrame(grid_df.copy(), geometry=polys, crs='EPSG:4326')

    # Prepare plotting
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))

    # Plot grid boxes colored by aggregated emissions (mean_ch4_kgh)
    grid_gdf.plot(column='mean_ch4_kgh', cmap='Reds', linewidth=0.2, edgecolor='gray', alpha=0.8, legend=True, ax=ax)

    # Overlay satellite points (small dots) colored by emission_rate_kgph
    sat_gdf.plot(ax=ax, column='emission_rate_kgph', cmap='viridis', markersize=6, alpha=0.6, legend=True)

    ax.set_title('MethaneSat points + 4x4 km aggregated boxes (PA Marcellus)')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')

    # Zoom to production bounds from grid
    minx, miny, maxx, maxy = grid_gdf.total_bounds
    padx = (maxx - minx) * 0.02
    pady = (maxy - miny) * 0.02
    ax.set_xlim(minx - padx, maxx + padx)
    ax.set_ylim(miny - pady, maxy + pady)

    out_png = out_dir / 'PA_MethaneSat_4x4km_overlay.png'
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    print('Saved overlay plot to', out_png)


if __name__ == '__main__':
    main()
