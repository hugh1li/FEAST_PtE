#!/usr/bin/env python3
"""
Fetch L4 MethaneSat (area) from Google Earth Engine, subset to Marcellus bbox,
download the 4 km raster, and plot/save the result.

Usage:
 - Place a service account JSON at `ee-service-account.json` in the repo, OR
 - Run interactively to authenticate (will open browser). 

Outputs saved to `ExampleRunScriptResults/`.
"""

from pathlib import Path
import json
import zipfile
import io
import requests
import sys
import os

try:
    import ee
except Exception as e:
    print('Earth Engine python package not found. Install with: pip install earthengine-api')
    raise

import rasterio
import matplotlib.pyplot as plt
import numpy as np


def ee_init():
    # Priority order for credentials:
    # 1. Path set in env var `MSAT_SERVICE_ACCOUNT_JSON` or standard `GOOGLE_APPLICATION_CREDENTIALS`
    # 2. Local file `ee-service-account.json` in repo root
    # 3. Any service account JSON placed under `MyCodetoRun/*.json` (first match)
    # 4. Interactive auth
    env_paths = [
        ("MSAT_SERVICE_ACCOUNT_JSON", None),
        ("GOOGLE_APPLICATION_CREDENTIALS", None),
    ]

    # Check environment variables first
    for var, _ in env_paths:
        p = os.environ.get(var)
        if p:
            ppath = Path(p)
            if ppath.exists():
                info = json.loads(ppath.read_text())
                service_account = info.get('client_email')
                print(f'Initializing Earth Engine with service account from env {var}: {service_account}')
                credentials = ee.ServiceAccountCredentials(service_account, str(ppath))
                ee.Initialize(credentials)
                return

    # Next, check common repo locations
    candidates = [Path('ee-service-account.json')]
    # look for JSONs in MyCodetoRun directory
    mycodetorun = Path('MyCodetoRun')
    if mycodetorun.exists() and mycodetorun.is_dir():
        for p in mycodetorun.glob('*.json'):
            candidates.append(p)

    for key_path in candidates:
        if key_path.exists():
            try:
                info = json.loads(key_path.read_text())
                service_account = info.get('client_email')
                print(f'Initializing Earth Engine with service account: {service_account} (from {key_path})')
                credentials = ee.ServiceAccountCredentials(service_account, str(key_path))
                ee.Initialize(credentials)
                return
            except Exception as e:
                print(f'Found key at {key_path} but failed to initialize: {e}')

    # Fallback to interactive auth
    print('No service account JSON found. Running interactive authentication (opens browser).')
    ee.Authenticate()
    ee.Initialize()


def download_image(img, region, scale, out_dir):
    params = {
        'scale': int(scale),
        'crs': 'EPSG:4326',
        'region': region.getInfo()['coordinates']
    }
    print('Requesting download URL from Earth Engine...')
    url = img.getDownloadURL(params)
    print('Download URL obtained; downloading zip...')
    r = requests.get(url, stream=True)
    r.raise_for_status()

    z = zipfile.ZipFile(io.BytesIO(r.content))
    z.extractall(out_dir)
    print(f'Extracted files to {out_dir}')
    # find first .tif
    for p in Path(out_dir).glob('**/*.tif'):
        return str(p)
    raise FileNotFoundError('No TIFF found in downloaded archive')


def plot_tif(tif_path, out_png):
    with rasterio.open(tif_path) as src:
        arr = src.read(1)
        arr = arr.astype(float)
        arr[arr == src.nodata] = np.nan
        extent = (src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top)

    plt.figure(figsize=(8, 6))
    im = plt.imshow(arr, origin='upper', extent=extent, cmap='viridis')
    plt.colorbar(im, label='MethaneSat L4 area (unit-dependent)')
    plt.title('MethaneSat L4 area — Marcellus subset')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    print(f'Plot saved to {out_png}')


def main():
    ee_init()

    # Asset ID from provided link (public preview)
    asset_id = 'projects/edf-methanesat-ee_assets_public-preview/L4area_v2'
    print('Loading asset:', asset_id)
    img = ee.Image(asset_id)

    print('Bands:', img.bandNames().getInfo())

    # Marcellus bounding box (approximate / recommended)
    lon_min, lat_min = -80.6, 39.5
    lon_max, lat_max = -75.8, 42.5

    region = ee.Geometry.Rectangle([lon_min, lat_min, lon_max, lat_max])

    # Dataset native resolution ~4 km; request scale 4000 m
    scale = 4000

    out_dir = Path('ExampleRunScriptResults/methanesat_download')
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        tif = download_image(img, region, scale, out_dir)
    except Exception as e:
        print('Error downloading image from Earth Engine:', e)
        sys.exit(1)

    out_png = Path('ExampleRunScriptResults/methanesat_marcellus.png')
    plot_tif(tif, out_png)


if __name__ == '__main__':
    main()
