# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 09:33:39 2025

@author: jo_ht
"""

# -*- coding: utf-8 -*-
"""
Recursive NetCDF Plotter with Multi-Variable Zonal Statistics
- Reads all .nc files recursively
- Masks rasters within shapefile polygons
- Calculates zonal statistics for: velocity_east, velocity_north, pressure, gust_speed
- Saves PNG plots ONLY for gust_speed, but comprehensive CSV with all variables
"""
import os
import xarray as xr
import numpy as np
import pandas as pd
import geopandas as gpd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from shapely.geometry import mapping
from rasterio.features import geometry_mask
from rasterio.transform import from_bounds

# =========================
# CONFIG
# =========================
root_dir = r"C:\Users\jo_ht\tcrm-3.1.16\output\UWAN_110825_0000\windfield\gust\instantaneous\evolution.001-00001"
output_root = r"D:\neu\CIM nov26\vizualized (single-time)\UWAN_110725_0000"   
shapefile_path = r"D:\neu\CIM oct28\shp\ph_admin_mun_boundaries\PHL_adm4_PSA_pn_2016Junprj_mun.shp"

os.makedirs(output_root, exist_ok=True)

# Geographic extent
min_lon, max_lon = 116.56, 126.83
min_lat, max_lat = 4.58, 21.12

# Figure size
fig_width_mm, fig_height_mm = 158, 225
fig_width_in = fig_width_mm / 25.4
fig_height_in = fig_height_mm / 25.4

# Variables to process for STATISTICS (all variables)
VARIABLES = ['velocity_east', 'velocity_north', 'pressure', 'gust_speed']

# Variable to PLOT (only gust_speed)
PLOT_VARIABLE = 'gust_speed'

# Variable-specific plotting configuration for gust_speed
VAR_CONFIG = {
    'gust_speed': {
        'label': 'Gust Speed (m/s)',
        'cmap': 'custom',
        'use_custom_cmap': True
    }
}

# =========================
# Load shapefile once
# =========================
print(f"📂 Loading shapefile: {shapefile_path}")
gdf = gpd.read_file(shapefile_path)
print(f"✅ Loaded {len(gdf)} polygons")
print(f"   Columns: {list(gdf.columns)}")

# =========================
# Colorbar function
# =========================
def wspd_lc():
    """Beaufort scale color scheme for wind speed (converted to m/s)"""
    # Original thresholds in kph, converted to m/s (divide by 3.6)
    thresholds_kph = [
        0.72, 5.4, 11.8, 19.44, 28.44, 38.52, 49.68, 61.56, 
        74.52, 87.84, 102.24, 117.36, 132.84, 149.04, 165.96, 
        183.24, 201.6, 300
    ]
    thresholds_ms = [t / 3.6 for t in thresholds_kph]
    
    colors = [
        "#ffffff", "#f2f2f2", "#dae9f8", "#c7eefd", "#61cbf3", "#119fd7",
        "#83e291", "#01ae50", "#fffc02", "#ffce01", "#ff7502", "#d01315",
        "#740f15", "#782172", "#175e85", "#0e2841", "#818180", "#000000"
    ]
    cmap = mcolors.ListedColormap(colors)
    return thresholds_ms, cmap

# =========================
# Multi-variable zonal statistics function
# =========================
def calculate_multi_zonal_stats(dataset, gdf, variables):
    """
    Calculate mean values within each polygon for multiple variables
    
    Parameters:
    - dataset: xarray Dataset with multiple variables
    - gdf: GeoDataFrame with polygons
    - variables: list of variable names to process
    
    Returns:
    - DataFrame with zonal statistics for all variables
    """
    # Get first available variable to extract coordinates
    first_var = None
    for var in variables:
        if var in dataset:
            first_var = var
            break
    
    if first_var is None:
        print("⚠️  No valid variables found in dataset")
        return None
    
    # Get coordinate info from first variable
    data_array = dataset[first_var]
    lats = data_array.lat.values
    lons = data_array.lon.values
    
    # Create affine transform for the raster
    height, width = len(lats), len(lons)
    transform = from_bounds(
        lons.min(), lats.min(), lons.max(), lats.max(),
        width, height
    )
    
    # Ensure CRS compatibility
    if gdf.crs is None:
        gdf = gdf.set_crs('EPSG:4326')
    elif gdf.crs != 'EPSG:4326':
        gdf = gdf.to_crs('EPSG:4326')
    
    results = []
    
    for idx, row in gdf.iterrows():
        try:
            geom = row.geometry
            
            # Create mask for this polygon (only once per polygon)
            mask = geometry_mask(
                [mapping(geom)],
                transform=transform,
                invert=True,
                out_shape=(height, width)
            )
            
            # Initialize result dict with shapefile attributes
            result = {col: row[col] for col in gdf.columns if col != 'geometry'}
            
            # Calculate statistics for each variable
            for var_name in variables:
                if var_name not in dataset:
                    continue
                
                var_data = dataset[var_name]
                
                # Apply mask to data
                masked_data = np.where(mask, var_data.values, np.nan)
                
                # Calculate statistics
                mean_val = np.nanmean(masked_data)
                max_val = np.nanmax(masked_data)
                min_val = np.nanmin(masked_data)
                std_val = np.nanstd(masked_data)
                count = np.sum(~np.isnan(masked_data))
                
                # Store results with variable prefix
                result.update({
                    f'{var_name}_mean': mean_val,
                    f'{var_name}_max': max_val,
                    f'{var_name}_min': min_val,
                    f'{var_name}_std': std_val,
                    f'{var_name}_pixel_count': count
                })
            
            results.append(result)
            
        except Exception as e:
            print(f"⚠️  Error processing polygon {idx}: {e}")
            continue
    
    return pd.DataFrame(results)

# =========================
# Plotting function for gust_speed
# =========================
def plot_gust_speed_with_shapefile(dataset, output_png, gdf):
    """Plot gust_speed variable with shapefile overlay"""
    
    var_name = PLOT_VARIABLE
    
    if var_name not in dataset:
        print(f"⚠️  Variable '{var_name}' not found, skipping plot")
        return
    
    print(f"📊 Plotting {var_name}: {output_png}")
    
    data_array = dataset[var_name]
    config = VAR_CONFIG[var_name]
    
    # Setup colormap
    lev, cmap = wspd_lc()
    norm = mcolors.BoundaryNorm(lev, cmap.N)
    
    fig = plt.figure(figsize=(fig_width_in, fig_height_in))
    ax = fig.add_axes([0.05, 0.05, 0.85, 0.9], projection=ccrs.PlateCarree())
    
    # Plot raster
    pc = data_array.plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap=cmap,
        norm=norm,
        add_colorbar=False
    )
    
    # Overlay shapefile boundaries
    gdf.boundary.plot(ax=ax, color='black', linewidth=0.5, alpha=0.7, transform=ccrs.PlateCarree())
    
    # Add VERTICAL colorbar with smaller font
    cbax = fig.add_axes([0.92, 0.1, 0.02, 0.8])  # [left, bottom, width, height]
    cb = plt.colorbar(
        pc,
        ticks=lev,
        orientation='vertical',  # Changed to vertical
        drawedges=True,
        extend='max',
        cax=cbax
    )
    
    cb.set_label(label=config['label'], size=10, labelpad=10)  # Reduced font size
    cb.ax.tick_params(length=0, direction='out', labelsize=8, pad=1)  # Reduced tick label size
    cb.outline.set_linewidth(1)
    
    # Map setup
    ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black', linewidth=0.8)
    ax.gridlines(draw_labels=False, linestyle='--', alpha=0.5)
    
    # Save PNG
    os.makedirs(os.path.dirname(output_png), exist_ok=True)
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Saved plot: {output_png}")

# =========================
# Main processing loop
# =========================
print("🚀 Starting recursive processing...")
print(f"📊 Plotting ONLY: {PLOT_VARIABLE}")
print(f"📋 CSV statistics for: {', '.join(VARIABLES)}\n")

for root, dirs, files in os.walk(root_dir):
    for f in files:
        if f.endswith(".nc"):
            nc_path = os.path.join(root, f)
            
            # Determine output paths
            rel_path = os.path.relpath(root, root_dir)
            output_folder = os.path.join(output_root, rel_path)
            os.makedirs(output_folder, exist_ok=True)
            
            base_name = os.path.splitext(f)[0]
            csv_path = os.path.join(output_folder, f"{base_name}_multi_zonal_stats.csv")
            png_path = os.path.join(output_folder, f"{base_name}_{PLOT_VARIABLE}.png")
            
            try:
                print(f"\n📂 Processing: {nc_path}")
                
                # Open dataset once for all operations
                ds = xr.open_dataset(nc_path)
                
                # 1. Create plot ONLY for gust_speed
                print(f"📊 Creating plot for {PLOT_VARIABLE}...")
                plot_gust_speed_with_shapefile(ds, png_path, gdf)
                
                # 2. Calculate zonal statistics for ALL variables
                print(f"📈 Calculating zonal statistics for {len(VARIABLES)} variables...")
                stats_df = calculate_multi_zonal_stats(ds, gdf, VARIABLES)
                
                if stats_df is not None and not stats_df.empty:
                    # 3. Save CSV
                    stats_df.to_csv(csv_path, index=False, encoding='utf-8-sig')
                    print(f"✅ Saved multi-variable zonal stats: {csv_path}")
                    print(f"   {len(stats_df)} municipalities processed")
                    print(f"   Columns: {len(stats_df.columns)} total")
                else:
                    print(f"⚠️  No statistics generated for {nc_path}")
                
                # Close dataset
                ds.close()
                
            except Exception as e:
                print(f"❌ Error processing {nc_path}: {e}")
                continue

print("\n✅ All .nc files processed!")
print(f"   📊 Plots created: gust_speed only")
print(f"   📋 CSV files contain: all {len(VARIABLES)} variables")