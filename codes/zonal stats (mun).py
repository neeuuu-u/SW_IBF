# -*- coding: utf-8 -*-
"""
Recursive NetCDF Plotter with Zonal Statistics
- Reads all .nc files recursively
- Masks rasters within shapefile polygons
- Calculates zonal statistics (mean per municipality)
- Saves both PNG plots and CSV zonal stats
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
root_dir = r"D:\neu\CIM oct28\tcrm\outputs (single_time)\opong (2025)"       # where .nc files are
output_root = r"C:\Users\jo_ht\OneDrive\Documents\neu\CIM oct28\outputs"   # output png folder
shapefile_path = r"C:\Users\jo_ht\OneDrive\Documents\neu\CIM oct28\shp\ph_admin_mun_boundaries\PHL_adm4_PSA_pn_2016Junprj_mun.shp"  # UPDATE THIS PATH

os.makedirs(output_root, exist_ok=True)

# Geographic extent
min_lon, max_lon = 116.56, 126.83
min_lat, max_lat = 4.58, 21.12

# Figure size
fig_width_mm, fig_height_mm = 158, 225
fig_width_in = fig_width_mm / 25.4
fig_height_in = fig_height_mm / 25.4

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
    """Beaufort scale color scheme"""
    thresholds = [
        0.72, 5.4, 11.8, 19.44, 28.44, 38.52, 49.68, 61.56, 
        74.52, 87.84, 102.24, 117.36, 132.84, 149.04, 165.96, 
        183.24, 201.6, 300
    ]
    colors = [
        "#ffffff", "#f2f2f2", "#dae9f8", "#c7eefd", "#61cbf3", "#119fd7",
        "#83e291", "#01ae50", "#fffc02", "#ffce01", "#ff7502", "#d01315",
        "#740f15", "#782172", "#175e85", "#0e2841", "#818180", "#000000"
    ]
    cmap = mcolors.ListedColormap(colors)
    return thresholds, cmap

# =========================
# Zonal statistics function
# =========================
def calculate_zonal_stats(data_array, gdf):
    """
    Calculate mean values within each polygon
    
    Parameters:
    - data_array: xarray DataArray with 'lat' and 'lon' coordinates
    - gdf: GeoDataFrame with polygons
    
    Returns:
    - DataFrame with zonal statistics
    """
    # Get coordinate info
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
            
            # Create mask for this polygon
            mask = geometry_mask(
                [mapping(geom)],
                transform=transform,
                invert=True,
                out_shape=(height, width)
            )
            
            # Apply mask to data
            masked_data = np.where(mask, data_array.values, np.nan)
            
            # Calculate statistics
            mean_val = np.nanmean(masked_data)
            max_val = np.nanmax(masked_data)
            min_val = np.nanmin(masked_data)
            std_val = np.nanstd(masked_data)
            count = np.sum(~np.isnan(masked_data))
            
            # Store results with all attributes from shapefile
            result = {col: row[col] for col in gdf.columns if col != 'geometry'}
            result.update({
                'mean_vmax': mean_val,
                'max_vmax': max_val,
                'min_vmax': min_val,
                'std_vmax': std_val,
                'pixel_count': count
            })
            
            results.append(result)
            
        except Exception as e:
            print(f"⚠️  Error processing polygon {idx}: {e}")
            continue
    
    return pd.DataFrame(results)

# =========================
# Plotting function
# =========================
def plot_nc_with_shapefile(nc_path, output_png, gdf):
    """Plot NetCDF with shapefile overlay"""
    print(f"📊 Plotting: {nc_path}")
    
    ds = xr.open_dataset(nc_path)
    
    if "vmax" not in ds:
        print("⚠️  Skipping (no 'vmax' variable):", nc_path)
        ds.close()
        return None
    
    gust = ds['vmax']
    
    # Get colorbar setup
    lev, cmap = wspd_lc()
    norm = mcolors.BoundaryNorm(lev, cmap.N)
    
    fig = plt.figure(figsize=(fig_width_in, fig_height_in))
    ax = fig.add_axes([0.05, 0.05, 0.9, 0.9], projection=ccrs.PlateCarree())
    
    # Plot raster
    pc = gust.plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap=cmap,
        norm=norm,
        add_colorbar=False
    )
    
    # Overlay shapefile boundaries
    gdf.boundary.plot(ax=ax, color='black', linewidth=0.5, alpha=0.7, transform=ccrs.PlateCarree())
    
    # Add colorbar
    cbax = fig.add_axes([0.1, 0.02, 0.8, 0.02])
    cb = plt.colorbar(
        pc,
        ticks=lev,
        orientation='horizontal',
        drawedges=True,
        extend='max',
        cax=cbax
    )
    cb.set_label(label='Wind Speed (kph)', size=14, labelpad=7)
    cb.ax.tick_params(length=0, direction='out', labelsize=13, pad=1)
    cb.outline.set_linewidth(1)
    
    # Map setup
    ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black', linewidth=0.8)
    ax.gridlines(draw_labels=False, linestyle='--', alpha=0.5)
    
    # Save PNG
    os.makedirs(os.path.dirname(output_png), exist_ok=True)
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.close()
    
    return gust

# =========================
# Main processing loop
# =========================
print("🚀 Starting recursive processing...")

for root, dirs, files in os.walk(root_dir):
    for f in files:
        if f.endswith(".nc"):
            nc_path = os.path.join(root, f)
            
            # Determine output paths
            rel_path = os.path.relpath(root, root_dir)
            output_folder = os.path.join(output_root, rel_path)
            os.makedirs(output_folder, exist_ok=True)
            
            base_name = os.path.splitext(f)[0]
            png_path = os.path.join(output_folder, f"{base_name}.png")
            csv_path = os.path.join(output_folder, f"{base_name}_zonal_stats.csv")
            
            try:
                # 1. Create plot
                gust_data = plot_nc_with_shapefile(nc_path, png_path, gdf)
                
                if gust_data is not None:
                    print(f"✅ Saved plot: {png_path}")
                    
                    # 2. Calculate zonal statistics
                    print(f"📈 Calculating zonal statistics...")
                    stats_df = calculate_zonal_stats(gust_data, gdf)
                    
                    # 3. Save CSV
                    stats_df.to_csv(csv_path, index=False, encoding='utf-8-sig')
                    print(f"✅ Saved zonal stats: {csv_path}")
                    print(f"   {len(stats_df)} municipalities processed\n")
                
                # Close dataset
                ds = xr.open_dataset(nc_path)
                ds.close()
                
            except Exception as e:
                print(f"❌ Error processing {nc_path}: {e}\n")
                continue

print("✅ All .nc files processed with plots and zonal statistics.")