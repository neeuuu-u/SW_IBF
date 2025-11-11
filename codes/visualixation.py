# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 14:27:49 2025

@author: jo_ht
"""

# -*- coding: utf-8 -*-
"""
Recursive NetCDF Plotter with Custom Wind Speed Colorbar and Boundaries
- Windfield raster (vmax) plotting
- Province (thick black) and municipality (thin gray) overlay
- CRS unified
- Safe filenames and recursive output folder creation
"""

import os
import re
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
import cartopy.crs as ccrs

# =========================
# CONFIG
# =========================
root_dir = r"C:\Users\jo_ht\tcrm-3.1.16\output\archive\opong (2025)"
output_root = r"C:\Users\jo_ht\OneDrive\Documents\neu\CIM oct28\outputs"
os.makedirs(output_root, exist_ok=True)

mun_shp_path = r"C:\Users\jo_ht\OneDrive\Documents\neu\CIM oct28\shp\ph_admin_mun_boundaries\PHL_adm4_PSA_pn_2016Junprj_mun.shp"
prov_shp_path = r"C:\Users\jo_ht\OneDrive\Documents\neu\CIM oct28\shp\ph_admin_prov_boundaries\PHL_adm4_PSA_pn_2016Junprj_prov.shp"

# =========================
# Colorbar function
# =========================
def wspd_lc():
    thresholds = [
        1.852, 5.556, 11.112, 18.52, 29.632, 38.892, 50.004, 61.116,
        74.08, 87.044, 101.86, 116.676, 131.492, 148.16, 164.828,
        183.348, 200.016, 277.8
    ]
    colors = [
        "#ffffff", "#f2f2f2", "#dae9f8", "#c7eefd", "#61cbf3", "#119fd7",
        "#83e291", "#01ae50", "#fffc02", "#ffce01", "#ff7502", "#d01315",
        "#740f15", "#782172", "#175e85", "#0e2841", "#818180", "#000000"
    ]
    cmap = mcolors.ListedColormap(colors)
    return thresholds, cmap

# =========================
# Load Shapefiles (once)
# =========================
mun_gdf = gpd.read_file(mun_shp_path)
prov_gdf = gpd.read_file(prov_shp_path)

# =========================
# Plotting function
# =========================
def plot_nc(nc_path, output_png):
    print(f"📊 Plotting: {nc_path}")
    ds = xr.open_dataset(nc_path)
    if "vmax" not in ds:
        print(f"⚠️ Skipping (no 'vmax' variable): {nc_path}")
        ds.close()
        return

    gust = ds['vmax']

    # Use windfield CRS if defined, else default PlateCarree
    try:
        raster_crs = ccrs.epsg(gust.rio.crs.to_epsg())
    except:
        raster_crs = ccrs.PlateCarree()

    # Reproject shapefiles to raster CRS if needed
    try:
        mun_plot = mun_gdf.to_crs(gust.rio.crs)
        prov_plot = prov_gdf.to_crs(gust.rio.crs)
    except:
        mun_plot = mun_gdf
        prov_plot = prov_gdf

    # Colorbar setup
    lev, cmap = wspd_lc()
    norm = mcolors.BoundaryNorm(lev, cmap.N)

    # Figure
    fig = plt.figure(figsize=(6.22,8.86))  # original mm -> inch
    ax = fig.add_axes([0.05,0.05,0.9,0.9], projection=raster_crs)

    # Windfield raster
    pc = gust.plot(ax=ax, transform=ccrs.PlateCarree(), cmap=cmap, norm=norm, add_colorbar=False)

    # Add colorbar
    cbax = fig.add_axes([0.1,0.02,0.8,0.02])
    cb = plt.colorbar(pc, ticks=lev, orientation='horizontal', drawedges=True, extend='max', cax=cbax)
    cb.set_label(label='Wind Speed (kph)', size=14, labelpad=7)
    cb.ax.tick_params(length=0, direction='out', labelsize=13, pad=1)
    cb.outline.set_linewidth(1)

    # =========================
    # Plot boundaries
    # =========================
    prov_plot.boundary.plot(ax=ax, edgecolor='black', linewidth=1.5, zorder=2)
    mun_plot.boundary.plot(ax=ax, edgecolor='gray', linewidth=0.5, zorder=3)

    # Save PNG
    os.makedirs(os.path.dirname(output_png), exist_ok=True)
    safe_title = re.sub(r'[<>:"/\\|?*]', '_', os.path.splitext(os.path.basename(output_png))[0])
    plt.savefig(os.path.join(os.path.dirname(output_png), safe_title + ".png"), dpi=300, bbox_inches='tight')
    plt.close()
    ds.close()

# =========================
# Main loop: search .nc recursively
# =========================
for root, dirs, files in os.walk(root_dir):
    for f in files:
        if f.endswith(".nc"):
            nc_path = os.path.join(root, f)
            rel_path = os.path.relpath(root, root_dir)
            output_folder = os.path.join(output_root, rel_path)
            os.makedirs(output_folder, exist_ok=True)
            png_filename = os.path.splitext(f)[0] + ".png"
            png_path = os.path.join(output_folder, png_filename)
            plot_nc(nc_path, png_path)

print("✅ All .nc files converted to PNG with boundaries & custom colorbar.")
