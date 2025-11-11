# -*- coding: utf-8 -*-
"""
Recursive NetCDF Plotter with Custom Wind Speed Colorbar (Vertical)
- Reads all .nc files recursively
- Uses vmax variable and your wspd_lc() color scheme
- Saves PNGs with same folder structure
"""

import os
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
from matplotlib.ticker import FuncFormatter

# =========================
# CONFIG
# =========================
root_dir = r"C:\Users\jo_ht\tcrm-3.1.16\output\UWAN_110825_0000\windfield\gust\instantaneous\evolution.001-00001"    # where .nc files are
output_root = r"Z:\Severe_Wind\Data\2025\UWAN_2025\PAGASA FCST track\instantaneous"   # output png folder
os.makedirs(output_root, exist_ok=True)

# Geographic extent
min_lon, max_lon = 116.56, 126.83
min_lat, max_lat = 4.58, 21.12

# Figure size (same as your script)
fig_width_mm, fig_height_mm = 158, 225
fig_width_in = fig_width_mm / 25.4
fig_height_in = fig_height_mm / 25.4

# =========================
# Colorbar function
# =========================
def wspd_lc():
    # Beaufort scale thresholds (in kph)
    thresholds = [
        0.72, 5.4, 11.8, 19.44, 28.44, 38.52, 49.68, 61.56,
        74.52, 87.84, 102.24, 117.36, 132.84, 149.04, 165.96,
        183.24, 201.6, 300
    ]

    # Hex colors per Beaufort level
    colors = [
        "#ffffff", "#f2f2f2", "#dae9f8", "#c7eefd", "#61cbf3", "#119fd7",
        "#83e291", "#01ae50", "#fffc02", "#ffce01", "#ff7502", "#d01315",
        "#740f15", "#782172", "#175e85", "#0e2841", "#818180", "#000000"
    ]

    cmap = mcolors.ListedColormap(colors)
    return thresholds, cmap


# =========================
# Plotting function
# =========================
def plot_nc(nc_path, output_png):
    print(f"📊 Plotting: {nc_path}")

    ds = xr.open_dataset(nc_path)
    if "vmax" not in ds:
        print("⚠️  Skipping (no 'vmax' variable):", nc_path)
        ds.close()
        return

    gust = ds['vmax']

    # Get colorbar setup
    lev, cmap = wspd_lc()
    norm = mcolors.BoundaryNorm(lev, cmap.N)

    # Create figure and axis
    fig = plt.figure(figsize=(fig_width_in, fig_height_in))
    ax = fig.add_axes([0.05, 0.05, 0.82, 0.9], projection=ccrs.PlateCarree())

    # Plot data
    pc = gust.plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap=cmap,
        norm=norm,
        add_colorbar=False
    )

    # Add vertical colorbar (right side)
    cbax = fig.add_axes([0.88, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
    cb = plt.colorbar(
        pc,
        ticks=lev,
        orientation='vertical',
        drawedges=True,
        extend='max',
        cax=cbax
    )
    cb.set_label(label='Wind Speed (kph)', size=14, labelpad=10)
    cb.ax.tick_params(length=0, direction='out', labelsize=12, pad=2)
    cb.outline.set_linewidth(1)

    # Map setup
    ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m', color='black', linewidth=0.8)
    ax.gridlines(draw_labels=False, linestyle='--', alpha=0.5)

    # Save PNG
    os.makedirs(os.path.dirname(output_png), exist_ok=True)
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
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

print("✅ All .nc files converted to PNG with custom vertical colorbar & ECMWF style.")
