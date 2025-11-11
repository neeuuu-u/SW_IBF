# -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 2025
Updated: Added dual colorbar logic for OBS/TCRM vs ANOM columns

@author: neu
"""

# ==========================================================
# 📦 IMPORTS
# ==========================================================
import os
from pathlib import Path
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from shapely.geometry import Point

# ==========================================================
# 🧭 USER CONFIG
# ==========================================================
TCRM_ROOT = Path(r"D:\neu\CIM oct28\hazard\png files\Opong [FST, with multipliers]")
OBS_ROOT  = Path(r"D:\neu\CIM oct28\observed [sycoder]")
OUTPUT_ROOT = Path(r"D:\neu\CIM oct28\zone3")
OUTPUT_SUFFIX = "_merged_obs"
KNOTS_TO_KPH = 1.852
OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

# ==========================================================
# 🎨 COLORBAR CONFIGURATIONS
# ==========================================================
def wspd_lc():
    thresholds = [
        1.852, 5.556, 11.112, 18.52, 29.632, 38.892,
        50.004, 61.116, 74.08, 87.044, 101.86, 116.676,
        131.492, 148.16, 164.828, 183.348, 200.016, 277.8
    ]
    colors = [
        "#ffffff", "#f2f2f2", "#dae9f8", "#c7eefd", "#61cbf3", "#119fd7",
        "#83e291", "#01ae50", "#fffc02", "#ffce01", "#ff7502", "#d01315",
        "#740f15", "#782172", "#175e85", "#0e2841", "#818180", "#000000"
    ]
    cmap = mcolors.ListedColormap(colors)
    return thresholds, cmap


def anomaly_lc():
    thresholds = [-50, -40, -30, -20, -10, -5, -1, 0, 1, 5, 10, 20, 30, 40, 50]
    colors = [
        "#440154", "#482878", "#3e4989", "#31688e", "#26828e", "#1f9e89",
        "#35b779", "#6ece58", "#b5de2b", "#fde725", "#ffd92f", "#ffa500",
        "#ff7500", "#ff4500", "#ff0000"
    ]
    cmap = mcolors.ListedColormap(colors)
    return thresholds, cmap

# ==========================================================
# 🗺️ SHAPEFILE LOADING (Municipalities)
# ==========================================================
MUNI_SHP = Path(r"D:\neu\CIM oct28\shp\ph_admin_mun_boundaries\PHL_adm4_PSA_pn_2016Junprj_mun.shp")
gdf_muni = gpd.read_file(MUNI_SHP)
gdf_muni["NAME_CLEAN"] = gdf_muni["Mun_Name"].str.strip().str.lower()

# ==========================================================
# ⚙️ HELPER FUNCTIONS
# ==========================================================
def clean_name(name):
    if pd.isna(name):
        return ""
    return str(name).strip().lower().replace("-", " ").replace(".", "").replace(",", "")

def match_municipality(df, col_name):
    df["NAME_CLEAN"] = df[col_name].apply(clean_name)
    merged = df.merge(gdf_muni, how="left", on="NAME_CLEAN")
    return merged

# ==========================================================
# 📊 LOAD AND MERGE TCRM & OBS DATA
# ==========================================================
def load_csvs():
    tcrm_files = list(TCRM_ROOT.glob("*.csv"))
    obs_files = list(OBS_ROOT.glob("*.csv"))
    tcrm_dfs, obs_dfs = {}, {}

    for f in tcrm_files:
        name = f.stem
        tcrm_dfs[name] = pd.read_csv(f)

    for f in obs_files:
        name = f.stem
        obs_dfs[name] = pd.read_csv(f)

    return tcrm_dfs, obs_dfs

def merge_data(tcrm_dfs, obs_dfs):
    merged_records = []
    for tname, tdf in tcrm_dfs.items():
        obs_key = next((o for o in obs_dfs.keys() if o in tname), None)
        if obs_key is None:
            continue

        odf = obs_dfs[obs_key]
        merged = pd.merge(tdf, odf, left_on="municipality", right_on="municipality", suffixes=("_TCRM", "_OBS"))
        merged_records.append(merged)

    if merged_records:
        merged_all = pd.concat(merged_records, ignore_index=True)
        merged_all = match_municipality(merged_all, "municipality")
        return merged_all
    else:
        return pd.DataFrame()

# ==========================================================
# 🗾 PLOTTING FUNCTION
# ==========================================================
def plot_maps(gdf, cols_to_plot):
    for col in cols_to_plot:
        if col not in gdf.columns:
            print(f"⚠️ Skipping missing column: {col}")
            continue

        # Select colorbar logic
        if col.startswith("anom_"):
            thresholds, cmap = anomaly_lc()
        else:
            thresholds, cmap = wspd_lc()

        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        gdf.plot(
            column=col,
            cmap=cmap,
            linewidth=0.5,
            ax=ax,
            edgecolor="gray",
            legend=True,
            legend_kwds={"label": col, "orientation": "vertical"},
        )

        ax.set_title(f"{col}", fontsize=13, fontweight="bold")
        ax.axis("off")

        output_file = OUTPUT_ROOT / f"{col}{OUTPUT_SUFFIX}.png"
        os.makedirs(OUTPUT_ROOT, exist_ok=True)
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"✅ Saved map: {output_file}")

# ==========================================================
# 🚀 MAIN EXECUTION
# ==========================================================
if __name__ == "__main__":
    print("🔹 Loading data...")
    tcrm_dfs, obs_dfs = load_csvs()

    print("🔹 Merging datasets...")
    merged = merge_data(tcrm_dfs, obs_dfs)

    if merged.empty:
        print("⚠️ No merged data generated.")
    else:
        print("🔹 Plotting maps...")

        COLS_TO_PLOT = [
            "TCRM_31", "TCRM_26", "OBS_9", "OBS_13", "anom_TCRM31_minus_OBS9",	"anom_TCRM26_minus_OBS13"
        ]

        plot_maps(merged, COLS_TO_PLOT)

    print("✅ Done.")
