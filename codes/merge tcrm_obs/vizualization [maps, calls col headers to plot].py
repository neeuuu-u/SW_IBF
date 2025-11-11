# -*- coding: utf-8 -*-
"""
Municipality-Based Visualization for TCRM + OBS CSVs
----------------------------------------------------
- Reads CSVs recursively (headerless or with header)
- Merges with municipality-level shapefile
- Plots selected columns by name
- Saves one PNG per variable in organized folders

@version: 2025-10-29
@author: jo_ht
"""

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
from pathlib import Path

# ==========================================================
# 🧭 USER CONFIG
# ==========================================================
SHP_PATH = Path(
    r"D:\neu\CIM oct28\shp\ph_admin_mun_boundaries\PHL_adm4_PSA_pn_2016Junprj_mun.shp"
)
CSV_ROOT = Path(r"D:\neu\CIM oct28\zones [final]")
OUTPUT_ROOT = Path(r"D:\neu\CIM oct28\zones_maps_muni")
OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

OUTPUT_SUFFIX = "_maps"

# ---------------------------
# Columns by name
# ---------------------------
TCRM_PROVINCE_COLNAME = "TCRM_2"
TCRM_MUNICIPALITY_COLNAME = "TCRM_7"

# Columns to plot (user-defined)
COLS_TO_PLOT = [
    "TCRM_31",     # gust anomaly
    "TCRM_26",    # pressure anomaly
    "OBS_9",                # obs u-component
    "OBS_13",                # obs v-component
]


COLS_TO_PLOT_anom = [
    "anom_TCRM31_minus_OBS9",
    "anom_TCRM26_minus_OBS13"
]


# Optional readable folder names
COL_FOLDERS = {
    "anom_TCRM31_minus_OBS9": "Anom_Gust",
    "anom_TCRM26_minus_OBS13": "Anom_Pressure",
    "OBS_7_Ucomp": "Obs_Ucomp",
    "OBS_7_Vcomp": "Obs_Vcomp",
}


# ==========================================================
# 🎨 COLOR MAPS
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


def generic_anom_lc():
    levels = list(range(-80, 85, 5))
    base_cmap = mpl.colormaps["RdBu_r"]
    colors = [base_cmap(i / (len(levels) - 1)) for i in range(len(levels) - 1)]
    cmap = mpl.colors.ListedColormap(colors)
    return levels, cmap


# ==========================================================
# 🧰 HELPERS
# ==========================================================
def normalize_text(s):
    """Normalize text for name matching."""
    if pd.isna(s):
        return ""
    return (
        str(s)
        .strip()
        .upper()
        .replace("-", " ")
        .replace("Ñ", "N")
        .replace("MUNICIPALITY OF ", "")
        .replace("CITY OF ", "")
        .replace("CITY", "")
        .strip()
    )


def plot_and_save_map(gdf, column, cmap, norm, folder_name, output_file):
    """Plot GeoDataFrame column and save PNG."""
    fig, ax = plt.subplots(figsize=(10, 10))
    gdf.plot(
        column=column,
        cmap=cmap,
        norm=norm,
        linewidth=0.2,
        edgecolor="black",
        legend=True,
        ax=ax,
    )
    ax.set_title(folder_name.replace("_", " "), fontsize=14)
    ax.axis("off")
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ==========================================================
# 🚀 MAIN SCRIPT
# ==========================================================
def main():
    print("📍 Loading municipality shapefile...")
    gdf = gpd.read_file(SHP_PATH)
    print(f"✅ Shapefile loaded: {gdf.shape}")

    # Normalize shapefile columns
    gdf["prov_norm"] = gdf["Pro_Name"].apply(normalize_text)
    gdf["muni_norm"] = gdf["Mun_Name"].apply(normalize_text)

    wspd_levels, wspd_cmap = wspd_lc()
    anom_levels, anom_cmap = generic_anom_lc()

    # ---- Find CSV files ----
    print(f"🔍 Searching CSV files in {CSV_ROOT} ...")
    csv_files = list(CSV_ROOT.rglob("*.csv"))
    if not csv_files:
        raise FileNotFoundError(f"No CSV files found under {CSV_ROOT}")
    print(f"📂 Found {len(csv_files)} CSV files")

    for csv_file in csv_files:
        try:
            df = pd.read_csv(csv_file)
            if df.empty:
                print(f"⚠️ Empty CSV skipped: {csv_file.name}")
                continue

            print(f"\n🗺 Processing: {csv_file.name}")

            # Normalize CSV location columns
            df["prov_norm"] = df[TCRM_PROVINCE_COLNAME].apply(normalize_text)
            df["muni_norm"] = df[TCRM_MUNICIPALITY_COLNAME].apply(normalize_text)

            # Merge shapefile + CSV data
            merged = gdf.merge(df, on=["prov_norm", "muni_norm"], how="left")

            # Track unmatched municipalities
            unmatched = merged[merged.isna().any(axis=1)]
            if len(unmatched) > 0:
                log_path = OUTPUT_ROOT / f"{csv_file.stem}_unmatched.txt"
                unmatched_names = unmatched[["Pro_Name", "Mun_Name"]]
                unmatched_names.to_csv(log_path, index=False)
                print(f"⚠️ Logged {len(unmatched_names)} unmatched municipalities → {log_path.name}")

            for col in COLS_TO_PLOT:
                if col not in df.columns:
                    print(f"⚠️ Column {col} not found → skip")
                    continue

                col_series = pd.to_numeric(merged[col], errors="coerce")
                if col_series.isna().all():
                    print(f"⚠️ Empty data in {col}, skipped.")
                    continue

                # Choose color map
                if "gust" in col.lower():
                    thresholds, cmap = wspd_levels, wspd_cmap
                    n_colors = len(cmap.colors)
                    t_min, t_max = thresholds[0], thresholds[-1]
                    new_thresholds = np.linspace(t_min, t_max, n_colors + 1)
                    norm = mcolors.BoundaryNorm(new_thresholds, ncolors=n_colors, clip=False)
                else:
                    thresholds, cmap = anom_levels, anom_cmap
                    norm = mcolors.BoundaryNorm(thresholds, ncolors=len(cmap.colors), clip=False)

                folder_name = COL_FOLDERS.get(col, col)
                out_dir = OUTPUT_ROOT / folder_name
                out_dir.mkdir(parents=True, exist_ok=True)
                out_file = out_dir / f"{csv_file.stem}_{folder_name}_map.png"

                merged[f"col_{col}"] = col_series
                plot_and_save_map(merged, f"col_{col}", cmap, norm, folder_name, out_file)
                print(f"✅ Saved: {out_file}")

        except Exception as e:
            print(f"❌ Error in {csv_file.name}: {e}")
            import traceback
            traceback.print_exc()

    print("\n🎉 Done generating all maps!")


# ==========================================================
# 🏁 ENTRY POINT
# ==========================================================
if __name__ == "__main__":
    main()
