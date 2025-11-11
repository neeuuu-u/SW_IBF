# -*- coding: utf-8 -*-
"""
TCRM + OBS Merger (Final Version)
---------------------------------
- Recursively reads TCRM gust zonal stats CSVs and OBS CSVs
- Matches files by datetime in filenames
- Keeps only gust-related columns from TCRM
- Merges with OBS data by Mun_Name + datetime
- Saves merged CSVs into output folder
- Skips files that don't match or lack required columns

@version: 2025-10-30
@author: neu
"""

import os
import re
import pandas as pd

# ==============================
# USER SETTINGS
# ==============================
TCRM_ROOT = r"D:\neu\CIM oct28\BST\nc\opong_bst (swath)\0923"
OBS_ROOT  = r"D:\neu\CIM oct28\BST\obs_csvs"
OUTPUT_DIR = r"D:\neu\CIM oct28\BST\merged_output"

# Columns to keep from TCRM
TCRM_USE_COLS = [
    "Mun_Name",
    "Pro_Name",
    "gust_speed_mean",
    "gust_speed_max",
    "gust_speed_min",
    "gust_speed_std",
    "gust_speed_pixel_count",
]

# Ensure output folder exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ==============================
# Helper Functions
# ==============================

def find_csv_files(root_dir):
    """Recursively find all CSV files under root_dir."""
    csvs = []
    for root, _, files in os.walk(root_dir):
        for f in files:
            if f.endswith(".csv"):
                csvs.append(os.path.join(root, f))
    return csvs


def extract_datetime_from_filename(filename):
    """
    Extract datetime from filename (e.g., '_2025-09-25_0600_')
    Returns string like '2025-09-25_0600' or None if not found.
    """
    match = re.search(r"\d{4}-\d{2}-\d{2}_\d{4}", filename)
    return match.group(0) if match else None


# ==============================
# Step 1: Load OBS CSVs
# ==============================

print("📂 Scanning OBS CSVs...")
obs_files = find_csv_files(OBS_ROOT)
obs_groups = {}

for obs_path in obs_files:
    dt_key = extract_datetime_from_filename(obs_path)
    if dt_key:
        obs_groups.setdefault(dt_key, []).append(obs_path)

print(f"✅ Found {len(obs_groups)} unique OBS datetime groups")

# ==============================
# Step 2: Load and process TCRM CSVs
# ==============================

tcrm_files = find_csv_files(TCRM_ROOT)
print(f"Found {len(tcrm_files)} TCRM CSVs to process\n")

merged_count = 0
skipped_no_obs = 0
skipped_bad_cols = 0

for tcrm_path in tcrm_files:
    tcrm_name = os.path.basename(tcrm_path)
    dt_key = extract_datetime_from_filename(tcrm_name)

    if not dt_key:
        print(f"⚠️  Skipping {tcrm_name}: no datetime found in filename")
        continue

    try:
        # Read with tab delimiter
        tcrm = pd.read_csv(tcrm_path, sep="\t", engine="python", encoding="utf-8", comment="#")

        # Debug: print columns if missing
        missing_cols = [c for c in TCRM_USE_COLS if c not in tcrm.columns]
        if missing_cols:
            print(f"❌ Error in {tcrm_name}: Missing columns {missing_cols}")
            print(f"   Detected columns: {tcrm.columns.tolist()}")
            skipped_bad_cols += 1
            continue

        tcrm_selected = tcrm[TCRM_USE_COLS].copy()
        tcrm_selected["datetime"] = dt_key

        # Find matching OBS file
        obs_match = obs_groups.get(dt_key)
        if not obs_match:
            print(f"⚠️  No OBS match for {tcrm_name}")
            skipped_no_obs += 1
            continue

        # Merge all OBS files for this datetime
        obs_dfs = []
        for obs_path in obs_match:
            try:
                obs_df = pd.read_csv(obs_path, encoding="utf-8")
                obs_dfs.append(obs_df)
            except Exception as e:
                print(f"   ⚠️  Skipping bad OBS file: {obs_path} ({e})")

        if not obs_dfs:
            print(f"⚠️  No valid OBS data for {tcrm_name}")
            skipped_no_obs += 1
            continue

        obs_all = pd.concat(obs_dfs, ignore_index=True)
        obs_all["datetime"] = dt_key

        # Merge on Mun_Name
        merged = pd.merge(obs_all, tcrm_selected, on="Mun_Name", how="inner", suffixes=("_obs", "_tcrm"))

        if merged.empty:
            print(f"⚠️  No matching municipalities for {tcrm_name}")
            continue

        # Save merged CSV
        out_name = f"merged_{dt_key}.csv"
        out_path = os.path.join(OUTPUT_DIR, out_name)
        merged.to_csv(out_path, index=False, encoding="utf-8-sig")
        print(f"✅ Created: {out_name}")
        merged_count += 1

    except Exception as e:
        print(f"❌ Error processing {tcrm_name}: {e}")

# ==============================
# Step 3: Summary
# ==============================

print("\n----- SUMMARY -----")
print(f"✅ Created merged files: {merged_count}")
print(f"⚠️ Skipped (no OBS match): {skipped_no_obs}")
print(f"⚠️ Skipped (bad columns): {skipped_bad_cols}")
print("Done.")
