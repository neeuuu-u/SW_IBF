# -*- coding: utf-8 -*-
"""
TCRM Zonal Stats + Observed Synop Merger (Fully index-based, configurable columns)
--------------------------------------------------------------------------------
- All columns referenced by index only (0-based)
- Specify columns to keep in output and columns for anomalies
- No reliance on CSV headers
"""

import os
import re
import pandas as pd
from pathlib import Path
import unicodedata

# ==========================================================
# 🧭 USER CONFIG
# ==========================================================
TCRM_ROOT = Path(r"D:\neu\CIM oct28\hazard\png files\Opong [FST, with multipliers]")
OBS_ROOT  = Path(r"D:\neu\CIM oct28\observed [sycoder]")
OUTPUT_ROOT = Path(r"D:\neu\CIM oct28\zonal stats")
OUTPUT_SUFFIX = "_merged_obs"
KNOTS_TO_KPH = 1.852
OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

# ---------------------------
# Column indices (0-based)
# ---------------------------
# TCRM columns
TCRM_PROVINCE_COL = 2
TCRM_MUNICIPALITY_COL = 7
TCRM_USE_COLS = [2, 7, 31, 16, 21, 26]  # province, municipality, gust, u, v, pressure

# OBS columns
OBS_SHORTNAME_COL = 0  # contains "Municipality, Province"
OBS_USE_COLS = [9, 7, 13]  # gust, direction, mslp

# ---------------------------
# Anomaly calculation pairs: (TCRM_col_index, OBS_col_index)
# ---------------------------
ANOMALY_PAIRS = [
    (31, 9),   # model gust (TCRM col 31) - obs gust (OBS col 9)
    (26, 13),  # model pressure (TCRM col 26) - obs mslp (OBS col 13)
]

# ==========================================================
# Helpers
# ==========================================================
def normalize_text(s):
    if pd.isna(s):
        return ''
    s = str(s).strip().lower()
    s = ''.join(c for c in unicodedata.normalize('NFD', s) if unicodedata.category(c) != 'Mn')
    s = re.sub(r'[^a-z0-9\s]', ' ', s)
    s = re.sub(r'\s+', ' ', s)
    return s.strip('_')

def parse_tcrm_datetime(fname):
    m = re.search(r'_(\d{4}-\d{2}-\d{2})_(\d{4})', fname)
    if not m:
        return None
    date_str = m.group(1).replace('-', '')
    time_str = m.group(2)
    return f"{date_str}_{time_str}"

def parse_obs_datetime(fname):
    m = re.search(r'obs_synop_(\d{8})_(\d{4})', fname)
    return f"{m.group(1)}_{m.group(2)}" if m else None

# ==========================================================
# Scan OBS CSVs
# ==========================================================
print("📂 Scanning OBS CSVs...")
obs_by_dt = {}
for obs_file in OBS_ROOT.rglob("*.csv"):
    dt_key = parse_obs_datetime(obs_file.name)
    if dt_key:
        obs_by_dt.setdefault(dt_key, []).append(obs_file)
print(f"✅ Found {len(obs_by_dt)} unique OBS datetime groups")

# ==========================================================
# Process TCRM CSVs
# ==========================================================
tcrm_files = list(TCRM_ROOT.rglob("*.csv"))
print(f"Found {len(tcrm_files)} TCRM CSVs to process")

created = skipped = 0

for tcrm_file in tcrm_files:
    try:
        dt_key = parse_tcrm_datetime(tcrm_file.name)
        if not dt_key:
            print(f"⚠️ Skipping (no date/time pattern): {tcrm_file.name}")
            continue

        obs_files = obs_by_dt.get(dt_key, [])
        if not obs_files:
            print(f"⚠️ No OBS match for {dt_key} → skip {tcrm_file.name}")
            skipped += 1
            continue

        tcrm = pd.read_csv(tcrm_file, header=None)
        if tcrm.empty:
            continue

        # Safety check: ensure all requested TCRM columns exist
        max_tcrm_col = max(TCRM_USE_COLS + [pair[0] for pair in ANOMALY_PAIRS])
        if max_tcrm_col >= tcrm.shape[1]:
            print(f"⚠️ Not enough columns in {tcrm_file.name} (have {tcrm.shape[1]}, need {max_tcrm_col + 1})")
            continue

        # Extract TCRM columns
        tcrm_selected = tcrm.iloc[:, TCRM_USE_COLS].copy()
        tcrm_selected.columns = [f"TCRM_{i}" for i in TCRM_USE_COLS]
        
        # Add normalized names for merge
        tcrm_selected["Province_norm"] = tcrm.iloc[:, TCRM_PROVINCE_COL].astype(str).apply(normalize_text)
        tcrm_selected["Municipality_norm"] = tcrm.iloc[:, TCRM_MUNICIPALITY_COL].astype(str).apply(normalize_text)

        # Add model columns for anomalies (keep original column index in name)
        for model_idx, _ in ANOMALY_PAIRS:
            tcrm_selected[f"model_col{model_idx}"] = pd.to_numeric(tcrm.iloc[:, model_idx], errors="coerce")

        # Combine OBS frames
        obs_frames = []
        for of in obs_files:
            obs = pd.read_csv(of, header=None)
            max_obs_col = max(OBS_USE_COLS + [OBS_SHORTNAME_COL] + [pair[1] for pair in ANOMALY_PAIRS])
            if obs.empty or max_obs_col >= obs.shape[1]:
                print(f"⚠️ Not enough columns in {of.name} (have {obs.shape[1]}, need {max_obs_col + 1})")
                continue

            obs_selected = obs.iloc[:, OBS_USE_COLS].copy()
            obs_selected.columns = [f"OBS_{i}" for i in OBS_USE_COLS]
            
            # Split short_name (column OBS_SHORTNAME_COL) into municipality + province
            short_series = obs.iloc[:, OBS_SHORTNAME_COL].astype(str)
            muni_list, prov_list = [], []
            for s in short_series:
                parts = [p.strip() for p in s.split(",")]
                if len(parts) == 1:
                    muni_list.append(parts[0])
                    prov_list.append('')
                else:
                    muni_list.append(parts[0])
                    prov_list.append(parts[-1])
            
            obs_selected["Municipality_norm"] = [normalize_text(m) for m in muni_list]
            obs_selected["Province_norm"] = [normalize_text(p) for p in prov_list]

            # Add observed columns for anomalies (keep original column index in name)
            for _, obs_idx in ANOMALY_PAIRS:
                col_name = f"obs_col{obs_idx}"
                val = pd.to_numeric(obs.iloc[:, obs_idx], errors="coerce")
                # Apply conversion if it's wind speed (col 9)
                if obs_idx == 9:
                    val = val * KNOTS_TO_KPH
                obs_selected[col_name] = val

            obs_frames.append(obs_selected)

        if not obs_frames:
            print(f"⚠️ No valid OBS data for {tcrm_file.name}")
            continue

        obs_all = pd.concat(obs_frames, ignore_index=True)

        # Merge strategy: Try exact match first, then province-only
        # 1. Exact municipality + province match
        merged = pd.merge(
            tcrm_selected, obs_all,
            on=["Municipality_norm", "Province_norm"], 
            how="left", 
            suffixes=('', '_obs')
        )

        # 2. For rows without match, try province-only
        no_match_mask = merged[[c for c in merged.columns if c.startswith('OBS_')]].isna().all(axis=1)
        if no_match_mask.any():
            prov_only = pd.merge(
                tcrm_selected[no_match_mask].drop(columns=['Municipality_norm']),
                obs_all.drop(columns=['Municipality_norm']),
                on="Province_norm",
                how="left",
                suffixes=('', '_prov')
            )
            # Update merged with province matches
            for col in [c for c in prov_only.columns if c.startswith('OBS_') or c.startswith('obs_col')]:
                if col in merged.columns:
                    merged.loc[no_match_mask, col] = prov_only[col].values

        # Compute anomalies
        for model_idx, obs_idx in ANOMALY_PAIRS:
            model_col = f"model_col{model_idx}"
            obs_col = f"obs_col{obs_idx}"
            anom_col = f"anom_TCRM{model_idx}_minus_OBS{obs_idx}"
            
            if model_col in merged.columns and obs_col in merged.columns:
                merged[anom_col] = merged[model_col] - merged[obs_col]
                print(f"✅ Added anomaly: {anom_col}")
            else:
                print(f"⚠️ Could not create anomaly for cols {model_idx}, {obs_idx}")

        # Drop normalized columns used for merging
        merged = merged.drop(columns=['Province_norm', 'Municipality_norm'], errors='ignore')

        # Save merged CSV
        rel = tcrm_file.relative_to(TCRM_ROOT)
        out_path = OUTPUT_ROOT / rel.parent / (tcrm_file.stem + OUTPUT_SUFFIX + ".csv")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        merged.to_csv(out_path, index=False)
        print(f"💾 Saved: {out_path.name}")
        created += 1

    except Exception as e:
        print(f"❌ Error in {tcrm_file.name}: {e}")
        import traceback
        traceback.print_exc()

# Summary
print("\n----- SUMMARY -----")
print(f"✅ Created merged files: {created}")
print(f"⚠️ Skipped (no OBS match): {skipped}")
print("Done.")