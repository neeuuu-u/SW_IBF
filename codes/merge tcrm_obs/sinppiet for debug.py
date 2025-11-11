# -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 18:40:13 2025

@author: jo_ht
"""

import pandas as pd

f = r"D:\neu\CIM oct28\zones [final]\evolution.001-00001_2025-09-26_0000_multi_zonal_stats_merged_obs.csv"
df = pd.read_csv(f)
print(df.head())
print(df.columns)
