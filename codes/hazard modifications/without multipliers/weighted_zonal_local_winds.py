## TCRM NC Gust Files to Municipal Weighted Mean Local Gusts
## Version 2025.05.06
## MHIBFEWS Project 

from osgeo import gdal, gdal_array
from scipy.stats import lognorm, mode
from pyproj import Transformer
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import csv, fnmatch, os, time, string, math
import pandas as pd
import shutil
from rasterstats import zonal_stats
from datetime import datetime as dt
from multiprocessing import Pool, cpu_count, set_start_method
import calendar
import glob
import re

from shapely.geometry import Polygon, MultiPolygon, shape, mapping


try:
    set_start_method("spawn", force=True)
except RuntimeError:
    pass

## Functions

def saveToTiff(outFile, arrayData, referenceData):
    if os.path.exists(outFile):
        os.remove(outFile)
        print(f'Deleted existing file: {outFile}')
    driver = gdal.GetDriverByName('GTiff')
    saveFile = driver.Create(outFile, referenceData.RasterXSize, referenceData.RasterYSize, 1, gdal.GDT_Float32,
                             options=['COMPRESS=LZW', 'TILED=YES'])
    saveFile.SetSpatialRef(referenceData.GetSpatialRef())
    saveFile.SetGeoTransform(referenceData.GetGeoTransform())
    saveFile.SetProjection(referenceData.GetProjection())
    saveFile.GetRasterBand(1).SetNoDataValue(-9999)
    saveFile.GetRasterBand(1).WriteArray(arrayData)
    saveFile.FlushCache()
    saveFile = None
    return


def resampleData(regionalGustData, multiplierData, outFile):
    sourceData = multiplierData
    xmin, xres, xskew, ymax, yskew, yres = sourceData.GetGeoTransform()
    xmax = xmin + (sourceData.RasterXSize * xres)
    ymin = ymax + (sourceData.RasterYSize * yres)
    resampledData = gdal.Warp(outFile, regionalGustData, format='GTiff', outputBounds=[xmin, ymin, xmax, ymax],
                              xRes=xres, yRes=-yres, dstSRS=sourceData.GetProjection(), resampleAlg='bilinear',
                              multithread=True, srcNodata=regionalGustData.GetRasterBand(1).GetNoDataValue(),
                              dstNodata=multiplierData.GetRasterBand(1).GetNoDataValue(),
                              creationOptions=['COMPRESS=LZW', 'TILED=YES'])
    return resampledData


def gdalDataToArray(gdalData):
    gdalArray = gdalData.ReadAsArray(0, 0, gdalData.RasterXSize, gdalData.RasterYSize)
    gdalBand = gdalData.GetRasterBand(1)
    gdalMin = gdalBand.GetMinimum()
    if gdalMin is None:
        gdalMin = 0
    gdalArray[gdalArray < gdalMin] = np.nan
    return gdalArray


def createDirNotExist(dirPath):
    if not os.path.isdir(dirPath):
        os.makedirs(dirPath)
    return dirPath


# def zonal_stats_worker(args):
#     row, pathtoLocGust, tif_files, prov_name = args
#     prov_name = row['clean_prov']
    
#     tif_file = next((f for f in tif_files if f.startswith(prov_name.lower())), None)
#     if tif_file:
#         tif_path = os.path.join(pathtoLocGust, tif_file)
#         stats = zonal_stats(row.geometry, tif_path, nodata=-9999, stats=["mean"])
#         mean_value = stats[0]['mean'] if stats[0]['mean'] is not None else 0
#     else:
#         mean_value = 0
    
#     return row.name, mean_value

def zonal_stats_worker(args):
    row, pathtoLocGust, tif_files, prov_name = args
    prov_name = row['clean_prov']
    
    tif_file = next((f for f in tif_files if f.startswith(prov_name.lower())), None)
    if not tif_file:
        return row.name, 0

    tif_path = os.path.join(pathtoLocGust, tif_file)
    geom = row.geometry

    if geom is None:
        print(f"[Warning] Invalid geometry in row {row.name}")
        return row.name, 0

    try:
        # Handle MultiPolygon by splitting and averaging results
        if isinstance(geom, MultiPolygon):
            means = []
            for poly in geom.geoms:
                stats = zonal_stats(mapping(poly), tif_path, nodata=-9999, stats=["mean"])
                if stats and stats[0]['mean'] is not None:
                    means.append(stats[0]['mean'])
            mean_value = sum(means) / len(means) if means else 0

        elif isinstance(geom, Polygon):
            stats = zonal_stats(mapping(geom), tif_path, nodata=-9999, stats=["mean"])
            mean_value = stats[0]['mean'] if stats and stats[0]['mean'] is not None else 0

        else:
            print(f"[Warning] Unsupported geometry type in row {row.name}")
            mean_value = 0

    except Exception as e:
        print(f"[Error] Zonal stats failed for row {row.name}: {e}")
        mean_value = 0

    return row.name, mean_value


def haversine(lat1, lon1, lat2, lon2):
    R = 6371.0
    lat1, lon1 = math.radians(lat1), math.radians(lon1)
    lat2, lon2 = math.radians(lat2), math.radians(lon2)
    dlat, dlon = lat2 - lat1, lon2 - lon1
    a = math.sin(dlat/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    return R * 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))


def get_provinces_within_buffer(eventName, ini, day_index, buffer_km, filename, pathtoHazard, metadataDir):
    YYYY = int(f"20{ini[4:6]}")
    track_csv_path = os.path.join(pathtoHazard, filename.replace(".nc", ".csv"))
    track_df = pd.read_csv(track_csv_path)

    provinces_in_buffer = []

    for file in os.listdir(metadataDir):
        if not file.endswith(".txt"):
            continue
        province_name = os.path.splitext(file)[0]
        file_path = os.path.join(metadataDir, file)

        metadata = {}
        with open(file_path, "r") as f:
            for line in f:
                if "=" in line:
                    key, value = map(str.strip, line.split("="))
                    try:
                        metadata[key] = float(value)
                    except ValueError:
                        continue

        if not all(k in metadata for k in ["xfirst", "xinc", "xsize", "yfirst", "yinc", "ysize"]):
            continue

        xfirst, xinc, xsize = metadata["xfirst"], metadata["xinc"], metadata["xsize"]
        yfirst, yinc, ysize = metadata["yfirst"], metadata["yinc"], metadata["ysize"]
        corners_proj = [
            (xfirst, yfirst),
            (xfirst + xinc * xsize, yfirst),
            (xfirst, yfirst + yinc * ysize),
            (xfirst + xinc * xsize, yfirst + yinc * ysize)
        ]

        transformer = Transformer.from_crs(32651, 4326, always_xy=True)
        corners_latlon = [transformer.transform(x, y)[::-1] for x, y in corners_proj]

        found = False
        for _, row in track_df.iterrows():
            tc_lat, tc_lon = row["LAT"], row["LON"]
            for corner_lat, corner_lon in corners_latlon:
                if haversine(tc_lat, tc_lon, corner_lat, corner_lon) <= buffer_km:
                    provinces_in_buffer.append(province_name)
                    found = True
                    break
            if found:
                break

    return sorted(set(provinces_in_buffer))

def regGustToLocGustMean(eventName, ini, day_index, new_year, mulFiles, filename, regGustData, gdf, brgy_gdf, pathtoMulti, pathtoHazard):

    # Filter provinces to process
    metadataDir = os.path.join(pathtoMulti, 'metadata')
    buffer_km = 200  # or make it a function argument
    provinces_to_process = get_provinces_within_buffer(eventName, ini, day_index, buffer_km, filename, pathtoHazard, metadataDir)
    print(f"Filtered provinces within buffer: {provinces_to_process}")

    pathtoLocGust = os.path.join(os.path.dirname(pathtoHazard),f'{eventName}_{new_year}', f'{eventName}_{ini}', f"{filename.split('.')[0]}_localgust")
    createDirNotExist(pathtoLocGust)
    pathtoResampled = os.path.join(os.path.dirname(pathtoHazard), f'{eventName}_{new_year}', f'{eventName}_{ini}', f"{filename.split('.')[0]}_resampled")
    createDirNotExist(pathtoResampled)
    
    # Limit multiplier files to only those matching filtered provinces
    mulFiles = [f for f in mulFiles if f.split("_")[0] in provinces_to_process]
    for multiplierFile in mulFiles:
        provinceName = multiplierFile.split("_")[0]
        mulData = gdal.Open(os.path.join(pathtoMulti, multiplierFile))
        multiplierArray = gdalDataToArray(mulData)

        print(f'Resampling {provinceName}...')
        regGustResampled = resampleData(regGustData, mulData, os.path.join(pathtoResampled, f'{provinceName}_resampled.tif'))
        regGustArray = gdalDataToArray(regGustResampled)
        del regGustResampled
        
        regGustkphArray = np.multiply(regGustArray, 3.6)
        regGustkphArray[np.isnan(multiplierArray)] = np.nan
        locGustkphArray = np.multiply(multiplierArray, regGustkphArray)
        del regGustkphArray
            
        saveToTiff(os.path.join(pathtoLocGust, f'{provinceName}_localgust.tif'), locGustkphArray, mulData)

    # Computing the average local gust and weighing it based on building population
    if brgy_gdf is not None:
        raster_files = os.listdir(pathtoLocGust)
        tif_files = [f for f in raster_files if f.lower().endswith('.tif')]
        
        toc = time.time()
        print("\nProcessing zonalstats...")
        
        results = []    

        for idx, row in brgy_gdf.iterrows():
            result = zonal_stats_worker((row, pathtoLocGust, tif_files, row['clean_prov']))
            results.append(result)

        for idx, mean_value in results:
            if brgy_gdf.at[idx, 'clean_prov'] in provinces_to_process:
                brgy_gdf.at[idx, 'mean_ctrl'] = mean_value
            else:
                brgy_gdf.at[idx, 'mean_ctrl'] = np.nan
        
        print(f'\nCompleted zonalstats after {"%.2f" % ((time.time() - toc) / 60)} minutes.\n')

    brgy_gdf['wtd_mean'] = brgy_gdf['mean_ctrl']
    brgy_gdf_filtered = brgy_gdf[brgy_gdf['clean_prov'].isin(provinces_to_process)]

    gdf['mean_ctrl'] = gdf['Mun_Name'].map(brgy_gdf_filtered.groupby('Mun_Name')['mean_ctrl'].mean())
    gdf['wtd_mean'] = gdf['Mun_Name'].map(brgy_gdf_filtered.groupby('Mun_Name')['wtd_mean'].sum())
    gdf.loc[~gdf['clean_prov'].isin(provinces_to_process), ['mean_ctrl', 'wtd_mean']] = np.nan
    
    shutil.rmtree(pathtoLocGust, ignore_errors=True)
    shutil.rmtree(pathtoResampled, ignore_errors=True)
    del brgy_gdf
    return gdf


def processInitialization(eventName, ini, fcast_days, pathtoHazard, pathtoMulti, mun_gdf, brgy_gdf, cluster_thresholds, filename_format):

    tic = time.time()
    
    print(f"\nSimulation initialized with data: {eventName}_{ini}")
    
    YYYY = int(f"20{ini[4:6]}")
    MM = int(ini[:2])
    DD = int(ini[2:4]) 
    hh = ini.split("_")[1][:2]
    mm = ini.split("_")[1][2:]

    start_day = 1 if f"{hh:02}{mm:02}" != "0000" else 0
    last_day = fcast_days + 1 if f"{hh:02}{mm:02}" != "0000" else fcast_days
    days_in_month = calendar.monthrange(YYYY, MM)[1]
    
    mulFiles = [os.path.basename(f) for f in glob.glob(f"{pathtoMulti}/*.tif")]

    for day in range(start_day, last_day):
        print(f"\nDay {day}")
        total_days = DD + day

        # Check if new_day exceeds the last day of the month
        if total_days > days_in_month:
            new_day = total_days - days_in_month 
            new_month = MM + 1  
            if new_month > 12:
                new_month = 1
                new_year = YYYY + 1
            month = dt(1900, new_month, 1).strftime('%b')
        else:
            new_day = DD + day
            new_month = MM
            new_year = YYYY
            month = dt(1900, MM, 1).strftime('%b')
            
        print(f"Processing forecast for {month} {new_day}, {new_year}...")
        pathtoProcessed = os.path.join(os.path.dirname(pathtoHazard),f'{eventName}_{new_year}', f'{eventName}_{ini}')
        createDirNotExist(pathtoProcessed)
        output_path = os.path.join(pathtoProcessed, f'{eventName}_{ini}_{month}{new_day}.csv')
        
        # Process the control track
        gdal.DontUseExceptions()
        filename = filename_format.format(eventName=eventName, ini=ini, day=day)
        print(f"Loading {os.path.join(pathtoHazard, filename)}...\n")
        regGustData = gdal.Open(f"NETCDF:{os.path.join(pathtoHazard, filename)}:vmax")

        # Initialize day_gdf once per day
        bldgs_sum_per_mun = brgy_gdf.groupby('Mun_Name')['bldgs'].sum()
        brgy_gdf['norm_bldgs'] = brgy_gdf.apply(lambda row: row['bldgs'] / bldgs_sum_per_mun[row['Mun_Name']], axis=1)
        brgy_gdf['mean_ctrl'] = np.nan
        brgy_gdf['wtd_mean'] = np.nan
        
        mun_gdf['cluster'] = mun_gdf['Mun_Name'].apply(lambda mun_name: mun_name if mun_name in cluster_thresholds else 'Generalized')
        mun_gdf['mean_ctrl'] = np.nan
        mun_gdf['wtd_mean'] = np.nan
        mun_gdf = regGustToLocGustMean(eventName, ini, day, new_year, mulFiles, filename, regGustData,
            mun_gdf, brgy_gdf, pathtoMulti, pathtoHazard)
        print(mun_gdf.head())
        
        mun_df = pd.DataFrame(mun_gdf)
        mun_df = mun_gdf.drop(columns=["geometry"])
        mun_df.to_csv(output_path, index=False)
        print(f"\nOutput CSV created: {output_path}")

        del mun_df
        
    print(f'\nProcessed in {"%.2f" % ((time.time() - tic) / 60)} minutes.')


############### For Single Swath
def processInitialization_single_swath(eventName, ini, fcast_days, pathtoHazard, pathtoMulti, mun_gdf, brgy_gdf, cluster_thresholds, filename_format):

    tic = time.time()
    
    print(f"\nSimulation initialized with data: {eventName}_{ini}")
    
    YYYY = int(f"20{ini[4:6]}")
    day = ''
    mulFiles = [os.path.basename(f) for f in glob.glob(f"{pathtoMulti}/*.tif")]

    pathtoProcessed = os.path.join(os.path.dirname(pathtoHazard),f'{eventName}_{YYYY}', f'{eventName}_{ini}')
    createDirNotExist(pathtoProcessed)
    output_path = os.path.join(pathtoProcessed, f'{eventName}_{ini}.csv')
    
    # Process the control track
    gdal.DontUseExceptions()
    filename = filename_format.format(eventName=eventName, ini=ini)
    print(f"Loading {os.path.join(pathtoHazard, filename)}...\n")
    regGustData = gdal.Open(f"NETCDF:{os.path.join(pathtoHazard, filename)}:vmax")

    # Initialize day_gdf once per day
    bldgs_sum_per_mun = brgy_gdf.groupby('Mun_Name')['bldgs'].sum()
    brgy_gdf['norm_bldgs'] = brgy_gdf.apply(lambda row: row['bldgs'] / bldgs_sum_per_mun[row['Mun_Name']], axis=1)
    brgy_gdf['mean_ctrl'] = np.nan
    brgy_gdf['wtd_mean'] = np.nan
    
    mun_gdf['cluster'] = mun_gdf['Mun_Name'].apply(lambda mun_name: mun_name if mun_name in cluster_thresholds else 'Generalized')
    mun_gdf['mean_ctrl'] = np.nan
    mun_gdf['wtd_mean'] = np.nan
    mun_gdf = regGustToLocGustMean(eventName, ini, day, YYYY, mulFiles, filename, regGustData,
        mun_gdf, brgy_gdf, pathtoMulti, pathtoHazard)
    print(mun_gdf.head())
    
    mun_df = pd.DataFrame(mun_gdf)
    mun_df = mun_gdf.drop(columns=["geometry"])
    mun_df.to_csv(output_path, index=False)
    print(f"\nOutput CSV created: {output_path}")

    del mun_df
        
    print(f'\nProcessed in {"%.2f" % ((time.time() - tic) / 60)} minutes.')
