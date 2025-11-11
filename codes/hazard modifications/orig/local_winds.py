# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 14:05:32 2025

SWIFTPh Hazard Module
Local winds
Input: TCRM Regional Wind Gust (.nc)
       TCRM-formatted TC tracks (.csv)
       Wind multipliers (provincial level) (.tiff)
       Philippine municipal & brgy boundaries
Output: Municipal Weighted Mean Local Gusts (.csv)
"""
import time, os, re, sys, shutil
import geopandas as gpd
os.chdir(os.path.dirname(os.path.abspath(__file__)))
from weighted_zonal_local_winds import processInitialization, processInitialization_single_swath

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from swiftph.prep.prep import load_config, update_tcinfo_config, extract_metadata
from swiftph.prep.create_tcrm_input_files_ini_csv import generate_date_ranges
import swiftph.prep.utils as utils

if __name__ == '__main__':
    update_tcinfo_config()
    config = load_config() 
    file_path = os.path.join(config['swiftphDir'], 'data', 'tc', f"{config['eventName']}_{config['ini']}.xlsx")
    if not os.path.isfile(file_path):
        print("TC Excel data file not found!")
        sys.exit(1)

    tc_name, year, initialization, filename = extract_metadata(file_path, prnt= False)
    
    #Load TC Info
    eventName = config['eventName']
    ini = config['ini']
    fcast_days = int(config['fcast_days'])
    
    #Define paths
    swiftphDir = config['swiftphDir']
    swiftph_outDir = os.path.join(swiftphDir,'output')
    tcrm_input_dataDir = os.path.join(config['tcrmDir'], 'input', 'tc')
    tcrm_output_dataDir = os.path.join(config['tcrmDir'], 'output')
    
    pathtoMulti = os.path.join(swiftphDir, 'data', 'multiplier')
    pathtoProvSHP = os.path.join(swiftphDir, 'data' , 'shp' , 'ph_admin_prov_boundaries', 'PHL_adm4_PSA_pn_2016Junprj_prov.shp')
    pathtoMunSHP = os.path.join(swiftphDir, 'data' , 'shp' , 'ph_admin_mun_boundaries', 'PHL_adm4_PSA_pn_2016Junprj_mun.shp')
    pathtoBrgySHP = os.path.join(swiftphDir, 'data' , 'shp' , 'ph_admin_brgy_boundaries', 'PHL_adm4_PSA_pn_2016Junprj_geomfix_bldgs_11262024.shp')
    
    # Create subdirectory for the TC event in swiftph output directory
    output_subdir = os.path.join(swiftph_outDir, f'{tc_name}_{year}')
    os.makedirs(output_subdir, exist_ok=True)
    
    pathtoHazard = output_subdir

    #Locate Regional Wind Gusts data
    if config['mode'] == 'single':
        csv_file = filename + '.csv'
        nc_file = filename + '.nc'
        
        # Source CSV path
        src_csv_path = os.path.join(tcrm_input_dataDir, csv_file)
        dst_csv_path = os.path.join(output_subdir, csv_file)
        
        # Copy CSV file
        if os.path.exists(src_csv_path):
            shutil.copy2(src_csv_path, dst_csv_path)
        else:
            print(f"[WARNING] CSV file not found: {src_csv_path}")

        # Locate windfield folder
        windfield_dir = os.path.join(tcrm_output_dataDir, filename, 'windfield')
        
        if os.path.exists(windfield_dir):
            # Find the only NetCDF file in the windfield folder
            nc_files = [f for f in os.listdir(windfield_dir) if f.endswith('.nc')]
            if len(nc_files) == 1:
                original_nc_path = os.path.join(windfield_dir, nc_files[0])
                renamed_nc_path = os.path.join(windfield_dir, nc_file)
        
                # Rename the .nc file
                os.rename(original_nc_path, renamed_nc_path)
        
                # Copy to output subdirectory
                dst_nc_path = os.path.join(output_subdir, nc_file)
                shutil.copy2(renamed_nc_path, dst_nc_path)
            else:
                print(f"[WARNING] Expected 1 .nc file in {windfield_dir}, found {len(nc_files)}.")
        else:
            print(f"[WARNING] Windfield folder not found: {windfield_dir}")
        
        tic = time.time()
        
        
        province_shapes = {}
        for meta_file in os.listdir(os.path.join(pathtoMulti, 'metadata')):
            if meta_file.endswith(".txt"):
                province = meta_file.replace(".txt", "")
                meta_path = os.path.join(pathtoMulti, 'metadata', meta_file)
        
                with open(meta_path, 'r') as f:
                    content = f.read()
                    ysize = int(re.search(r'ysize\s*=\s*(\d+)', content).group(1))
                    xsize = int(re.search(r'xsize\s*=\s*(\d+)', content).group(1))
                    province_shapes[province] = (ysize, xsize)
        
        mun_gdf = gpd.read_file(pathtoMunSHP)
        mun_gdf["clean_prov"] = mun_gdf["Pro_Name"].apply(utils.clean_prov_name)
        print(f"Successfully read and cleaned {pathtoMunSHP}.\n")
        print(mun_gdf.head())
        
        brgy_gdf = gpd.read_file(pathtoBrgySHP)
        brgy_gdf['clean_prov'] = brgy_gdf['Pro_Name'].apply(utils.clean_prov_name)
        print(f"Successfully read and cleaned {pathtoBrgySHP}.\n")
        print(brgy_gdf.head())
        
        # Processing TCRM NC Gust Files to Municipal Weighted Mean Local Gusts
        processInitialization_single_swath(eventName, ini, fcast_days, pathtoHazard, pathtoMulti, mun_gdf, brgy_gdf,
                              utils.cluster_thresholds, filename_format="{eventName}_{ini}.nc")
    

    
    elif config['mode'] == 'time-segmented':
        date_ranges = generate_date_ranges(initialization, num_days=fcast_days)   
        for i, (start_date, end_date) in enumerate(date_ranges, start=0):
            csv_file = filename + f'_Day{i}.csv'
            nc_file = filename + f'_Day{i}.nc'
        
            # Source CSV path
            src_csv_path = os.path.join(tcrm_input_dataDir, csv_file)
            dst_csv_path = os.path.join(output_subdir, csv_file)
        
            # Copy CSV file
            if os.path.exists(src_csv_path):
                shutil.copy2(src_csv_path, dst_csv_path)
                print(f'Moved CSV to {dst_csv_path}')
            else:
                print(f"[WARNING] CSV file not found: {src_csv_path}")
        
            # Locate windfield folder
            windfield_dir = os.path.join(tcrm_output_dataDir, filename + f'_Day{i}', 'windfield')
        
            if os.path.exists(windfield_dir):
                # Find the only NetCDF file in the windfield folder
                nc_files = [f for f in os.listdir(windfield_dir) if f.endswith('.nc')]
                if len(nc_files) == 1:
                    original_nc_path = os.path.join(windfield_dir, nc_files[0])
                    renamed_nc_path = os.path.join(windfield_dir, nc_file)
        
                    # Rename the .nc file
                    os.rename(original_nc_path, renamed_nc_path)
        
                    # Copy to output subdirectory
                    dst_nc_path = os.path.join(output_subdir, nc_file)
                    shutil.copy2(renamed_nc_path, dst_nc_path)
                    print(f'Moved NC to {dst_nc_path}')
                else:
                    print(f"[WARNING] Expected 1 .nc file in {windfield_dir}, found {len(nc_files)}.")
            else:
                print(f"[WARNING] Windfield folder not found: {windfield_dir}")
            
            tic = time.time()
            
            
        province_shapes = {}
        for meta_file in os.listdir(os.path.join(pathtoMulti, 'metadata')):
            if meta_file.endswith(".txt"):
                province = meta_file.replace(".txt", "")
                meta_path = os.path.join(pathtoMulti, 'metadata', meta_file)
        
                with open(meta_path, 'r') as f:
                    content = f.read()
                    ysize = int(re.search(r'ysize\s*=\s*(\d+)', content).group(1))
                    xsize = int(re.search(r'xsize\s*=\s*(\d+)', content).group(1))
                    province_shapes[province] = (ysize, xsize)
        
        mun_gdf = gpd.read_file(pathtoMunSHP)
        mun_gdf["clean_prov"] = mun_gdf["Pro_Name"].apply(utils.clean_prov_name)
        print(f"Successfully read and cleaned {pathtoMunSHP}.\n")
        print(mun_gdf.head())
        
        brgy_gdf = gpd.read_file(pathtoBrgySHP)
        brgy_gdf['clean_prov'] = brgy_gdf['Pro_Name'].apply(utils.clean_prov_name)
        print(f"Successfully read and cleaned {pathtoBrgySHP}.\n")
        print(brgy_gdf.head())
        
        # Processing TCRM NC Gust Files to Municipal Weighted Mean Local Gusts
        processInitialization(eventName, ini, fcast_days, pathtoHazard, pathtoMulti, mun_gdf, brgy_gdf,
                              utils.cluster_thresholds, filename_format="{eventName}_{ini}_Day{day}.nc")
    

    