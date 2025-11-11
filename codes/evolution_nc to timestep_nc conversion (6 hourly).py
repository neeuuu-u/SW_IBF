# -*- coding: utf-8 -*-
"""
Created on Fri Oct 17 11:29:42 2025

@author: jo_ht
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 11:15:24 2025

@author: jo_ht

6-HOURLY VERSION - Splits NetCDF files into 6-hourly timesteps (00, 06, 12, 18)
"""

import netCDF4 as nc
import numpy as np
from pathlib import Path
import argparse
from netCDF4 import num2date

def split_nc_by_time(input_file, output_dir, preserve_structure=True):
    """
    Split a single .nc file into multiple files, one per 6-hourly time step.
    Only keeps timesteps that fall exactly on 0, 6, 12, 18 hours.
    
    Args:
        input_file: Path to input .nc file
        output_dir: Directory to save output files
        preserve_structure: If True, creates subfolder for each file's outputs
    """
    # Convert to Path objects if they aren't already
    output_dir = Path(output_dir)
    input_file = Path(input_file)
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Open the input file
    with nc.Dataset(input_file, 'r') as src:
        # Get dimensions
        n_times = len(src.dimensions['time'])
        n_lat = len(src.dimensions['lat'])
        n_lon = len(src.dimensions['lon'])
        
        # Read coordinate variables
        lat = src.variables['lat'][:]
        lon = src.variables['lon'][:]
        time_var = src.variables['time']
        time_units = time_var.units
        time_values = time_var[:]
        
        # Convert time values to datetime objects
        try:
            calendar = getattr(time_var, 'calendar', 'standard')
            dates = num2date(time_values, units=time_units, calendar=calendar)
        except Exception as e:
            print(f"  ERROR: Could not parse time units: {e}")
            print(f"  Time units: {time_units}")
            return
        
        # Filter to get only 6-hourly timesteps (hour divisible by 6, minute=0, second=0)
        hourly_indices = []
        hourly_dates = []
        
        for t_idx, dt in enumerate(dates):
            if dt.hour % 6 == 0 and dt.minute == 0 and dt.second == 0:
                hourly_indices.append(t_idx)
                hourly_dates.append(dt)
        
        print(f"Processing {input_file.name}:")
        print(f"  Time units: {time_units}")
        print(f"  Time calendar: {calendar}")
        print(f"  Total timesteps: {n_times}")
        print(f"  6-hourly timesteps: {len(hourly_indices)}")
        
        if len(hourly_indices) == 0:
            print(f"  WARNING: No 6-hourly timesteps found in {input_file.name}")
            return
        
        # DEBUG: Print first 20 6-hourly timesteps that will be created
        print(f"  First {min(20, len(hourly_indices))} 6-hourly timesteps to be created:")
        for i in range(min(20, len(hourly_indices))):
            t_idx = hourly_indices[i]
            dt = hourly_dates[i]
            print(f"    [{i}] Original index {t_idx}: {dt} (raw: {time_values[t_idx]:.2f})")
        
        # Show time intervals if more than 1 timestep
        if len(hourly_indices) > 1:
            print(f"  Time intervals between 6-hourly timesteps:")
            for i in range(min(5, len(hourly_indices)-1)):
                delta = hourly_dates[i+1] - hourly_dates[i]
                print(f"    {hourly_dates[i]} -> {hourly_dates[i+1]} = {delta}")
        print()
        
        # Process each 6-hourly time step
        for i, (t_idx, dt) in enumerate(zip(hourly_indices, hourly_dates)):
            # Read time value
            time_val = time_values[t_idx]
            
            # Create output filename with datetime
            input_stem = input_file.stem
            date_str = dt.strftime('%Y-%m-%d_%H%M')
            output_file = output_dir / f"{input_stem}_{date_str}.nc"
            
            # Create new NetCDF file for this time step
            with nc.Dataset(output_file, 'w', format='NETCDF4') as dst:
                # Create dimensions
                dst.createDimension('time', 1)
                dst.createDimension('lat', n_lat)
                dst.createDimension('lon', n_lon)
                
                # Create coordinate variables
                time_var_out = dst.createVariable('time', 'f8', ('time',))
                time_var_out.units = time_units
                time_var_out[:] = time_val
                
                lat_var = dst.createVariable('lat', 'f4', ('lat',))
                lat_var[:] = lat
                
                lon_var = dst.createVariable('lon', 'f4', ('lon',))
                lon_var[:] = lon
                
                # Create and populate data variables
                for var_name in ['pressure', 'velocity_east', 'velocity_north', 'gust_speed']:
                    src_var = src.variables[var_name]
                    
                    # Create variable with compression
                    dst_var = dst.createVariable(
                        var_name, 
                        'f4', 
                        ('time', 'lat', 'lon'),
                        fill_value=src_var._FillValue,
                        zlib=True,
                        complevel=4
                    )
                    
                    # Copy data for this time step
                    dst_var[0, :, :] = src_var[t_idx, :, :]
                
                # Copy global attributes
                dst.setncattr('description', src.getncattr('description'))
                dst.setncattr('source_file', str(input_file.name))
                dst.setncattr('time_index', t_idx)
                dst.setncattr('datetime', date_str)
            
            if (i + 1) % 50 == 0:
                print(f"  Processed {i + 1}/{len(hourly_indices)} 6-hourly timesteps")
        
        print(f"  Completed: {len(hourly_indices)} files created")

def process_directory(input_dir, output_base_dir, pattern="*.nc", subfolder_per_file=True):
    """
    Process all NetCDF files in a directory recursively.
    
    Args:
        input_dir: Input directory containing .nc files
        output_base_dir: Base output directory for split files
        pattern: File pattern to match (default: *.nc)
        subfolder_per_file: If True, creates a subfolder for each file's outputs
                           If False, outputs directly to the same directory level
    """
    input_path = Path(input_dir)
    output_base = Path(output_base_dir)
    
    # Find all .nc files recursively
    nc_files = list(input_path.rglob(pattern))
    
    if not nc_files:
        print(f"No .nc files found in {input_dir}")
        return
    
    print(f"Found {len(nc_files)} .nc files to process\n")
    
    for i, nc_file in enumerate(nc_files, 1):
        print(f"[{i}/{len(nc_files)}] Processing {nc_file}")
        
        # Create output directory preserving relative structure
        rel_path = nc_file.parent.relative_to(input_path)
        
        if subfolder_per_file:
            # Create subfolder for each file: output/subdir/filename/
            output_dir = output_base / rel_path / nc_file.stem
        else:
            # Keep same directory level: output/subdir/
            output_dir = output_base / rel_path
        
        try:
            split_nc_by_time(nc_file, output_dir)
        except Exception as e:
            print(f"  ERROR: Failed to process {nc_file}: {e}")
        
        print()

def main():
    parser = argparse.ArgumentParser(
        description="Split NetCDF files by 6-hourly time steps (keeps only 00, 06, 12, 18 hour timestamps)"
    )
    parser.add_argument(
        'input_dir',
        help='Input directory containing .nc files (will search recursively)'
    )
    parser.add_argument(
        'output_dir',
        help='Output directory for split files (will mirror input directory structure)'
    )
    parser.add_argument(
        '--pattern',
        default='*.nc',
        help='File pattern to match (default: *.nc)'
    )
    parser.add_argument(
        '--no-subfolder',
        action='store_true',
        help='Output files directly to same directory level instead of creating subfolders per file'
    )
    
    args = parser.parse_args()
    
    process_directory(
        args.input_dir, 
        args.output_dir, 
        args.pattern,
        subfolder_per_file=not args.no_subfolder
    )
    print("All files processed successfully!")

if __name__ == "__main__":
    # ========== EDIT THESE PATHS ==========
    INPUT_FOLDER = r"C:\Users\jo_ht\tcrm-3.1.16\output\UWAN_110825_0000\windfield\evolution"
    OUTPUT_FOLDER = r"C:\Users\jo_ht\tcrm-3.1.16\output\UWAN_110825_0000\windfield\gust\instantaneous"
    # ======================================
    
    # Process all .nc files recursively
    process_directory(INPUT_FOLDER, OUTPUT_FOLDER)
    print("All files processed successfully!")