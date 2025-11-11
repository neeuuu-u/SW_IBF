# -*- coding: utf-8 -*-
"""
Created on Fri Oct 17 11:30:27 2025

@author: jo_ht
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 11:15:24 2025

@author: jo_ht

DATETIME RANGE VERSION - Extracts NetCDF timesteps within a specified datetime range
to use the script in command line, input []
"""

import netCDF4 as nc
import numpy as np
from pathlib import Path
import argparse
from netCDF4 import num2date
from datetime import datetime

def split_nc_by_time_range(input_file, output_dir, start_datetime, end_datetime):
    """
    Extract timesteps from a .nc file that fall within a specified datetime range.
    
    Args:
        input_file: Path to input .nc file
        output_dir: Directory to save output files
        start_datetime: Start datetime (datetime object or string in format 'YYYY-MM-DD HH:MM:SS')
        end_datetime: End datetime (datetime object or string in format 'YYYY-MM-DD HH:MM:SS')
    """
    # Convert to Path objects if they aren't already
    output_dir = Path(output_dir)
    input_file = Path(input_file)
    
    # Parse datetime strings if necessary
    if isinstance(start_datetime, str):
        start_datetime = datetime.strptime(start_datetime, '%Y-%m-%d %H:%M:%S')
    if isinstance(end_datetime, str):
        end_datetime = datetime.strptime(end_datetime, '%Y-%m-%d %H:%M:%S')
    
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
        
        # Filter timesteps within the specified range
        range_indices = []
        range_dates = []
        
        for t_idx, dt in enumerate(dates):
            if start_datetime <= dt <= end_datetime:
                range_indices.append(t_idx)
                range_dates.append(dt)
        
        print(f"Processing {input_file.name}:")
        print(f"  Time units: {time_units}")
        print(f"  Time calendar: {calendar}")
        print(f"  Total timesteps: {n_times}")
        print(f"  Requested range: {start_datetime} to {end_datetime}")
        print(f"  Timesteps in range: {len(range_indices)}")
        
        if len(range_indices) == 0:
            print(f"  WARNING: No timesteps found within the specified range in {input_file.name}")
            return
        
        # Show first and last timesteps in range
        print(f"  First timestep: {range_dates[0]}")
        print(f"  Last timestep: {range_dates[-1]}")
        
        # DEBUG: Print first 20 timesteps that will be created
        print(f"  First {min(20, len(range_indices))} timesteps to be created:")
        for i in range(min(20, len(range_indices))):
            t_idx = range_indices[i]
            dt = range_dates[i]
            print(f"    [{i}] Original index {t_idx}: {dt} (raw: {time_values[t_idx]:.2f})")
        print()
        
        # Process each time step in range
        for i, (t_idx, dt) in enumerate(zip(range_indices, range_dates)):
            # Read time value
            time_val = time_values[t_idx]
            
            # Create output filename with datetime
            input_stem = input_file.stem
            date_str = dt.strftime('%Y-%m-%d_%H%M%S')
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
                dst.setncattr('extraction_range', f"{start_datetime} to {end_datetime}")
            
            if (i + 1) % 50 == 0:
                print(f"  Processed {i + 1}/{len(range_indices)} timesteps")
        
        print(f"  Completed: {len(range_indices)} files created")

def process_directory(input_dir, output_base_dir, start_datetime, end_datetime, 
                     pattern="*.nc", subfolder_per_file=True):
    """
    Process all NetCDF files in a directory recursively, extracting timesteps within range.
    
    Args:
        input_dir: Input directory containing .nc files
        output_base_dir: Base output directory for split files
        start_datetime: Start datetime for extraction
        end_datetime: End datetime for extraction
        pattern: File pattern to match (default: *.nc)
        subfolder_per_file: If True, creates a subfolder for each file's outputs
    """
    input_path = Path(input_dir)
    output_base = Path(output_base_dir)
    
    # Find all .nc files recursively
    nc_files = list(input_path.rglob(pattern))
    
    if not nc_files:
        print(f"No .nc files found in {input_dir}")
        return
    
    print(f"Found {len(nc_files)} .nc files to process")
    print(f"Extracting timesteps from {start_datetime} to {end_datetime}\n")
    
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
            split_nc_by_time_range(nc_file, output_dir, start_datetime, end_datetime)
        except Exception as e:
            print(f"  ERROR: Failed to process {nc_file}: {e}")
        
        print()

def main():
    parser = argparse.ArgumentParser(
        description="Extract NetCDF timesteps within a specified datetime range",
        epilog="""
Examples:
  python %(prog)s input/ output/ "2024-10-01 00:00:00" "2024-10-03 23:59:59"
  python %(prog)s input/ output/ "2024-10-01 00:00:00" "2024-10-03 23:59:59" --pattern "wind_*.nc"
  python %(prog)s input/ output/ "2024-10-01 00:00:00" "2024-10-03 23:59:59" --no-subfolder
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        'input_dir',
        help='Input directory containing .nc files (will search recursively)'
    )
    parser.add_argument(
        'output_dir',
        help='Output directory for extracted files (will mirror input directory structure)'
    )
    parser.add_argument(
        'start_datetime',
        help='Start datetime in format "YYYY-MM-DD HH:MM:SS"'
    )
    parser.add_argument(
        'end_datetime',
        help='End datetime in format "YYYY-MM-DD HH:MM:SS"'
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
        args.start_datetime,
        args.end_datetime,
        args.pattern,
        subfolder_per_file=not args.no_subfolder
    )
    print("All files processed successfully!")

if __name__ == "__main__":
    import sys
    
    # If no command line arguments provided, use hardcoded settings
    if len(sys.argv) == 1:
        print("=" * 70)
        print("Running in HARDCODED MODE (no command line arguments provided)")
        print("Edit the settings below in the script to change paths/dates")
        print("=" * 70)
        print()
        
        # ========== EDIT THESE SETTINGS ==========
        INPUT_FOLDER = r"C:\Users\jo_ht\tcrm-3.1.16\output\opong_bst\windfield"
        OUTPUT_FOLDER = r"Z:\Severe_Wind\CIM oct28\SAR\OPONG_0925_21_32_31\BST"
        
        # Define your datetime range (format: 'YYYY-MM-DD HH:MM:SS')
        START_DATETIME = '2025-09-25 21:00:00'
        END_DATETIME = '2025-09-25 22:00:00'
        # =========================================
        
        # Process all .nc files recursively within the datetime range
        process_directory(INPUT_FOLDER, OUTPUT_FOLDER, START_DATETIME, END_DATETIME)
        print("All files processed successfully!")
    else:
        # Use command line arguments
        main()