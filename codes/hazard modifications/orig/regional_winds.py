# -*- coding: utf-8 -*-
"""
Created on Thu May  1 14:45:01 2025

SWIFTPh Hazard Module
Regional winds
Input: TCRM-formatted Files (.ini and .csv)
Output: TCRM Regional Wind Gust (.nc)
"""

import os
import sys
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import shutil, time
import subprocess
import multiprocessing

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from swiftph.prep.prep import load_config

def single_tc_swath(config):
    pagasadata_file = config['pagasadata_file']
    if pagasadata_file is not None:
        base_filename = pagasadata_file[:-5]
        tcrm_work_dir = config['tcrmDir']
        print(f'TC Single Swath Simulation for {base_filename}')
        
        base_file_folder = os.path.join(tcrm_work_dir, 'output', base_filename)
        if os.path.exists(base_file_folder):
            shutil.rmtree(base_file_folder)
        
        ini_file = os.path.join(tcrm_work_dir, 'input', 'tc', f'{base_filename}.ini')
        tcevent_script = os.path.join(tcrm_work_dir, 'tcevent.py')
    
        if not os.path.exists(ini_file):
            print(f"Error: {ini_file} does not exist.\n")
            return
    
        if not os.path.exists(tcevent_script):
            print(f"Error: {tcevent_script} does not exist.\n")
            return
    
        # conda_activate = config['condaDir'] + "\\Scripts\\activate.bat"
        # activate_conda = f'call {conda_activate} tcrm'
        # command = f'{activate_conda} && python {tcevent_script} -c {ini_file} -v'
        command = f'python {tcevent_script} -c {ini_file} -v'
        subprocess.run(command, shell=True)
    else:
        print('Invalid TC Data or TCRM input files missing')
# def run_single_forecast_day(args):
#     base_filename, config = args
#     tcrm_work_dir = config['tcrmDir']
#     tcevent_script = os.path.join(tcrm_work_dir, 'tcevent.py')
#     conda_activate = config['condaDir'] + "\\Scripts\\activate.bat"
#     activate_conda = f'call {conda_activate} tcrm'
    
#     base_file_folder = os.path.join(tcrm_work_dir, 'output', base_filename)
#     if os.path.exists(base_file_folder):
#         shutil.rmtree(base_file_folder)

#     ini_file = os.path.join(tcrm_work_dir, 'input', 'tc', f'{base_filename}.ini')
#     if not os.path.exists(ini_file):
#         print(f"Warning: {ini_file} does not exist. Skipping...\n")
#         return

#     command = f'{activate_conda} && python {tcevent_script} -c {ini_file} -v'
#     subprocess.run(command, shell=True)

def run_single_forecast_day(args):
    """
    Runs a single forecast day process, compatible with both Windows and Linux.

    Args:
        args (tuple): A tuple containing:
            - base_filename (str): The base filename for the forecast.
            - config (dict): A dictionary containing configuration details,
                             including 'tcrmDir' and 'condaDir'.
    """
    base_filename, config = args
    tcrm_work_dir = config['tcrmDir']
    tcevent_script = os.path.join(tcrm_work_dir, 'tcevent.py')

    base_file_folder = os.path.join(tcrm_work_dir, 'output', base_filename)
    if os.path.exists(base_file_folder):
        shutil.rmtree(base_file_folder)

    ini_file = os.path.join(tcrm_work_dir, 'input', 'tc', f'{base_filename}.ini')
    if not os.path.exists(ini_file):
        print(f"Warning: {ini_file} does not exist. Skipping...\n")
        return

    if sys.platform.startswith('win'):
        # Windows specific commands
        # conda_activate = os.path.join(config['condaDir'], 'Scripts', 'activate.bat')
        # activate_conda = f'call "{conda_activate}" tcrm'
        # command = f'{activate_conda} && python "{tcevent_script}" -c "{ini_file}" -v'
        command = f'python {tcevent_script} -c {ini_file} -v'
        subprocess.run(command, shell=True)
    elif sys.platform.startswith('linux') or sys.platform.startswith('darwin'):
        # Linux and macOS specific commands (macOS uses 'darwin' for sys.platform)
        conda_activate = os.path.join(config['condaDir'], 'activate')
        # Use single quotes around the command to ensure proper execution in bash
        command = f"bash -c 'source \"{conda_activate}\" tcrm && python \"{tcevent_script}\" -c \"{ini_file}\" -v'"
        subprocess.run(command, shell=True, executable="/bin/bash")
    else:
        print(f"Error: Unsupported operating system: {sys.platform}")
        return

def time_segmented_swath(config):
    pagasadata_file = config['pagasadata_file']
    if pagasadata_file is None:
        print('Please select TC data and create its TCRM input files')
        return

    base_filename_base = pagasadata_file[:-5]
    print(f'Time-Segmented TC Swaths Simulation for {base_filename_base}')

    tasks = [(f"{base_filename_base}_Day{n}", config) for n in range(5)]

    with multiprocessing.Pool(processes=min(5, multiprocessing.cpu_count())) as pool:
        pool.map(run_single_forecast_day, tasks)

def run_command(command):
    def execute_command():
        print("\nModel Simulation finished.\n")
        pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
        if config['mode'] == 'single':
            pool.map(single_tc_swath)
        elif  config['mode'] == 'time-segmented':
            pool.map(time_segmented_swath)
        pool.close()
        pool.join()   
        
if __name__ == '__main__':
    config = load_config()   
    start_time = time.time()
    
    if config['mode'] == 'single':
        single_tc_swath(config)
    elif config['mode'] == 'time-segmented':
        time_segmented_swath(config)
    else:
        print("Invalid mode of simulation (Options: 'single' or 'time-segmented').")
        
    
    end_time = time.time()
    duration = end_time - start_time
    minutes, seconds = divmod(duration, 60)
    print(f"Total time: {int(minutes)} min {seconds:.2f} sec\n")