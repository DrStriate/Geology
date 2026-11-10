import pandas as pd
import numpy as np
import pmagpy.pmag as pmag
import os
from pathlib import Path
import arcpy
# Python.exe needs to be in ArcGIS:  C:\Program Files\ArcGIS\Pro\bin\Python\envs\arcgispro-py3\python.exe


def parse_magic_file(filepath):
    """Parses multi-section MagIC v3 files into a dictionary of DataFrames."""
    tables = {}
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    sections = content.split('>>>>>>>>>>')
    for section in sections:
        lines = [line.strip() for line in section.strip().split('\n') if line.strip()]
        if not lines:
            continue
        
        table_name = None
        data_start_idx = 0
        for i, line in enumerate(lines):
            if line.startswith('tab delimited'):
                table_name = line.split('tab delimited')[1].strip()
                data_start_idx = i + 1
                break
        
        if table_name and data_start_idx < len(lines):
            header = lines[data_start_idx].split('\t')                  
            rows = [line.split('\t') for line in lines[data_start_idx+1:]]
            valid_rows = [r for r in rows if len(r) == len(header)]
            if valid_rows:
                tables[table_name] = pd.DataFrame(valid_rows, columns=header)
    return tables

# 0. Get the list of MagIC files
magICPath = Path.cwd() / "MagIC"
fileNames = os.listdir(magICPath)
site_results = []
header_sent = False
cleaned_csv = os.path.abspath("cleaned_data.csv")
for fileName in fileNames:
    # 1. Load the MagIC file
    magic_tables = parse_magic_file(magICPath / fileName)
    df = magic_tables['locations']

    # 2. Clean numerical fields
    df['lat_n'] = pd.to_numeric(df['lat_n'], errors='coerce')
    df['lon_w'] = pd.to_numeric(df['lon_w'], errors='coerce')
    df['age_high'] = pd.to_numeric(df['age_high'], errors='coerce')
    df['age_low'] = pd.to_numeric(df['age_low'], errors='coerce')
    df['dir_dec'] = pd.to_numeric(df['dir_dec'], errors='coerce')
    df['dir_inc'] = pd.to_numeric(df['dir_inc'], errors='coerce')
    df['age_av'] = (df['age_low'] + df['age_high']) / 2.0
    df['dir_tilt_correction'] = pd.to_numeric(df['dir_tilt_correction'], errors='coerce')

    # 3. Filter for fully tilt-corrected data (100% structural restoration)[cite: 1, 2]
    df_corrected = df[df['dir_tilt_correction'] == 100].copy()

    # Convert MagIC 0-360 longitude to standard -180 to 180 for ArcGIS
    df_corrected['x_lon'] = df['lon_w'].apply(lambda lon: lon - 360 if lon > 180 else lon)
    df_corrected['y_lat'] = df['lat_n']

    # 4. Handle polarity reversals (flip negative inclinations and their declinations by 180°)
    def normalize_polarity(row):
        dec = row['dir_dec']
        inc = row['dir_inc']
        if inc < 0:
            return (dec + 180) % 360, -inc
        return dec, inc

    df_corrected[['norm_dec', 'norm_inc']] = df_corrected.apply(
        lambda row: pd.Series(normalize_polarity(row)), axis=1
    )

    # 5. Group by geological location/site and compute Fisher Mean
    for location, group in df_corrected.groupby('location'):
        # Prepare list of [declination, inclination] pairs for PmagPy
        dir_list = group[['norm_dec', 'norm_inc']].values.tolist()
        
        # Calculate Fisher mean using PmagPy
        fisher_mean = pmag.fisher_mean(dir_list)
        
        mean_dec = fisher_mean['dec']
        mean_inc = fisher_mean['inc']
        alpha95 = fisher_mean['alpha95']
        n_samples = len(dir_list)

        av_age = group['age_av'].mean()
        lat_n = group['y_lat'].mean()
        long_e = group['x_lon'].mean()
        # Calculate net tectonic rotation relative to stable North America (ref_dec = 0.0)
        rotation_angle = (mean_dec - 0.0 + 180) % 360 - 180
        rotation_rate = rotation_angle / av_age
        
        site_results.append({
            'Location': location,
            'N_Entries': n_samples,
            'lat': lat_n,
            'long': long_e,
            'Mean_Declination': round(mean_dec, 3),
            'Mean_Inclination': round(mean_inc, 3),
            'Av Age': round(av_age, 3),
            'Alpha_95': round(alpha95, 3),
            'Net_Rot_Angle': round(rotation_angle, 3),
            'Rot_Rate': round(rotation_rate, 3)
        })
        
    df_summary = pd.DataFrame(site_results)
    print(df_summary.to_string(index=False))

    # Save df to csv file
    if not header_sent:
        df_summary.to_csv(cleaned_csv, index=False)
        header_sent = True
    else:
        df_summary.to_csv(cleaned_csv, mode='a', header=False, index=False) 

# Explicitly define a file geodatabase in your script's directory so it always has a valid location
script_dir = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
gdb_name = "PaleomagData.gdb"
gdb_path = os.path.join(script_dir, gdb_name)

# Create the file geodatabase if it doesn't already exist
if not arcpy.Exists(gdb_path):
    arcpy.management.CreateFileGDB(script_dir, gdb_name)
    print(f"Created new file geodatabase at: {gdb_path}")

output_feature_class = os.path.join(gdb_path, "Processed_Site_Rotations")


print(f"Successfully processed {len(df)} records. Saved to: {cleaned_csv}")

# Run XYTableToPoint with the fully qualified path
arcpy.env.overwriteOutput = True
try:
    arcpy.management.XYTableToPoint(
        in_table=cleaned_csv,
        out_feature_class=output_feature_class,
        x_field= 'long',
        y_field= 'lat',
        coordinate_system=arcpy.SpatialReference(4326) # WGS 1984
    )
    print(f"\nSuccess! Feature class created inside geodatabase:")
    print(output_feature_class)
except Exception as e:
    print(f"\nArcPy Error encountered: {e}")
