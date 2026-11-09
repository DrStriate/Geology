import pandas as pd
import numpy as np
import pmagpy.pmag as pmag

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

# 1. Load the MagIC file
magic_data = parse_magic_file('magic_contribution_18469.txt')
df = magic_data['locations']

# 2. Clean numerical fields
df['dir_dec'] = pd.to_numeric(df['dir_dec'], errors='coerce')
df['dir_inc'] = pd.to_numeric(df['dir_inc'], errors='coerce')
df['age_av'] = (pd.to_numeric(df['age_low'], errors='coerce') + pd.to_numeric(df['age_high'], errors='coerce')) / 2.0
df['dir_tilt_correction'] = pd.to_numeric(df['dir_tilt_correction'], errors='coerce')

# 3. Filter for fully tilt-corrected data (100% structural restoration)[cite: 1, 2]
df_corrected = df[df['dir_tilt_correction'] == 100].copy()

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
site_results = []
for location, group in df_corrected.groupby('location'):
    # Prepare list of [declination, inclination] pairs for PmagPy
    dir_list = group[['norm_dec', 'norm_inc']].values.tolist()
    
    # Calculate Fisher mean using PmagPy
    fisher_mean = pmag.fisher_mean(dir_list)
    
    mean_dec = fisher_mean['dec']
    mean_inc = fisher_mean['inc']
    alpha95 = fisher_mean['alpha95']
    n_samples = len(dir_list)
    
    # Calculate net tectonic rotation relative to stable North America (ref_dec = 0.0)
    rotation_angle = (mean_dec - 0.0 + 180) % 360 - 180
    
    site_results.append({
        'Location': location,
        'N_Entries': n_samples,
        'Mean_Declination': round(mean_dec, 2),
        'Mean_Inclination': round(mean_inc, 2),
        #'Av Age': df_corrected.age_av,
        'Alpha_95': round(alpha95, 2),
        'Net_Rotation_Angle': round(rotation_angle, 2)
    })

df_summary = pd.DataFrame(site_results)
print(df_summary.to_string(index=False))