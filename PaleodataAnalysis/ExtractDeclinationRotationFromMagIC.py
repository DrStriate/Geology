import pandas as pd
import numpy as np

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

# 1. Load your MagIC file
magic_data = parse_magic_file('magic_contribution_18104.txt')
df = magic_data['locations']

# 2. Convert key columns to numeric types
df['lat_n'] = pd.to_numeric(df['lat_n'], errors='coerce')
df['lon_w'] = pd.to_numeric(df['lon_w'], errors='coerce')
df['age_high'] = pd.to_numeric(df['age_high'], errors='coerce')
df['age_low'] = pd.to_numeric(df['age_low'], errors='coerce')
df['dir_dec'] = pd.to_numeric(df['dir_dec'], errors='coerce')
df['dir_inc'] = pd.to_numeric(df['dir_inc'], errors='coerce')
df['dir_tilt_correction'] = pd.to_numeric(df['dir_tilt_correction'], errors='coerce')

# Convert MagIC 0-360 longitude to standard -180 to 180 for ArcGIS
df['x_lon'] = df['lon_w'].apply(lambda lon: lon - 360 if lon > 180 else lon)
df['y_lat'] = df['lat_s']
df['age_av'] = (df['age_high'] + df['age_low']) / 2.0


# 3. Filter for fully tilt-corrected data (100% structural restoration)
df_corrected = df[df['dir_tilt_correction'] == 100].copy()

# 4. Normalize declination for polarity and compute rotation angle
def compute_tectonic_rotation(row, ref_dec=0.0):
    dec = row['dir_dec']
    inc = row['dir_inc']
    
    # If reversed polarity (negative inclination), flip 180 degrees to common frame
    if inc < 0:
        normalized_dec = (dec + 180) % 360
    else:
        normalized_dec = dec
        
    # Calculate net clockwise/counter-clockwise rotation relative to reference
    rotation_angle = (normalized_dec - ref_dec + 180) % 360 - 180
    return normalized_dec, rotation_angle

df_corrected[['normalized_dec', 'rotation_angle']] = df_corrected.apply(
    lambda row: pd.Series(compute_tectonic_rotation(row)), axis=1
)

# 5. Display summary results
results = df_corrected[['location', 'x_lon', 'y_lat', 'dir_dec', 'dir_inc', 'normalized_dec', 'rotation_angle', 'age_av']]
print(results.to_string(index=False))