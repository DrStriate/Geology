import os
import pandas as pd
import pmagpy.pmag as pmag
import arcpy
from arcgis.features import GeoAccessor

# Python.exe needs to be in ArcGIS:  C:\Program Files\ArcGIS\Pro\bin\Python\envs\arcgispro-py3\python.exe

file_name = "Magill_1986_Table_1_PO.csv"
location = "PO"

file_path = "PaleomagDataSources/" + file_name
df = pd.read_csv(file_path)
print (f"\n{df}")

# Clean numerical fields
df['lat'] = pd.to_numeric(df['lat'], errors='coerce')
df['long'] = pd.to_numeric(df['long'], errors='coerce')
df['Ma'] = pd.to_numeric(df['Ma'], errors='coerce')
df['norm_dec'] = pd.to_numeric(df['norm_dec'], errors='coerce')
df['norm_inc'] = pd.to_numeric(df['norm_inc'], errors='coerce')

# Drop any rows with dir_dec or dir_inc not set
df = df.dropna(subset=['norm_dec', 'norm_inc'])

# Convert MagIC 0-360 longitude to standard -180 to 180 for ArcGIS
df['long'] = df['long'].apply(lambda lon: lon - 360 if lon > 180 else lon)

# Handle polarity reversals (flip negative inclinations and their declinations by 180°)
def normalize_polarity(row):
    dec = row['norm_dec']
    inc = row['norm_inc']
    if inc < 0:
        return (dec + 180) % 360, -inc
    return dec, inc

df[['norm_dec', 'norm_inc']] = df.apply(
    lambda row: pd.Series(normalize_polarity(row)), axis=1
)

site_results = []
for location, group in df.groupby('location'):

    dir_list = group[['norm_dec', 'norm_inc']].values.tolist()

    # Calculate Fisher mean using PmagPy
    fisher_mean = pmag.fisher_mean(dir_list)
    if 'alpha95' not in fisher_mean:
        print(f"Can't get Fisher mean stats in {file_path}")
        
    mean_dec = fisher_mean['dec']
    mean_inc = fisher_mean['inc']
    alpha95 = fisher_mean['alpha95']
    n_samples = len(dir_list)

    av_age = group['Ma'].mean()
    lat_n = group['lat'].mean()
    long_e = group['long'].mean()
    # Calculate net tectonic rotation relative to stable North America (ref_dec = 0.0)
    rotation_angle = (mean_dec - 0.0 + 180) % 360 - 180
    rotation_rate = rotation_angle / av_age

    site_entry = {
        'file_name': file_name,
        'location': location,
        'N': n_samples,
        'lat': lat_n,
        'long': long_e,
        'D': round(mean_dec, 3),
        'I': round(mean_inc, 3),
        'ma_av': round(av_age, 3),
        'alpha_95': round(alpha95, 3),
        'rot_Angle': round(rotation_angle, 3),
        'rot_Rate': round(rotation_rate, 3)
    }
    site_results.append(site_entry)

df_summary = pd.DataFrame(site_results)

# Explicitly define a file geodatabase in your script's directory so it always has a valid location
script_dir = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
gdb_name = "PaleomagData.gdb"
gdb_path = os.path.join(script_dir, gdb_name)

# Create the file geodatabase if it doesn't already exist
if not arcpy.Exists(gdb_path):
    arcpy.management.CreateFileGDB(script_dir, gdb_name)
    print(f"Created new file geodatabase at: {gdb_path}")

output_feature_class = os.path.join(gdb_path, "Processed_Paleo_Rotations")

# Convert DataFrame into a Spatially Enabled DataFrame (SEDF)
sedf = GeoAccessor.from_xy(
    df_summary, # df
    x_column="long",  
    y_column="lat",  
    sr=4326 # WGS 1984
)

# Export directly to a File Geodatabase Feature Class
sedf.spatial.to_featureclass(location=output_feature_class)

